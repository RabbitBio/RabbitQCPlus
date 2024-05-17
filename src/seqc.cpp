//
// Created by ylf9811 on 2021/7/6.
//

#include "seqc.h"

using namespace std;


extern "C" {
#include <athread.h>
#include <pthread.h>
    void slave_tgsfunc();
    void slave_ngsfunc();
    void slave_formatfunc();
    void slave_SeFormatQC();
    void slave_semergefunc();
    void slave_seallfunc();
    void slave_decompressfunc();
    void slave_compressfunc();
}


void SkipToLineEnd(char *data_, int64_t &pos_, const int64_t size_) {
    // cerr << "pos: " << pos_ << " size: " << size_ << endl;
    ASSERT(pos_ < size_);

    while (data_[pos_] != '\n' && data_[pos_] != '\r' && pos_ < size_) ++pos_;

    if (data_[pos_] == '\r' && pos_ < size_) {
        if (data_[pos_ + 1] == '\n') {
            //usesCrlf = true;
            ++pos_;
        }
    }
}


int64_t GetNextFastq(char *data_, int64_t pos_, const int64_t size_) {
    SkipToLineEnd(data_, pos_, size_);
    ++pos_;

    // find beginning of the next record
    while (data_[pos_] != '@') {
        SkipToLineEnd(data_, pos_, size_);
        ++pos_;
    }
    int64_t pos0 = pos_;

    SkipToLineEnd(data_, pos_, size_);
    ++pos_;

    if (data_[pos_] == '@')// previous one was a quality field
        return pos_;
    //-----[haoz:] is the following code necessary??-------------//
    SkipToLineEnd(data_, pos_, size_);
    ++pos_;
    if (data_[pos_] != '+') cout << "core dump is pos: " << pos_ << " char: " << data_[pos_] << endl;
    ASSERT(data_[pos_] == '+');// pos0 was the start of tag
    return pos0;
}

static bool use_in_mem = 1;
static bool use_out_mem = 0;

static bool use_swidx_file = 0;

#define use_mpi_file

char* OutMemData;


struct format_data2 {
    uint64_t *res;
    int my_rank;
};


/**
 * @brief Construct function
 * @param cmd_info1 : cmd information
 */
SeQc::SeQc(CmdInfo *cmd_info1, int my_rank_, int comm_size_) {
    my_rank = my_rank_;
    comm_size = comm_size_;

    part_sizes = new int64_t[comm_size];
    now_pos_ = 0;
    zip_now_pos_ = 0;
    cmd_info_ = cmd_info1;
    filter_ = new Filter(cmd_info1);
    done_thread_number_ = 0;
    int out_block_nums = int(1.0 * cmd_info1->in_file_size1_ / cmd_info1->out_block_size_);
    out_queue_ = NULL;
    p_out_queue_ = NULL;
    in_is_zip_ = cmd_info1->in_file_name1_.find(".gz") != string::npos;
    out_is_zip_ = cmd_info1->out_file_name1_.find(".gz") != string::npos;

    if(in_is_zip_ && comm_size > 1) {
#ifdef USE_LIBDEFLATE
#else
        if(my_rank == 0) fprintf(stderr, "compress file input for multi-process is not support, please use -DUSE_LIBDEFLATE in makefile.\n");
        exit(0);
#endif
    }
    ifstream gFile;
    gFile.open(cmd_info1->in_file_name1_.c_str());
    gFile.seekg(0, ios_base::end);
    long long real_file_size = gFile.tellg();
    gFile.close();
    int64_t start_pos, end_pos;
    if(in_is_zip_) {
#ifdef USE_CC_GZ
        for(int i = 0; i < 64; i++) {
            cc_gz_in_buffer[i] = new char[BLOCK_SIZE];
        }
#endif
        vector<size_t> block_sizes;
        if(use_swidx_file) {
            ifstream iff_idx;
            iff_idx.open(cmd_info1->in_file_name1_ + ".swidx");
            size_t block_size = 0;
            while (iff_idx >> block_size) {
                block_sizes.push_back(block_size);
            }
            iff_idx.close();
        } else {
            ifstream iff_idx;
            iff_idx.open(cmd_info1->in_file_name1_, ios::binary | ios::ate);
            if (!iff_idx.is_open()) {
                cerr << "Failed to open file!" << endl;
                exit(0);
            }
            size_t file_size = iff_idx.tellg();
            size_t blocknum = 0;
            iff_idx.seekg(-static_cast<int>(sizeof(size_t)), ios::end);
            iff_idx.read(reinterpret_cast<char*>(&blocknum), sizeof(size_t));

            real_file_size -= (blocknum + 1) * sizeof(size_t);
            block_sizes.reserve(blocknum);
            iff_idx.seekg(-static_cast<int>((blocknum + 1) * sizeof(size_t)), ios::end);

            size_t block_size = 0;
            for (size_t i = 0; i < blocknum; ++i) {
                iff_idx.read(reinterpret_cast<char*>(&block_size), sizeof(size_t));
                block_sizes.push_back(block_size);
            }
            iff_idx.close();
            
        }
        if (block_sizes.empty()) {
            // Handle error or empty file scenario
            cerr << "Error: No data read from index file." << endl;
            exit(0);
        }
        int pre_size = ceil(1.0 * block_sizes.size() / comm_size);
        int start_line = pre_size * my_rank;
        int end_line = min(start_line + pre_size, (int) block_sizes.size());

        for (int i = 1; i < block_sizes.size(); i++) {
            block_sizes[i] += block_sizes[i - 1];
        }
        start_pos = (start_line == 0) ? 0 : block_sizes[start_line - 1];
        end_pos = (end_line == block_sizes.size()) ? real_file_size : block_sizes[end_line - 1];
        assert(block_sizes.back() == real_file_size);
        start_line_ = start_line;
        end_line_ = end_line;
    } else {
        int64_t pre_size = ceil(1.0 * real_file_size / comm_size);
        start_pos = pre_size * my_rank;
        end_pos = start_pos + pre_size;
        if(end_pos > real_file_size) end_pos = real_file_size;
    }
    if(in_is_zip_) {
        fprintf(stderr, "rank%d line: [%d %d]\n", my_rank, start_line_, end_line_);
    }

    int64_t right_pos;
    if(in_is_zip_) {
        right_pos = 0;
    } else {
        FILE *pre_fp;
        pre_fp = fopen(cmd_info1->in_file_name1_.c_str(), "rb");
        fseek(pre_fp, start_pos, SEEK_SET);
        char *tmp_chunk = new char[1 << 20];
        int res_size = 0;
        res_size = fread(tmp_chunk, sizeof(char), 1 << 20, pre_fp);
        MPI_Barrier(MPI_COMM_WORLD);
        //cerr << "res_size" << my_rank << " " << res_size << endl;
        if(my_rank == 0) right_pos = 0;
        else right_pos = GetNextFastq(tmp_chunk, 0, res_size);
        //cerr << "right_pos" << my_rank << " " << right_pos << endl;
        fclose(pre_fp);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    long long now_pos = right_pos + start_pos;
    long long now_poss[comm_size];
    now_poss[0] = now_pos;
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank) {
        MPI_Send(&now_pos, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
    } else {
        for(int ii = 1; ii < comm_size; ii++) {
            MPI_Recv(&(now_poss[ii]), 1, MPI_LONG, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } 
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0) {
        for(int ii = 1; ii < comm_size; ii++) {
            MPI_Send(now_poss, comm_size, MPI_LONG, ii, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(now_poss, comm_size, MPI_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for(int i = 0; i < comm_size; i++) {
        if(i == comm_size - 1) part_sizes[i] = real_file_size - now_poss[i];
        else part_sizes[i] = now_poss[i + 1] - now_poss[i];
    }

    start_pos_ = now_poss[my_rank];
    if(my_rank == 0 && start_pos_ != 0) fprintf(stderr, "GG size division\n");
    if(my_rank == comm_size - 1) end_pos_ = real_file_size;
    else end_pos_ = now_poss[my_rank + 1];



    p_out_queue_ = new vector<rabbit::fq::FastqDataChunk *>[1 << 20];
    p_queueP1 = 0;
    p_queueP2 = 0;
    p_queueNumNow = 0;
    p_queueSizeLim = Q_lim_se;

     
    if (cmd_info1->write_data_) {
        out_queue_ = new CIPair[1 << 20];
//        out_queue_ = new QChunkItem[1 << 20];
        queueP1 = 0;
        queueP2 = 0;
        queueNumNow = 0;
        queueSizeLim = Q_lim_se * 64;
//        queueSizeLim = Q_lim_se;
        if(use_out_mem) OutMemData = new char[128 << 20];
        if (out_is_zip_) {
            if (cmd_info1->use_pigz_) {
                fprintf(stderr, "use pigz TODO...\n");
                exit(0);
#ifdef Verbose
                printf("now use pigz to compress output data\n");
#endif
            } else {
#ifdef Verbose
                printf("open gzip stream %s\n", cmd_info1->out_file_name1_.c_str());
#endif

#ifdef USE_LIBDEFLATE
                //string idx_name = cmd_info1->out_file_name1_ + ".swidx";
                //fprintf(stderr, "idx_name %s\n", idx_name.c_str());
                //off_idx.open(idx_name);

                ofstream file(cmd_info1->out_file_name1_.c_str());

                if (file.is_open()) {
                    file << "Hello, file!" << endl;
                    file.close();
                } else {
                    cerr << "Failed to open file: " << cmd_info1->out_file_name1_ << endl;
                }

                int err = MPI_File_open(MPI_COMM_WORLD, cmd_info1->out_file_name1_.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

                if (err != MPI_SUCCESS) {
                    char error_string[BUFSIZ];
                    int length_of_error_string, error_class;

                    MPI_Error_class(err, &error_class);
                    MPI_Error_string(error_class, error_string, &length_of_error_string);

                    fprintf(stderr, "%3d: Couldn't open file, MPI Error: %s\n", my_rank, error_string);

                    MPI_Abort(MPI_COMM_WORLD, err);
                    exit(0);
                }
#else 
                zip_out_stream = gzopen(cmd_info1->out_file_name1_.c_str(), "w");
                gzsetparams(zip_out_stream, cmd_info1->compression_level_, Z_DEFAULT_STRATEGY);
                gzbuffer(zip_out_stream, 1024 * 1024);
#endif
            }
        } else {
#ifdef Verbose
            printf("open stream %s\n", cmd_info1->out_file_name1_.c_str());
#endif

 

#ifdef use_mpi_file

            ofstream file(cmd_info1->out_file_name1_.c_str());

            if (file.is_open()) {
                file << "Hello, file!" << endl;
                file.close();
                //cout << "File " << cmd_info1->out_file_name1_ << " has been created and written successfully." << endl;
            } else {
                //cerr << "Failed to open file: " << cmd_info1->out_file_name1_ << endl;
            }

            int err = MPI_File_open(MPI_COMM_WORLD, cmd_info1->out_file_name1_.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

            if (err != MPI_SUCCESS) {
                char error_string[BUFSIZ];
                int length_of_error_string, error_class;

                MPI_Error_class(err, &error_class);
                MPI_Error_string(error_class, error_string, &length_of_error_string);

                fprintf(stderr, "%3d: Couldn't open file, MPI Error: %s\n", my_rank, error_string);

                MPI_Abort(MPI_COMM_WORLD, err);
                exit(0);
            }


#else
            if(my_rank == 0) {
                printf("pre real file size %lld\n", real_file_size);
                int fd = open(cmd_info1->out_file_name1_.c_str(), O_CREAT | O_TRUNC | O_RDWR | O_EXCL, 0644);
                ftruncate(fd, sizeof(char) * real_file_size);
                out_stream_ = fdopen(fd, "w");
            } else {
                int found = 0;
                do {
                    if(-1 == access(cmd_info1->out_file_name1_.c_str(), F_OK)) {
                        if(ENOENT == errno) {
                            fprintf(stderr, "rank%d waiting file1\n");
                            //cerr << "waiting file be created..." << endl;
                            usleep(1000000);
                        } else {
                            fprintf(stderr, "rank%d waiting file2\n");
                            //cerr << "file open GG" << endl;
                            usleep(1000000);
                            //exit(0);
                        }
                    } else {
                        out_stream_ = fopen(cmd_info1->out_file_name1_.c_str(), "r+b");   
                        found = 1;
                    }
                } while(found == 0);
            }
#endif
            MPI_Barrier(MPI_COMM_WORLD);
            //if((out_stream_ = fopen(cmd_info1->out_file_name1_.c_str(),"w")) == NULL) {
            //    printf("File cannot be opened\n");
            //}
            //out_stream_.open(cmd_info1->out_file_name1_);
        }
    }

    //for(int i = 0; i < 64; i++) {
    //    duplicate_[i] = NULL;
    //}
    duplicate_ = NULL;
    if (cmd_info1->state_duplicate_) {
        //for(int i = 0; i < 64; i++) {
        //    duplicate_[i] = new Duplicate(cmd_info1);
        //}
        duplicate_ = new Duplicate(cmd_info1);
    }
    umier_ = NULL;
    if (cmd_info1->add_umi_) {
        umier_ = new Umier(cmd_info1);
    }
    if (cmd_info1->use_pugz_) {
        fprintf(stderr, "use pugz TODO...\n");
        exit(0);
    }
    if (cmd_info1->use_pigz_) {
        fprintf(stderr, "use pigz TODO...\n");
        exit(0);
    }
    pugzDone = 0;
//    producerDone = 0;
    writerDone = 0;
//    consumerCommDone = 0;
//    producerStop = 0;
//    now_chunks = 0;
//    mx_chunks = 0;
}

SeQc::~SeQc() {

    delete[] part_sizes;
    delete filter_;
    delete[] p_out_queue_;
    if (cmd_info_->write_data_) {
        delete[] out_queue_;
    }
    if(in_is_zip_) {
#ifdef USE_CC_GZ
        for(int i = 0; i < 64; i++) {
            delete []cc_gz_in_buffer[i];
        }
#endif
    }

    if (cmd_info_->state_duplicate_) {
        //for(int i = 0; i <64; i++) {
        //    delete duplicate_[i];
        //}
        delete duplicate_;
    }
    if (cmd_info_->add_umi_) {
        delete umier_;
    }
}


/**
 * @brief get fastq data chunk from fastq_data_pool and put it into data queue
 * @param file : fastq file name, which is also the input file
 * @param fastq_data_pool : fastq data pool
 * @param dq : a data queue
 */

static rabbit::fq::FastqFileReader *fqFileReader;

void SeQc::ProducerSeFastqTask64(string file, rabbit::fq::FastqDataPool *fastq_data_pool) {



    double t0 = GetTime();
    rabbit::uint32 tmpSize = SWAP1_SIZE;
    if (cmd_info_->seq_len_ <= 200) tmpSize = SWAP1_SIZE;
    int64_t n_chunks = 0;

    double t_sum1 = 0;
    double t_sum2 = 0;
    double t_sum2_1 = 0;
    double t_sum2_2 = 0;
    double t_sum2_3 = 0;
    double t_sum3 = 0;
    bool is_first = 1;
    int64_t tot_size = 0;
    vector<rabbit::fq::FastqDataChunk *> tmp_chunks;
    while (true) {
        double tt0 = GetTime();
        rabbit::fq::FastqDataChunk *fqdatachunk;
        int64_t offset;
        if(is_first) {
            offset = 0;
            for(int i = 0; i < my_rank; i++)
                offset += part_sizes[i];
        } else {
            offset = -1;
        }
        is_first = 0;
        fqdatachunk = fqFileReader->readNextChunk();
        t_sum1 += GetTime() - tt0;

        tt0 = GetTime();
        //if(fqdatachunk != NULL) {
        //    fprintf(stderr, "ppp chunk is not null\n");
        //    fprintf(stderr, "ppp chunk size %llu\n", fqdatachunk->size);
        //}
        if ((fqdatachunk == NULL) || (fqdatachunk != NULL && fqdatachunk->size == 1ll << 32)) {
            double tt1 = GetTime();
            if(tmp_chunks.size()) {
                //fprintf(stderr, "ppp put last %d\n", tmp_chunks.size());
                //p_mylock.lock();
                while (p_queueNumNow >= p_queueSizeLim) {
                    usleep(100);
                }
                p_out_queue_[p_queueP2++] = tmp_chunks;
                p_queueNumNow++;
                //p_mylock.unlock();
                tmp_chunks.clear();
                n_chunks++;
            }
            t_sum2_1 += GetTime() - tt1;
//            producerStop = 1;
            //fprintf(stderr, "ppp producer rank%d stop, tot chunk size %d\n", my_rank, n_chunks);
            //if(cmd_info_->write_data_) {
            tmp_chunks.clear();
            while (p_queueNumNow >= p_queueSizeLim) {
                usleep(100);
            }
            p_out_queue_[p_queueP2++] = tmp_chunks;
            p_queueNumNow++;
            //} 
            break;
        }
        t_sum2 += GetTime() - tt0;

        tt0 = GetTime();
        tot_size += fqdatachunk->size;
        tmp_chunks.push_back(fqdatachunk); 
        if(tmp_chunks.size() == 64) {
            //fprintf(stderr, "ppp put 64\n");
            //p_mylock.lock();
            while (p_queueNumNow >= p_queueSizeLim) {
                usleep(100);
            }
            p_out_queue_[p_queueP2++] = tmp_chunks;
            p_queueNumNow++;
            //p_mylock.unlock();
            tmp_chunks.clear();
            n_chunks++;
        }
        t_sum3 += GetTime() - tt0;
    }
    //printf("totsize %lld\n", tot_size);

    fprintf(stderr, "producer%d sum1 cost %lf\n", my_rank, t_sum1);
    fprintf(stderr, "producer%d sum2 cost %lf [%lf %lf %lf]\n", my_rank, t_sum2, t_sum2_1, t_sum2_2, t_sum2_3);
    fprintf(stderr, "producer%d sum3 cost %lf\n", my_rank, t_sum3);
    fprintf(stderr, "producer %d cost %lf\n", my_rank, GetTime() - t0);
    fprintf(stderr, "producer %d done in func\n", my_rank);
    //dq.SetCompleted();
//    producerDone = 1;
}

string SeQc::Read2String(neoReference &ref) {
    return string((char *) ref.base + ref.pname, ref.lname) + "\n" +
        string((char *) ref.base + ref.pseq, ref.lseq) + "\n" +
        string((char *) ref.base + ref.pstrand, ref.lstrand) + "\n" +
        string((char *) ref.base + ref.pqual, ref.lqual) + "\n";
}

void SeQc::Read2Chars(neoReference &ref, char *out_data, int &pos) {
    memcpy(out_data + pos, ref.base + ref.pname, ref.lname);
    pos += ref.lname;
    out_data[pos++] = '\n';
    memcpy(out_data + pos, ref.base + ref.pseq, ref.lseq);
    pos += ref.lseq;
    out_data[pos++] = '\n';
    memcpy(out_data + pos, ref.base + ref.pstrand, ref.lstrand);
    pos += ref.lstrand;
    out_data[pos++] = '\n';
    memcpy(out_data + pos, ref.base + ref.pqual, ref.lqual);
    pos += ref.lqual;
    out_data[pos++] = '\n';
}

struct format_data {
    vector <neoReference> *data;
    rabbit::fq::FastqDataChunk *fqdatachunk;
    uint64_t *res;
    int fqSize;
    int my_rank;
};

struct SeMerge_data {
    char* out_data[64];
    int out_lens[64] = {0};
    vector <neoReference> *pass_data[64];
    int pass_data_size[64];
};

struct SeAll_data{
    format_data para1[64];
    qc_data *para2;
    SeMerge_data *para3;
    int out_len_slave[64];
    bool write_data;
};


double tsum1 = 0;
double tsum2 = 0;
double tsum3 = 0;
double tsum4 = 0;
double tsum4_1 = 0;
double tsum4_1_1 = 0;
double tsum4_1_2 = 0;
double tsum4_1_3 = 0;
double tsum4_2 = 0;
double tsum4_2_1 = 0;
double tsum4_2_2 = 0;
double tsum4_2_3 = 0;
double tsum4_3 = 0;
double tsum4_3_1 = 0;
double tsum4_3_2 = 0;
double tsum4_3_3 = 0;
double tsum4_3_4 = 0;
double tsum4_3_5 = 0;
double tsum4_4 = 0;
double tsum5 = 0;
double tsum6 = 0;
double tsum7 = 0;

double t_wait_producer = 0;
double t_format = 0;
double t_decom = 0;
double t_resize = 0;
double t_ngsfunc = 0;
double t_slave_gz = 0;
double t_slave_gz1 = 0;
double t_slave_gz2 = 0;
double t_slave_gz3_1 = 0;
double t_slave_gz3_2 = 0;
double t_slave_gz3_3 = 0;
double t_slave_gz4 = 0;
double t_slave_gz5 = 0;
double t_push_q = 0;
double t_copy_data = 0;

void PrintRef(neoReference &ref) {
    printf("%.*s\n", ref.lname, (char *) ref.base + ref.pname);
    printf("%.*s\n", ref.lseq, (char *) ref.base + ref.pseq);
    printf("%.*s\n", ref.lstrand, (char *) ref.base + ref.pstrand);
    printf("%.*s\n", ref.lqual, (char *) ref.base + ref.pqual);
}

struct ConsumerWriteInfo{
    char* buffer;
    int buffer_len;
    long long file_offset;
};

vector<ConsumerWriteInfo> writeInfos;

static int se_pre_out_len_slave[64];
static int se_pre_pass_data_size[64];

void SeQc::ProcessFormatQCWrite(bool &allIsNull, vector <neoReference> *data, vector <neoReference> *pass_data, vector <neoReference> *pre_pass_data,
                           vector <rabbit::fq::FastqDataChunk *> fqdatachunks, vector <rabbit::fq::FastqDataChunk *> pre_fqdatachunks, qc_data *para,
                           rabbit::fq::FastqDataPool *fastq_data_pool) {

    double tt0 = GetTime();
    format_data para2[64];
    uint64_t counts[64] = {0};
    int fq_size = fqdatachunks.size();
    for(int i = 0; i < 64; i++) {
        //data[i].clear();
        para2[i].fqSize = fq_size;
        para2[i].my_rank = my_rank;
    }
    for(int i = 0; i < fq_size; i++) {
        //if(fqdatachunks[i] != NULL) {
        //    data[i].resize(fqdatachunks[i]->size / (cmd_info_->seq_len_ * 2));
        //    if(cmd_info_->write_data_) {
        //        pass_data[i].resize(data[i].size());
        //    }
        //} else {
        //    data[i].resize(0);
        //    if(cmd_info_->write_data_) {
        //        pass_data[i].resize(0);
        //    }
        //}
        para2[i].data = &(data[i]);
        para2[i].fqdatachunk = fqdatachunks[i];
        para2[i].res = counts;
        para2[i].fqSize = fq_size;
    }
    for(int i = 0; i < 64; i++) {
        para->data1_[i] = &data[i];
        para->pass_data1_[i] = &pass_data[i];
    }
    tsum1 += GetTime() - tt0;


    tt0 = GetTime();
    SeMerge_data se_merge_data;
    SeAll_data se_all_data;
    se_all_data.write_data = cmd_info_->write_data_;
    int out_lens[64] = {0};
    char *out_data[64];
    for(int i = 0; i < 64; i++) {
        se_all_data.para1[i] = para2[i];
        se_all_data.out_len_slave[i] = 0;
        if (cmd_info_->write_data_) {
            out_lens[i] = se_pre_out_len_slave[i];
            if(out_lens[i]) out_data[i] = new char[out_lens[i]];
            else out_data[i] = (char *)(NULL);
        }
    }
    for(int i = 0; i < 64; i++) {
        se_merge_data.out_data[i] = out_data[i];
        se_merge_data.out_lens[i] = out_lens[i];
        se_merge_data.pass_data[i] = &(pre_pass_data[i]);
        se_merge_data.pass_data_size[i] = se_pre_pass_data_size[i];
    }
    se_all_data.para2 = para;
    se_all_data.para3 = &se_merge_data;
    tsum2 += GetTime() - tt0;

    tt0 = GetTime();
    __real_athread_spawn((void *)slave_seallfunc, &se_all_data, 1);
    athread_join();
    tsum3 += GetTime() - tt0;


    tt0 = GetTime();
    if (cmd_info_->write_data_) {
        long long now_pos_base = 0;
        bool use_consumer_libdeflate = 0;

#ifdef CONSUMER_USE_LIBDEFLATE
        use_consumer_libdeflate = 1;
#endif
        if(!out_is_zip_ || use_consumer_libdeflate == 0) {
            int out_len = 0;
            for(int i = 0; i < 64; i++) {
                out_len += out_lens[i];
            }
            int now_size = out_len;
            int now_sizes[comm_size];
            now_sizes[0] = now_size;

            MPI_Barrier(MPI_COMM_WORLD);

            if(my_rank) {
                MPI_Send(&now_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            } else {
                for(int ii = 1; ii < comm_size; ii++) {
                    MPI_Recv(&(now_sizes[ii]), 1, MPI_INT, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            if(my_rank == 0) {
                for(int ii = 1; ii < comm_size; ii++) {
                    MPI_Send(now_sizes, comm_size, MPI_INT, ii, 0, MPI_COMM_WORLD);
                }
            } else {
                MPI_Recv(now_sizes, comm_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            int pre_sizes[comm_size];
            pre_sizes[0] = 0;
            for(int ii = 1; ii < comm_size; ii++) {
                pre_sizes[ii] = pre_sizes[ii - 1] + now_sizes[ii - 1];
            }
            //pre_sizes[comm_size - 1] = 0;
            //for(int ii = comm_size - 2; ii >= 0; ii--) {
            //    pre_sizes[ii] = pre_sizes[ii + 1] + now_sizes[ii + 1];
            //}
            now_pos_base = now_pos_ + pre_sizes[my_rank];
            MPI_Barrier(MPI_COMM_WORLD);
            for(int ii = 0; ii < comm_size; ii++) {
                now_pos_ += now_sizes[ii];
            }
        }
        writeInfos.clear();
        for(int i = 0; i < 64; i++) {
            writeInfos.push_back({out_data[i], out_lens[i], now_pos_base});
        }
    }
    for(int i = 0; i < 64; i++) {
        se_pre_out_len_slave[i] = se_all_data.out_len_slave[i];
        se_pre_pass_data_size[i] = se_merge_data.pass_data_size[i];
    }
    tsum4 += GetTime() - tt0;


    tt0 = GetTime();
    int isNull = 1;
    int all_stop = 1;
    {
        for(int i = 0; i < 64; i++) {
            //data[i].resize(counts[i]);
            //if(cmd_info_->write_data_) {
            //    pass_data[i].resize(data[i].size());
            //}
            //if(data[i].size()) isNull = 0;
            if(counts[i]) isNull = 0;
        }
        int isNulls[comm_size];
        isNulls[0] = isNull;
        MPI_Barrier(MPI_COMM_WORLD);
        if(my_rank) {
            MPI_Send(&isNull, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        } else {
            for(int ii = 1; ii < comm_size; ii++) {
                MPI_Recv(&(isNulls[ii]), 1, MPI_INT, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if(my_rank == 0) {
            for(int ii = 1; ii < comm_size; ii++) {
                MPI_Send(isNulls, comm_size, MPI_INT, ii, 0, MPI_COMM_WORLD);
            }
        } else {
            MPI_Recv(isNulls, comm_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for(int ii = 0; ii < comm_size; ii++) {
            if(isNulls[ii] == 0) all_stop = 0;
        }
        if(all_stop) allIsNull = 1;
    }
    tsum5 += GetTime() - tt0;


    tt0 = GetTime();
    if(!all_stop && isNull) {
        fprintf(stderr, "consumer%d push null\n", my_rank);
        vector<rabbit::fq::FastqDataChunk *> tmp_chunks;
        tmp_chunks.clear();
        while (p_queueNumNow >= p_queueSizeLim) {
            usleep(100);
        }
        p_out_queue_[p_queueP2++] = tmp_chunks;
        p_queueNumNow++;
    }
    tsum6 += GetTime() - tt0;


    tt0 = GetTime();
    for(int i = 0; i < 64; i++) {
        if(pre_fqdatachunks[i] != NULL) fastq_data_pool->Release(pre_fqdatachunks[i]);
    }
    tsum7 += GetTime() - tt0;

}

inline void gather_and_sort_vectors_se(int rank, int comm_size, int out_round, std::vector<std::pair<int, size_t>>& local_vec, int root) {
    int local_size = local_vec.size();
    std::vector<int> send_data(local_size * 2);
    for (int i = 0; i < local_size; ++i) {
        send_data[2 * i] = local_vec[i].first;
        send_data[2 * i + 1] = static_cast<int>(local_vec[i].second);  // Assuming size_t can be safely cast to int
    }

    // Root process initializes a vector to hold all received data
    std::vector<int> recv_data;
    if (rank == root) {
        recv_data.reserve(local_size * 2 * comm_size);  // Reserve enough space
    }

    // Root process receives data from all other processes
    if (rank == root) {
        std::vector<int> buffer(local_size * 2);
        for (int i = 0; i < comm_size; ++i) {
            if (i == root) {
                // Include the root's own data
                recv_data.insert(recv_data.end(), send_data.begin(), send_data.end());
            } else {
                MPI_Recv(buffer.data(), local_size * 2, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                recv_data.insert(recv_data.end(), buffer.begin(), buffer.end());
            }
        }
    } else {
        // Other processes send their data to the root
        MPI_Send(send_data.data(), local_size * 2, MPI_INT, root, 0, MPI_COMM_WORLD);
    }

    // Root process sorts all received data
    if (rank == root) {
        std::vector<std::pair<int, size_t>> all_data;
        for (int i = 0; i < recv_data.size() / 2; ++i) {
            all_data.emplace_back(recv_data[2 * i], static_cast<size_t>(recv_data[2 * i + 1]));
        }
        std::sort(all_data.begin(), all_data.end());

        //// Print sorted data for verification
        //for (const auto& p : all_data) {
        //    std::cout << "(" << p.first << ", " << p.second << ") ";
        //}
        //std::cout << std::endl;

        // Optionally clear and repopulate local_vec
        local_vec.clear();
        for (const auto& item : all_data) {
            local_vec.push_back(item);
        }
    }
}


inline void gather_and_sort_vectors_se2(int rank, int comm_size, int out_round, vector<pair<int, size_t>>& local_vec, int root) {

    int local_size = local_vec.size();
    vector<int> sizes(comm_size, 0);
    MPI_Gather(&local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, root, MPI_COMM_WORLD);
    vector<int> displs(comm_size, 0);
    for (int i = 1; i < comm_size; ++i) {
        displs[i] = displs[i - 1] + sizes[i - 1];
    }
    vector<int> send_data(local_size * 2);
    for (int i = 0; i < local_size; ++i) {
        send_data[2 * i] = local_vec[i].first;
        send_data[2 * i + 1] = local_vec[i].second;
    }
    vector<int> recv_data;
    if (rank == root) {
        recv_data.resize(2 * displs[comm_size - 1] + 2 * sizes[comm_size - 1]);
        for (int i = 0; i < comm_size; ++i) {
            fprintf(stderr, "Size[%d] = %d, Displ[%d] = %d\n", i, sizes[i], i, displs[i]);
        }
        if (rank == root) {
            fprintf(stderr, "Final recv_data size expected: %d, Actual: %lu\n", 2 * displs[comm_size - 1] + 2 * sizes[comm_size - 1], recv_data.size());
        }

    }
    fprintf(stderr, "gather%d == %d %d\n", rank, 2 * local_size, 2 * displs[comm_size - 1] + 2 * sizes[comm_size - 1]);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(send_data.data(), 2 * local_size, MPI_INT,
                recv_data.data(), sizes.data(), displs.data(), MPI_INT, root, MPI_COMM_WORLD);
    if (rank == root) {
        vector<pair<int, size_t>> all_data;
        for (int i = 0; i < recv_data.size() / 2; ++i) {
            all_data.emplace_back(recv_data[2 * i], recv_data[2 * i + 1]);
        }

        // Sort the data
        sort(all_data.begin(), all_data.end());
        local_vec.clear();
        for(auto item : all_data) local_vec.push_back(item);

    }
} 


void SeQc::ConsumerSeFastqTask64(ThreadInfo **thread_infos, rabbit::fq::FastqDataPool *fastq_data_pool) {


    rabbit::int64 id = 0;
    rabbit::fq::FastqDataChunk *fqdatachunk;
    qc_data para;
    para.cmd_info_ = cmd_info_;
    para.thread_info_ = thread_infos;
    para.bit_len = 0;
    para.duplicate_ = duplicate_;

    bool allIsNull = 0;
    if(cmd_info_->state_duplicate_) {
        para.bit_len = duplicate_->key_len_base_;
    }
    athread_init();

    double t0 = GetTime();

    vector<rabbit::fq::FastqDataChunk *> fqdatachunks;
    vector<rabbit::fq::FastqDataChunk *> pre_fqdatachunks;
    int out_round = 0;
    vector<neoReference> data[64];
    vector<neoReference> pass_data[64];
    vector<neoReference> pre_pass_data[64];
    int max_data_size = BLOCK_SIZE / (cmd_info_->seq_len_ * 2);
    for(int i = 0; i < 64; i++) {
        data[i].clear();
        pass_data[i].clear();
        pre_pass_data[i].clear();
        pre_fqdatachunks.push_back(NULL);
        data[i].resize(max_data_size);
        if(cmd_info_->write_data_) {
            pass_data[i].resize(max_data_size);
            pre_pass_data[i].resize(max_data_size);
        }
    }
    while(true) {
        double tt0 = GetTime();
        while (p_queueNumNow == 0) {
            usleep(1000);
        }
        t_wait_producer += GetTime() - tt0;

        tt0 = GetTime();
        fqdatachunks = p_out_queue_[p_queueP1++];
        p_queueNumNow--;
        for(int i = fqdatachunks.size(); i < 64; i++) {
            fqdatachunks.push_back(NULL);
        }
#ifdef USE_CC_GZ
#ifdef USE_LIBDEFLATE
        if(in_is_zip_) {
            Para degz_paras[64];
            size_t out_size[64] = {0};
            for (int i = 0; i < 64; i++) {
                if(fqdatachunks[i] == NULL) {
                    degz_paras[i].in_buffer = NULL;
                    degz_paras[i].in_size = 0;
                    //fprintf(stderr, "in size %d\n", degz_paras[i].in_size);
                    continue;
                }
                memcpy(cc_gz_in_buffer[i], fqdatachunks[i]->data.Pointer(), fqdatachunks[i]->size);
                degz_paras[i].in_buffer = cc_gz_in_buffer[i];
                degz_paras[i].out_buffer = (char*)fqdatachunks[i]->data.Pointer();
                degz_paras[i].in_size = fqdatachunks[i]->size;
                degz_paras[i].out_size = &(out_size[i]);
                //fprintf(stderr, "in size %d\n", degz_paras[i].in_size);
            }

            {
                std::lock_guard <std::mutex> guard(globalMutex);
                __real_athread_spawn((void *) slave_decompressfunc, degz_paras, 1);
                athread_join();
            }

            for (int i = 0; i < 64; i++) {
                if(fqdatachunks[i] == NULL) {
                } else {
                    if(out_size[i]) fqdatachunks[i]->size = out_size[i] - 1;
                    else fqdatachunks[i]->size = out_size[i];
                }
                //fprintf(stderr, "ccc decom size %d\n", out_size[i]);
            }

        }
#endif
#endif
        t_decom += GetTime() - tt0;

        tt0 = GetTime();
        ProcessFormatQCWrite(allIsNull, data, pass_data, pre_pass_data, fqdatachunks, pre_fqdatachunks, &para, fastq_data_pool);
        t_ngsfunc += GetTime() - tt0;

        //fprintf(stderr, "rank%d pending write size %d\n", my_rank, writeInfos.size());
#ifdef CONSUMER_USE_LIBDEFLATE

        tt0 = GetTime();
        if (cmd_info_->write_data_ && out_is_zip_) {
            double tt00 = GetTime();
            Para paras[64];
            size_t out_size[64] = {0};
            for(int i = 0; i < writeInfos.size(); i++) {
                paras[i].in_buffer = writeInfos[i].buffer;
                paras[i].out_buffer = new char[BLOCK_SIZE];
                paras[i].in_size = writeInfos[i].buffer_len;
                paras[i].out_size = &(out_size[i]);
                paras[i].level = 1;
            }
            for(int i = writeInfos.size(); i < 64; i++) {
                paras[i].in_size = 0;
            }
            t_slave_gz1 += GetTime() - tt00;

            tt00 = GetTime();
            {
                lock_guard<mutex> guard(globalMutex);
                __real_athread_spawn((void *)slave_compressfunc, paras, 1);
                athread_join();

            }
            t_slave_gz2 += GetTime() - tt00;

            tt00 = GetTime();
            for(int i = 0; i < 64; i++) {
                out_gz_block_sizes.push_back(make_pair(out_round * comm_size * 64 + my_rank * 64 + i, out_size[i]));
            }
            out_round++;
            int sum_buffer_len = 0;
            t_slave_gz3_1 += GetTime() - tt00;

            tt00 = GetTime();
            for(int i = 0; i < writeInfos.size(); i++) {
                delete[] writeInfos[i].buffer;
            }
            t_slave_gz3_2 += GetTime() - tt00;

            tt00 = GetTime();
            for(int i = 0; i < writeInfos.size(); i++) {
                writeInfos[i].buffer = paras[i].out_buffer;
                writeInfos[i].buffer_len = out_size[i];
                sum_buffer_len += out_size[i];
            }
            int now_size = sum_buffer_len;
            int now_sizes[comm_size];
            now_sizes[0] = now_size;
            t_slave_gz3_3 += GetTime() - tt00;

            tt00 = GetTime();
            MPI_Barrier(MPI_COMM_WORLD);
            t_slave_gz4 += GetTime() - tt00;

            tt00 = GetTime();
            if(my_rank) {
                MPI_Send(&now_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            } else {
                for(int ii = 1; ii < comm_size; ii++) {
                    MPI_Recv(&(now_sizes[ii]), 1, MPI_INT, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            if(my_rank == 0) {
                for(int ii = 1; ii < comm_size; ii++) {
                    MPI_Send(now_sizes, comm_size, MPI_INT, ii, 0, MPI_COMM_WORLD);
                }
            } else {
                MPI_Recv(now_sizes, comm_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            int pre_sizes[comm_size];
            pre_sizes[0] = 0;
            for(int ii = 1; ii < comm_size; ii++) {
                pre_sizes[ii] = pre_sizes[ii - 1] + now_sizes[ii - 1];
            }
            long long now_pos_base = zip_now_pos_ + pre_sizes[my_rank];
            MPI_Barrier(MPI_COMM_WORLD);
            for(int ii = 0; ii < comm_size; ii++) {
                zip_now_pos_ += now_sizes[ii];
            }
            for(int i = 0; i < writeInfos.size(); i++) {
                writeInfos[i].file_offset = now_pos_base;
            }
            t_slave_gz5 += GetTime() - tt00;

        }
        t_slave_gz += GetTime() - tt0;
        
#endif
        
        if (cmd_info_->write_data_) {
            tt0 = GetTime();
            assert(writeInfos.size() == 64);
            long long real_pre_chunk_pos = writeInfos[0].file_offset;
            for(int i = 0; i < 64; i++) {
                //mylock.lock();
                while (queueNumNow >= queueSizeLim) {
#ifdef Verbose
                    //printf("waiting to push a chunk to out queue %d\n",out_len);
#endif
                    usleep(1000);
                }
                //out_queue_[queueP2++] = Qitem;
                out_queue_[queueP2++] = make_pair(writeInfos[i].buffer, make_pair(writeInfos[i].buffer_len, real_pre_chunk_pos));
                queueNumNow++;
                real_pre_chunk_pos += writeInfos[i].buffer_len;
                //mylock.unlock();
            }
            t_push_q += GetTime() - tt0;
        }

        tt0 = GetTime();
        pre_fqdatachunks.clear();
        for(int i = 0; i < 64; i++) {
            if(cmd_info_->write_data_) {
                pass_data[i].swap(pre_pass_data[i]);
            }
            pre_fqdatachunks.push_back(fqdatachunks[i]);
        }
        t_copy_data += GetTime() - tt0;
        if(allIsNull) break;

    }

    double tt0 = GetTime();
    gather_and_sort_vectors_se(my_rank, comm_size, out_round, out_gz_block_sizes, 0);

    printf("comm gz block size cost %lf\n", GetTime() - tt0);



    fprintf(stderr, "consumer%d NGSnew tot cost %lf\n", my_rank, GetTime() - t0);
    fprintf(stderr, "consumer%d NGSnew wait producer cost %lf\n", my_rank, t_wait_producer);
    fprintf(stderr, "consumer%d NGSnew format cost %lf\n", my_rank, t_format);
    fprintf(stderr, "consumer%d NGSnew resize cost %lf\n", my_rank, t_resize);
    fprintf(stderr, "consumer%d NGSnew ngsfunc cost %lf\n", my_rank, t_ngsfunc);
    fprintf(stderr, "consumer%d NGSnew      pre_qc1 cost %lf\n", my_rank, tsum1);
    fprintf(stderr, "consumer%d NGSnew      pre_qc2 cost %lf\n", my_rank, tsum2);
    fprintf(stderr, "consumer%d NGSnew      slave cost %lf\n", my_rank, tsum3);
    fprintf(stderr, "consumer%d NGSnew      write cost %lf (%lf [%lf %lf %lf], %lf [%lf %lf %lf], %lf [%lf %lf %lf %lf %lf], %lf)\n", my_rank, tsum4, tsum4_1, tsum4_1_1, tsum4_1_2, tsum4_1_3, tsum4_2, tsum4_2_1, tsum4_2_2, tsum4_2_3, tsum4_3, tsum4_3_1, tsum4_3_2, tsum4_3_3, tsum4_3_4, tsum4_3_5, tsum4_4);
    fprintf(stderr, "consumer%d NGSnew      isnull comm cost %lf\n", my_rank, tsum5);
    fprintf(stderr, "consumer%d NGSnew      push null cost %lf\n", my_rank, tsum6);
    fprintf(stderr, "consumer%d NGSnew      release cost %lf\n", my_rank, tsum7);
    fprintf(stderr, "consumer%d NGSnew gz slave cost %lf\n", my_rank, t_slave_gz);
    fprintf(stderr, "consumer%d NGSnew gz slave sub cost [%lf %lf %lf %lf %lf %lf %lf]\n", my_rank, t_slave_gz1, t_slave_gz2, t_slave_gz3_1, t_slave_gz3_2, t_slave_gz3_3, t_slave_gz4, t_slave_gz5);
    fprintf(stderr, "consumer%d NGSnew push to queue cost %lf\n", my_rank, t_push_q);
    fprintf(stderr, "consumer%d NGSnew copy data cost %lf\n", my_rank, t_copy_data);
    done_thread_number_++;
    athread_halt();
    fprintf(stderr, "consumer %d done in func\n", my_rank);
}


/**
 * @brief a function to write data from out_data queue to file
 */
void SeQc::WriteSeFastqTask() {
#ifdef Verbose
    double t0 = GetTime();
#endif
    int cnt = 0;
    long long tot_size = 0;
    bool overWhile = 0;
    //QChunkItem Qitem;
    CIPair now;
    double t_wait = 0;
    double t_gz_slave = 0;
    double t_write = 0;
    double t_del = 0;
    double t_free = 0;

    while (true) {
        double tt0 = GetTime();
        while (queueNumNow == 0) {
            if (done_thread_number_ == cmd_info_->thread_number_) {
                overWhile = 1;
                fprintf(stderr, "writer rank %d writer while break\n", my_rank);
                break;
            }
            
            //fprintf(stderr, "writer rank %d write wait\n", my_rank);
            usleep(10000);
        }
        if (overWhile) break;
        //Qitem = out_queue_[queueP1++];
        now = out_queue_[queueP1++];
        queueNumNow--;
        t_wait += GetTime() - tt0;

        //fprintf(stderr, "writer rank %d write get %p\n", my_rank, now.first);

        if (out_is_zip_) {
            if (cmd_info_->use_pigz_) {
                fprintf(stderr, "use pigz TODO...\n");
                exit(0);
            } else {
                tt0 = GetTime();
                if(use_out_mem) {
                    memcpy(OutMemData, now.first, now.second.first);
                } else {
#ifdef USE_LIBDEFLATE

#ifdef use_mpi_file

                    MPI_File_seek(fh, now.second.second, MPI_SEEK_SET);
                    MPI_File_write(fh, now.first, now.second.first, MPI_CHAR, &status);
#else
                    fseek(out_stream_, now.second.second, SEEK_SET);
                    fwrite(now.first, sizeof(char), now.second.first, out_stream_);
#endif


#else
                    int written = gzwrite(zip_out_stream, now.first, now.second.first);
                    if (written != now.second.first) {
                        printf("gzwrite error\n");
                        exit(0);
                    }
#endif
                }

                t_write += GetTime() - tt0;

                tt0 = GetTime();
                delete[] now.first;
                t_del += GetTime() - tt0;
            }
        } else {
            tot_size += now.second.first;
            if(now.second.first) {
                tt0 = GetTime();
                if(use_out_mem) {
                    memcpy(OutMemData, now.first, now.second.first);
                } else {
                    //fprintf(stderr, "rank %d write seek %lld\n", my_rank, now.second.second);

#ifdef use_mpi_file

                    MPI_File_seek(fh, now.second.second, MPI_SEEK_SET);
                    //fprintf(stderr, "rank %d write ww %lld\n", my_rank, now.second.first);
                    MPI_File_write(fh, now.first, now.second.first, MPI_CHAR, &status);
                    //fprintf(stderr, "rank %d write ww done\n", my_rank);
#else
                    fseek(out_stream_, now.second.second, SEEK_SET);
                    //fprintf(stderr, "rank %d write ww %lld\n", my_rank, now.second.first);
                    fwrite(now.first, sizeof(char), now.second.first, out_stream_);
                    //fprintf(stderr, "rank %d write ww done\n", my_rank);
#endif
                }
                t_write += GetTime() - tt0;

                tt0 = GetTime();
                delete[] now.first;
                t_del += GetTime() - tt0;
                //fprintf(stderr, "rank %d write del done\n", my_rank);
            }
            //fprintf(stderr, "writer %d %d done\n", my_rank, cnt++);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    //fprintf(stderr, "writer rank%d write donedone\n", my_rank);
    double tt0 = GetTime();
    if (out_is_zip_) {
        if (cmd_info_->use_pigz_) {
            fprintf(stderr, "use pigz TODO...\n");
            exit(0);
        } else {
#ifdef USE_LIBDEFLATE

#ifdef use_mpi_file
            MPI_File_close(&fh);
#else
            fclose(out_stream_);
            if(my_rank == 0) {
                truncate(cmd_info_->out_file_name1_.c_str(), sizeof(char) * zip_now_pos_);
            }
#endif
            //off_idx.close();

#else
            if (zip_out_stream) {
                gzflush(zip_out_stream, Z_FINISH);
                gzclose(zip_out_stream);
                zip_out_stream = NULL;
            }
#endif
        }
    } else {
#ifdef use_mpi_file
        MPI_File_close(&fh);
#else
        fclose(out_stream_);
        if(my_rank == 0) {
            truncate(cmd_info_->out_file_name1_.c_str(), sizeof(char) * now_pos_);
        }
#endif

    }
    t_free += GetTime() - tt0;

    if(my_rank == 0 && out_is_zip_) {
        tt0 = GetTime();
        ofstream ofs(cmd_info_->out_file_name1_, ios::binary | ios::app);
        for (const auto& pair : out_gz_block_sizes) {
            ofs.write(reinterpret_cast<const char*>(&pair.second), sizeof(size_t));
        }
        size_t vector_size = out_gz_block_sizes.size();
        ofs.write(reinterpret_cast<const char*>(&vector_size), sizeof(size_t));
        ofs.close();
        printf("writer final cost %lf\n", GetTime() - tt0);
    }


#ifdef Verbose
    printf("writer wait queue cost %.4f\n", t_wait);
    printf("writer gz slave cost %.4f\n", t_gz_slave);
    printf("writer write cost %.4f\n", t_write);
    printf("writer del cost %.4f\n", t_del);
    printf("writer free cost %.4f\n", t_free);
    //printf("writer %d cost %lf, tot size %lld\n", my_rank, GetTime() - t0, tot_size);
    fprintf(stderr, "writer %d cost %lf, tot size %lld\n", my_rank, GetTime() - t0, tot_size);
    printf("writer %d done in func\n", my_rank);
    //printf("write %d cost %.5f --- %lld\n", my_rank, GetTime() - t0, tot_size);

#endif
}


bool checkStates(State* s1, State* s2) {
    int ok = 1;
    if(s1->q20bases_ != s2->q20bases_) {
        ok = 0;
        cerr << "GG1" << endl;
    }

    if(s1->q30bases_ != s2->q30bases_) {cerr << "GG2" << endl;ok=0;}
    if(s1->lines_ != s2->lines_) {cerr << "GG3" << endl;ok=0;}

    if(s1->malloc_seq_len_ != s2->malloc_seq_len_) {cerr << "GG4" << endl;ok=0;}

    if(s1->qul_range_ != s2->qul_range_) {cerr << "GG5" << endl;ok=0;}

    if(s1->real_seq_len_ != s2->real_seq_len_) {cerr << "GG6" << endl;ok=0;}

    if(s1->kmer_buf_len_ != s2->kmer_buf_len_) {cerr << "GG7" << endl;ok=0;}

    if(s1->tot_bases_ != s2->tot_bases_) {cerr << "GG8" << endl; ok=0;}

    if(s1->gc_bases_ != s2->gc_bases_) {cerr << "GG9" << endl;ok=0;}

    if(fabs(s1->avg_len - s2->avg_len) / s2->avg_len > 1e-2) {cerr << "GG10" << endl;ok=0;}

    if(s1->pass_reads_ != s2->pass_reads_) {cerr << "GG11" << endl;ok=0;}

    if(s1->fail_short_ != s2->fail_short_) {cerr << "GG12" << endl;ok=0;}

    if(s1->fail_long_ != s2->fail_long_) {cerr << "GG13" << endl;ok=0;}

    if(s1->fail_N_ != s2->fail_N_) {cerr << "GG14" << endl;ok=0;}

    if(s1->fail_lowq_!= s2->fail_lowq_) {cerr << "GG15" << endl;ok=0;}

    if(s1->trim_adapter_ != s2->trim_adapter_) {cerr << "GG16" <<endl;ok=0;}

    if(s1->trim_adapter_bases_ != s2->trim_adapter_bases_) {cerr << "GG17" << endl;ok=0;}

    //pos_qul
    for(int i = 0; i < s1->real_seq_len_; i++) {
        if(s1->pos_qul_[i] != s2->pos_qul_[i]) {
            cerr << "GG18 " << s1->pos_qul_[i] << " " << s2->pos_qul_[i] << endl;
            ok = 0;
        }
    }

    //len_cnt
    for(int i = 0; i < s1->real_seq_len_; i++) {
        if(s1->len_cnt_[i] != s2->len_cnt_[i]) {
            cerr << "GG19 " << s1->len_cnt_[i] << " " << s2->len_cnt_[i] << endl;
            ok = 0;
        }
    }

    //pos_cnt
    for(int i = 0; i < s1->real_seq_len_; i++) { 
        for(int j = 0; j < 4; j++) {
            if(s1->pos_cnt_[i * 4 + j] != s2->pos_cnt_[i * 4 + j]) {
                cerr << "GG20 " << s1->pos_cnt_[i * 4 + j] << " " << s2->pos_cnt_[i * 4 + j] << endl;
                ok = 0;
            }
        }
    }

    //qul_cnt
    for(int i = 0; i < s1->qul_range_; i++) {
        if(s1->qul_cnt_[i] != s2->qul_cnt_[i]) {
            cerr << "GG21" << endl;
            ok = 0;
        }
    }
    //gc_cnt
    for(int i = 0; i < 100; i++) {
        if(s1->gc_cnt_[i] != s2->gc_cnt_[i]) {
            cerr << "GG22" << endl;
            ok = 0;
        }
    }



    if(ok) return 1;
    else return 0;
}


/**
 * @brief do QC for single-end data
 */

void SeQc::ProcessSeFastq() {

    //TODO
    //if(my_rank == 0) block_size = 6 * (1 << 20);
    //else block_size = 4 * (1 << 20);
    auto *fastqPool = new rabbit::fq::FastqDataPool(Q_lim_se * 64, BLOCK_SIZE);
    //rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> queue1(Q_lim_se * 64, 1);
    auto **p_thread_info = new ThreadInfo *[slave_num];
    for (int t = 0; t < slave_num; t++) {
        p_thread_info[t] = new ThreadInfo(cmd_info_, false);
    }
    rabbit::uint32 tmpSize = SWAP1_SIZE;
    if (cmd_info_->seq_len_ <= 200) tmpSize = SWAP2_SIZE;

    if(in_is_zip_) {
        fqFileReader = new rabbit::fq::FastqFileReader(cmd_info_->in_file_name1_, *fastqPool, "", in_is_zip_, tmpSize, start_line_, end_line_, use_in_mem);
    } else {
        fqFileReader = new rabbit::fq::FastqFileReader(cmd_info_->in_file_name1_, *fastqPool, "", in_is_zip_, tmpSize, start_pos_, end_pos_, use_in_mem);
    }
    if(use_in_mem) {
        MPI_Barrier(MPI_COMM_WORLD);
    }

    double ttt = GetTime();

    

    //thread *write_thread;
    //if (cmd_info_->write_data_) {
    //    write_thread = new thread(bind(&SeQc::WriteSeFastqTask, this));
    //}
    //tagg
    //thread producer(bind(&SeQc::ProducerSeFastqTask, this, cmd_info_->in_file_name1_, fastqPool, ref(queue1)));
    thread producer(bind(&SeQc::ProducerSeFastqTask64, this, cmd_info_->in_file_name1_, fastqPool));
    //thread consumer(bind(&SeQc::ConsumerSeFastqTask, this, p_thread_info, fastqPool, ref(queue1)));
    thread consumer(bind(&SeQc::ConsumerSeFastqTask64, this, p_thread_info, fastqPool));
    if (cmd_info_->write_data_) {
        WriteSeFastqTask();
        printf("write%d done\n", my_rank);
    }

    consumer.join();
    printf("consumer done\n");

    producer.join();
    printf("producer done\n");

    //if (cmd_info_->write_data_) {
    //    write_thread->join();
    //    writerDone = 1;
    //}
    printf("all pro write done1\n");

    MPI_Barrier(MPI_COMM_WORLD);
    printf("all pro write done2\n");
    printf("TOT TIME1 %lf\n", GetTime() - ttt);

    double tt00 = GetTime();
#ifdef Verbose
    printf("all thrad done\n");
    printf("now merge thread info\n");
#endif
    vector < State * > pre_vec_state;
    vector < State * > aft_vec_state;

    for (int t = 0; t < slave_num; t++) {
        pre_vec_state.push_back(p_thread_info[t]->pre_state1_);
        aft_vec_state.push_back(p_thread_info[t]->aft_state1_);
    }

    State* pre_state_tmp;
    State* aft_state_tmp;
    if(cmd_info_->do_overrepresentation_) {
    //if(0) {
        pre_state_tmp = State::MergeStatesSlave(pre_vec_state);
        aft_state_tmp = State::MergeStatesSlave(aft_vec_state);
    } else {
         pre_state_tmp = State::MergeStates(pre_vec_state);
         aft_state_tmp = State::MergeStates(aft_vec_state);
    }
    printf("merge1 done\n");
    printf("merge cost %lf\n", GetTime() - tt00);


    tt00 = GetTime();
    vector<State *> pre_state_mpis;
    pre_state_mpis.push_back(pre_state_tmp);
    vector<State *> aft_state_mpis;
    aft_state_mpis.push_back(aft_state_tmp);
    printf("merge2 done\n");

    MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank == 0) {
        for(int i = 1; i < comm_size; i++) {
            int now_size = 0;
            MPI_Recv(&now_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            char* info = new char[now_size];
            MPI_Recv(info, now_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            State* tmp_state = new State(info, now_size, cmd_info_, cmd_info_->seq_len_, cmd_info_->qul_range_, false);
            delete[] info;
            pre_state_mpis.push_back(tmp_state);
        }
    } else {
        string pre_state_is = pre_state_tmp->ParseString();
        State* tmp_state = new State(pre_state_is.c_str(), pre_state_is.length(), cmd_info_, cmd_info_->seq_len_, cmd_info_->qul_range_, false);
        //bool res = checkStates(pre_state_tmp, tmp_state);
        //if(!res) cerr << "GGGGGGGG" << endl;
        int now_size = pre_state_is.length();
        MPI_Send(&now_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(pre_state_is.c_str(), now_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    printf("pre merge done\n");

    if(my_rank == 0) {
        for(int i = 1; i < comm_size; i++) {
            int now_size = 0;
            MPI_Recv(&now_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            char* info = new char[now_size];
            MPI_Recv(info, now_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            State* tmp_state = new State(info, now_size, cmd_info_, cmd_info_->seq_len_, cmd_info_->qul_range_, false);
            delete[] info;
            aft_state_mpis.push_back(tmp_state);
        }
    } else {
        string aft_state_is = aft_state_tmp->ParseString();
        int now_size = aft_state_is.length();
        MPI_Send(&now_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(aft_state_is.c_str(), now_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printf("aft merge done\n");

    printf("TOT TIME3 %lf\n", GetTime() - ttt);

    if(my_rank) {
        delete pre_state_tmp;
        delete aft_state_tmp;
        delete fastqPool;
        for (int t = 0; t < slave_num; t++) {
            delete p_thread_info[t];
        }
        delete[] p_thread_info;
        delete fqFileReader;
        if(use_out_mem) {
            delete[] OutMemData;
        }
        return;
    }


    auto pre_state = State::MergeStates(pre_state_mpis);
    auto aft_state = State::MergeStates(aft_state_mpis);
    printf("merge2 cost %lf\n", GetTime() - tt00);


    tt00 = GetTime();
#ifdef Verbose
    if (cmd_info_->do_overrepresentation_) {
        printf("orp cost %lf\n", pre_state->GetOrpCost() + aft_state->GetOrpCost());
    }
    printf("merge done\n");
#endif
    printf("\nprint read (before filter) info :\n");
    State::PrintStates(pre_state);
    printf("\nprint read (after filter) info :\n");
    State::PrintStates(aft_state);
    printf("\n");
    if (cmd_info_->print_what_trimmed_)
        State::PrintAdapterToFile(aft_state);
    State::PrintFilterResults(aft_state);
    printf("\n");


    if (cmd_info_->do_overrepresentation_ && cmd_info_->print_ORP_seqs_) {
        auto hash_graph1 = pre_state->GetHashGraph();
        int hash_num1 = pre_state->GetHashNum();

        auto hash_graph2 = aft_state->GetHashGraph();
        int hash_num2 = aft_state->GetHashNum();

        int spg = cmd_info_->overrepresentation_sampling_;

        ofstream ofs;
        string srr_name = cmd_info_->in_file_name1_;
        srr_name = PaseFileName(srr_name);

        string out_name = "se_" + srr_name + "_before_ORP_sequences.txt";
        ofs.open(out_name, ifstream::out);
        ofs << "sequence count"
            << "\n";
        int cnt1 = 0;
        for (int i = 0; i < hash_num1; i++) {
            if (!overRepPassed(hash_graph1[i].seq, hash_graph1[i].cnt, spg)) continue;
            ofs << hash_graph1[i].seq << " " << hash_graph1[i].cnt << "\n";
            cnt1++;
        }
        ofs.close();
        printf("in %s (before filter) find %d possible overrepresented sequences (store in %s)\n", srr_name.c_str(),
                cnt1, out_name.c_str());


        out_name = "se_" + srr_name + "_after_ORP_sequences.txt";
        ofs.open(out_name, ifstream::out);
        ofs << "sequence count"
            << "\n";
        int cnt2 = 0;
        for (int i = 0; i < hash_num2; i++) {
            if (!overRepPassed(hash_graph2[i].seq, hash_graph2[i].cnt, spg)) continue;
            ofs << hash_graph2[i].seq << " " << hash_graph2[i].cnt << "\n";
            cnt2++;
        }
        ofs.close();
        printf("in %s (after filter) find %d possible overrepresented sequences (store in %s)\n", srr_name.c_str(),
                cnt2, out_name.c_str());
        printf("\n");
    }


    int *dupHist = NULL;
    double *dupMeanGC = NULL;
    double dupRate = 0.0;
    int histSize = 32;
    if (cmd_info_->state_duplicate_) {
        dupHist = new int[histSize];
        memset(dupHist, 0, sizeof(int) * histSize);
        dupMeanGC = new double[histSize];
        memset(dupMeanGC, 0, sizeof(double) * histSize);
        dupRate = duplicate_->statAll(dupHist, dupMeanGC, histSize);
        printf("Duplication rate (may be overestimated since this is SE data): %.5f %%\n", dupRate * 100.0);
        delete[] dupHist;
        delete[] dupMeanGC;
    }
    printf("print and other cost %lf\n", GetTime() - tt00);

    tt00 = GetTime();
    string srr_name = cmd_info_->in_file_name1_;
    srr_name = PaseFileName(srr_name);
    if(pre_state->GetLines() == 0) {
        printf("read number is 0, don't print html reporter, please check input file\n");
    } else {
        Repoter::ReportHtmlSe(srr_name + "_RabbitQCPlus.html", pre_state, aft_state, cmd_info_->in_file_name1_,
                dupRate * 100.0);
    }

    printf("report cost %lf\n", GetTime() - tt00);
    if(use_in_mem) {
        fprintf(stderr, "TOT TIME %lf\n", GetTime() - ttt);
    }

    tt00 = GetTime();
    delete pre_state_tmp;
    delete aft_state_tmp;

    delete pre_state;
    delete aft_state;

    delete fastqPool;


    for (int t = 0; t < slave_num; t++) {
        delete p_thread_info[t];
    }

    delete[] p_thread_info;
    printf("delete cost %lf\n", GetTime() - tt00);
    delete fqFileReader;
    if(use_out_mem) {
        delete[] OutMemData;
    }
}


/**
 * @brief get fastq data chunk from fastq_data_pool and put it into data queue
 * @param file : fastq file name, which is also the input file
 * @param fastq_data_pool : fastq data pool
 * @param dq : a data queue
 */

void SeQc::ProducerSeFastqTask(string file, rabbit::fq::FastqDataPool *fastq_data_pool,
        rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {

    double t0 = GetTime();
    //cerr << "part_sizes ";
    //for(int i = 0; i < comm_size; i++) cerr << part_sizes[i] << " ";
    //cerr << endl;
    rabbit::fq::FastqFileReader *fqFileReader;
    rabbit::uint32 tmpSize = SWAP1_SIZE;
    if (cmd_info_->seq_len_ <= 200) tmpSize = SWAP2_SIZE;
    fqFileReader = new rabbit::fq::FastqFileReader(file, *fastq_data_pool, "", in_is_zip_, tmpSize);
    int64_t n_chunks = 0;

    if (cmd_info_->use_pugz_) {
        fprintf(stderr, "use pugz TODO...\n");
        exit(0);
    } else {
        double t_sum1 = 0;
        double t_sum2 = 0;
        double t_sum3 = 0;
        bool is_first = 1;
        int64_t tot_size = 0;
        while (true) {
            double tt0 = GetTime();
            rabbit::fq::FastqDataChunk *fqdatachunk;
            int64_t offset;
            if(is_first) {
                offset = 0;
                for(int i = 0; i < my_rank; i++)
                    offset += part_sizes[i];
            } else {
                offset = -1;
            }
            is_first = 0;
            fqdatachunk = fqFileReader->readNextChunk();
            t_sum1 += GetTime() - tt0;
            tt0 = GetTime();
            if (fqdatachunk == NULL) {
                //printf("null\n");
                if(cmd_info_->write_data_) {
                    now_chunks = n_chunks;
                    producerStop = 1;
                    while(consumerCommDone == 0) {
                        usleep(100);
                    }
                    //printf("producer%d get val done %lld %lld\n", my_rank, now_chunks, mx_chunks);
                    int bu_chunks = mx_chunks - now_chunks;
                    if(bu_chunks < 0) {
                        fprintf(stderr, "error bu chunks number\n");
                        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                    }
                    cerr << "bu rank" << my_rank << " : " << bu_chunks << " of (" << now_chunks << ", " << mx_chunks << ")" << endl;
                    if(bu_chunks) {
                        for(int i = 0; i < bu_chunks; i++) {
                            dq.Push(n_chunks, fqdatachunk);
                            n_chunks++;
                        }
                    }
                } 
                break;
            }
            tot_size += fqdatachunk->size;
            if(fqdatachunk->size <= 0 || tot_size <= 0) break;
            //printf("producer%d read %lld\n", my_rank, tot_size);
            dq.Push(n_chunks, fqdatachunk);


            //printf("producer%d %lld done\n", my_rank, n_chunks);
            n_chunks++;
            t_sum2 += GetTime() - tt0;
        }
        printf("producer sum1 cost %lf\n", t_sum1);
        printf("producer sum2 cost %lf\n", t_sum2);
    }

    printf("producer%d cost %lf\n", my_rank, GetTime() - t0);
    dq.SetCompleted();
    delete fqFileReader;
//    producerDone = 1;
}


void SeQc::ConsumerSeFastqTask(ThreadInfo **thread_infos, rabbit::fq::FastqDataPool *fastq_data_pool,
        rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {
    rabbit::int64 id = 0;
    rabbit::fq::FastqDataChunk *fqdatachunk;
    qc_data_tgs para;
    para.cmd_info_ = cmd_info_;
    para.thread_info_ = thread_infos;
    para.bit_len = 0;

    bool proDone = 0;
    if(cmd_info_->state_duplicate_) {
        para.bit_len = duplicate_->key_len_base_;
    }
    athread_init();
    if (cmd_info_->is_TGS_) {
        double t0 = GetTime();
        double tsum1 = 0;
        double tsum2 = 0;
        double tsum3 = 0;
        while (dq.Pop(id, fqdatachunk)) {
            double tt0 = GetTime();
            tsum1 += GetTime() - tt0;
            tt0 = GetTime();
            vector <neoReference> data;
            rabbit::fq::chunkFormat(fqdatachunk, data, true);
            tsum2 += GetTime() - tt0;
            tt0 = GetTime();
            para.data1_ = &data;
            {
                lock_guard<mutex> guard(globalMutex); 
                __real_athread_spawn((void *)slave_tgsfunc, &para, 1);
                athread_join();
            }
            tsum3 += GetTime() - tt0;
            fastq_data_pool->Release(fqdatachunk);
        }
        printf("TGSnew tot cost %lf\n", GetTime() - t0);
        printf("TGSnew producer cost %lf\n", tsum1);
        printf("TGSnew format cost %lf\n", tsum2);
        printf("TGSnew slave cost %lf\n", tsum3);
    } else {
        
    }
    done_thread_number_++;
    athread_halt();
    //printf("consumer done\n");
}


void SeQc::ProcessSeTGS() {

    auto *fastqPool = new rabbit::fq::FastqDataPool(256, BLOCK_SIZE);
    rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> queue1(256, 1);
    auto **p_thread_info = new ThreadInfo *[slave_num];
    for (int t = 0; t < slave_num; t++) {
        p_thread_info[t] = new ThreadInfo(cmd_info_, false);
    }
    thread producer(bind(&SeQc::ProducerSeFastqTask, this, cmd_info_->in_file_name1_, fastqPool, ref(queue1)));
    thread consumer(bind(&SeQc::ConsumerSeFastqTask, this, p_thread_info, fastqPool, ref(queue1)));
    producer.join();
    printf("producer done\n");
    consumer.join();
    printf("consumer done\n");

#ifdef Verbose
    printf("all thrad done\n");
    printf("now merge thread info\n");
#endif
    vector <TGSStats*> vec_state;

    for (int t = 0; t < slave_num; t++) {
        vec_state.push_back(p_thread_info[t]->TGS_state_);
    }
    auto mer_state = TGSStats::merge(vec_state);


#ifdef Verbose
    printf("merge done\n");
#endif
    printf("\nprint TGS state info :\n");

    //report3(mer_state);
    mer_state->CalReadsLens();

    mer_state->print();
    string srr_name = cmd_info_->in_file_name1_;
    srr_name = PaseFileName(srr_name);
    string command = cmd_info_->command_;
    Repoter::ReportHtmlTGS(srr_name + "_RabbitQCPlus.html", command, mer_state, cmd_info_->in_file_name1_);

    delete fastqPool;
    for (int t = 0; t < slave_num; t++) {
        delete p_thread_info[t];
    }
    delete[] p_thread_info;
}

