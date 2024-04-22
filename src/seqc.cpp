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
        ifstream iff_idx;
        iff_idx.open(cmd_info1->in_file_name1_ + ".swidx");
        size_t block_size = 0;
        vector<size_t> block_sizes;
        while (iff_idx >> block_size) {
            block_sizes.push_back(block_size);
        }
        iff_idx.close();
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
                string idx_name = cmd_info1->out_file_name1_ + ".swidx";
                //fprintf(stderr, "idx_name %s\n", idx_name.c_str());
                off_idx.open(idx_name);

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
        if (fqdatachunk == NULL) {
            double tt1 = GetTime();
            if(tmp_chunks.size()) {
                //fprintf(stderr, "put last %d\n", tmp_chunks.size());
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
            fprintf(stderr, "producer rank%d stop, tot chunk size %d\n", my_rank, n_chunks);
            if(cmd_info_->write_data_) {
                tmp_chunks.clear();
                while (p_queueNumNow >= p_queueSizeLim) {
                    usleep(100);
                }
                p_out_queue_[p_queueP2++] = tmp_chunks;
                p_queueNumNow++;
            } 
            break;
        }
        t_sum2 += GetTime() - tt0;

        tt0 = GetTime();
        tot_size += fqdatachunk->size;
        if(fqdatachunk->size <= 0 || tot_size <= 0) {
            if(tmp_chunks.size()) {
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
            break;
        }
        tmp_chunks.push_back(fqdatachunk); 
        if(tmp_chunks.size() == 64) {
            //fprintf(stderr, "put 64\n");
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

size_t se_rank_size = 0;

//void SeQc::ProcessNgsData(bool &proDone, vector <neoReference> &data, rabbit::fq::FastqDataChunk *fqdatachunk, qc_data *para, rabbit::fq::FastqDataPool *fastq_data_pool) {
void SeQc::ProcessNgsData(bool &allIsNull, vector <neoReference> data[64], vector <rabbit::fq::FastqDataChunk *> fqdatachunks, qc_data *para, rabbit::fq::FastqDataPool *fastq_data_pool) {

    //fprintf(stderr, "rank%d data size %d\n", my_rank, data.size());
    double tt0 = GetTime();
    vector <neoReference> pass_data[64];
    int isNull = 1;
    for(int i = 0; i < 64; i++) {
        se_rank_size += data[i].size();
        if(data[i].size()) isNull = 0;
        if(cmd_info_->write_data_) {
            pass_data[i].resize(data[i].size());
        }
        para->data1_[i] = &data[i];
        para->pass_data1_[i] = &pass_data[i];
    }
    fprintf(stderr, "consumer%d isNull %d\n", my_rank, isNull);
    tsum1 += GetTime() - tt0;

    tt0 = GetTime();
    {
        lock_guard<mutex> guard(globalMutex);
        __real_athread_spawn((void *)slave_ngsfunc, para, 1);
        athread_join();
    }
    tsum2 += GetTime() - tt0;


    tt0 = GetTime();
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
    int all_stop = 1;
    for(int ii = 0; ii < comm_size; ii++) {
        if(isNulls[ii] == 0) all_stop = 0;
    }
    tsum3 += GetTime() - tt0;

    if(all_stop) {
        allIsNull = 1;
        return;
    }


    tt0 = GetTime();
    if (cmd_info_->write_data_) {
        double tt1 = GetTime();
        double tt2 = GetTime();
        int out_lens[64] = {0};
        for(int i = 0; i < 64; i++) {
            for(auto item : pass_data[i]){
                if(item.lname == 0) continue;
                out_lens[i] += item.lname + item.lseq + item.lstrand + item.lqual + 4;
            }
        }

        tsum4_1_1 += GetTime() - tt2;


        tt2 = GetTime();
        char *out_data[64];
        for(int i = 0; i < 64; i++) {
            if(out_lens[i]) out_data[i] = new char[out_lens[i]];
            else out_data[i] = (char *)(NULL);
//            if(out_lens[i]) memset(out_data[i], 0, sizeof(char) * out_lens[i]);
        }
        tsum4_1_2 += GetTime() - tt2;
        
        tt2 = GetTime();
        //OPT: reduce memcpy times, merge contiguous datas.

        for(int ii = 0; ii < 64; ii++) {
            int tmp_len = 0;
            int pos = 0;
            char* start_pos = NULL;
            for(int i = 0; i < pass_data[ii].size(); i++) {
                if(pass_data[ii][i].lname == 0 || pass_data[ii][i].pname + pass_data[ii][i].lname + pass_data[ii][i].lseq + pass_data[ii][i].lstrand + 3 != pass_data[ii][i].pqual) {
                    //this record has been trimmed
                    if(tmp_len) memcpy(out_data[ii] + pos, start_pos, tmp_len);
                    pos += tmp_len;
                    if(pass_data[ii][i].lname) Read2Chars(pass_data[ii][i], out_data[ii], pos);
                    tmp_len = 0;
                    if(i < pass_data[ii].size() - 1) {
                        start_pos = (char*)pass_data[ii][i + 1].base + pass_data[ii][i + 1].pname;
                    }
                    continue;
                }
                if((char*)pass_data[ii][i].base + pass_data[ii][i].pname != start_pos + tmp_len) {
                    //record is complete, but can extend to last record
                    if(tmp_len) memcpy(out_data[ii] + pos, start_pos, tmp_len);
                    pos += tmp_len;
                    tmp_len = 0;
                    start_pos = (char*)pass_data[ii][i].base + pass_data[ii][i].pname;
                }
                tmp_len += pass_data[ii][i].lname + pass_data[ii][i].lseq + pass_data[ii][i].lstrand + pass_data[ii][i].lqual + 4;
            }
            if(tmp_len) {
                memcpy(out_data[ii] + pos, start_pos, tmp_len);
                pos += tmp_len;
                tmp_len == 0;
            }
            ASSERT(pos == out_lens[ii]);
        }



//        for(int i = 0; i < 64; i++) {
//            int pos = 0;
//            for (auto item: pass_data[i]) {
//                if(item.lname == 0) continue;
//                Read2Chars(item, out_data[i], pos);
//            }
//            ASSERT(pos == out_lens[i]);
//        }


        tsum4_1_3 += GetTime() - tt2;

        tsum4_1 += GetTime() - tt1;

        long long now_pos_base = 0;
        bool use_consumer_libdeflate = 0;

#ifdef CONSUMER_USE_LIBDEFLATE
        use_consumer_libdeflate = 1;
#endif
        if(!out_is_zip_ || use_consumer_libdeflate == 0) {
            tt1 = GetTime();

            tt2 = GetTime();
            int out_len = 0;
            for(int i = 0; i < 64; i++) {
                out_len += out_lens[i];
            }
            int now_size = out_len;
            int now_sizes[comm_size];
            now_sizes[0] = now_size;
            tsum4_2_1 += GetTime() - tt2;

            tt2 = GetTime();
            MPI_Barrier(MPI_COMM_WORLD);
            tsum4_2_2 += GetTime() - tt2;

            tt2 = GetTime();
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
            tsum4_2_3 += GetTime() - tt2;
            tsum4_2 += GetTime() - tt1;
        }



        tt1 = GetTime();



        tsum4_3 += GetTime() - tt1;


        tt1 = GetTime();
        for(int i = 0; i < 64; i++) {
            writeInfos.push_back({out_data[i], out_lens[i], now_pos_base});
        }
//
//        //mylock.lock();
//        while (queueNumNow >= queueSizeLim) {
//#ifdef Verbose
//            //printf("waiting to push a chunk to out queue %d\n",out_len);
//#endif
//            usleep(100);
//        }
//        out_queue_[queueP2++] = make_pair(out_data, make_pair(out_len, now_pos_base));
//        queueNumNow++;
//        //mylock.unlock();
        tsum4_4 += GetTime() - tt1;

    }
    tsum4 += GetTime() - tt0;

    tt0 = GetTime();
    if(isNull) {
        fprintf(stderr, "consumer%d push null\n", my_rank);
        vector<rabbit::fq::FastqDataChunk *> tmp_chunks;
        tmp_chunks.clear();
        while (p_queueNumNow >= p_queueSizeLim) {
            usleep(100);
        }
        p_out_queue_[p_queueP2++] = tmp_chunks;
        p_queueNumNow++;
    }
    tsum5 += GetTime() - tt0;

    tt0 = GetTime();
    for(int i = 0; i < 64; i++) {
        if(fqdatachunks[i] != NULL) fastq_data_pool->Release(fqdatachunks[i]);
    }
    tsum6 += GetTime() - tt0;

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
    double t_wait_producer = 0;
    double t_format = 0;
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
    vector<rabbit::fq::FastqDataChunk *> fqdatachunks;
    while(true) {
        double tt0 = GetTime();
        while (p_queueNumNow == 0) {
            usleep(1000);
        }
        t_wait_producer += GetTime() - tt0;

        tt0 = GetTime();
        fqdatachunks = p_out_queue_[p_queueP1++];
        p_queueNumNow--;
        vector<neoReference> data[64];
        format_data para2[64];
        uint64_t counts[64] = {0};
        for(int i = fqdatachunks.size(); i < 64; i++) {
            fqdatachunks.push_back(NULL);
        }
        int fq_size = fqdatachunks.size();
        for(int i = 0; i < 64; i++) {
            para2[i].fqSize = fq_size;
            para2[i].my_rank = my_rank;
        }
//        fprintf(stderr, "consumer rank%d fq_size %d\n", my_rank, fq_size);
        for(int i = 0; i < fq_size; i++) {
            //if(fqdatachunks[i] != NULL) fprintf(stderr, "rank%d addr %p, size %d\n", my_rank, fqdatachunks[i], fqdatachunks[i]->size);
            if(fqdatachunks[i] != NULL) data[i].resize(fqdatachunks[i]->size / (cmd_info_->seq_len_ * 2));
            para2[i].data = &(data[i]);
            para2[i].fqdatachunk = fqdatachunks[i];
            para2[i].res = counts;
            para2[i].fqSize = fq_size;
        }
        {
            lock_guard<mutex> guard(globalMutex);
            __real_athread_spawn((void *)slave_formatfunc, para2, 1);
            athread_join();
        }
        t_format += GetTime() - tt0;

        tt0 = GetTime();
        writeInfos.clear();
        for(int i = 0; i < fq_size; i++) {
            //fprintf(stderr, "rank%d count %d\n", my_rank, counts[i]);
            data[i].resize(counts[i]);
//            ProcessNgsData(allProDone, data[i], fqdatachunks[i], &para, fastq_data_pool);
        }
        t_resize += GetTime() - tt0;

        tt0 = GetTime();
        ProcessNgsData(allIsNull, data, fqdatachunks, &para, fastq_data_pool);
        fqdatachunks.clear();
        t_ngsfunc += GetTime() - tt0;

        if(allIsNull) break;

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
            //TODO
//            for(int i = 0; i < 64; i++) {
//                off_idx << out_size[i] << endl;
//            }
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
            //QChunkItem Qitem;
            //for(int i = 0; i < writeInfos.size(); i++) {
            //    Qitem.buffer[i] = writeInfos[i].buffer;
            //    Qitem.buffer_len[i] = writeInfos[i].buffer_len;
            //    Qitem.file_offset[i] = writeInfos[i].file_offset;
            //}

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
    }


    fprintf(stderr, "consumer%d NGSnew tot cost %lf\n", my_rank, GetTime() - t0);
    fprintf(stderr, "consumer%d NGSnew wait producer cost %lf\n", my_rank, t_wait_producer);
    fprintf(stderr, "consumer%d NGSnew format cost %lf\n", my_rank, t_format);
    fprintf(stderr, "consumer%d NGSnew resize cost %lf\n", my_rank, t_resize);
    fprintf(stderr, "consumer%d NGSnew ngsfunc cost %lf\n", my_rank, t_ngsfunc);
    fprintf(stderr, "consumer%d NGSnew resize cost %lf\n", my_rank, tsum1);
    fprintf(stderr, "consumer%d NGSnew slave cost %lf\n", my_rank, tsum2);
    fprintf(stderr, "consumer%d NGSnew isnull comm cost %lf\n", my_rank, tsum3);
    fprintf(stderr, "consumer%d NGSnew write cost %lf (%lf [%lf %lf %lf], %lf [%lf %lf %lf], %lf [%lf %lf %lf %lf %lf], %lf)\n", my_rank, tsum4, tsum4_1, tsum4_1_1, tsum4_1_2, tsum4_1_3, tsum4_2, tsum4_2_1, tsum4_2_2, tsum4_2_3, tsum4_3, tsum4_3_1, tsum4_3_2, tsum4_3_3, tsum4_3_4, tsum4_3_5, tsum4_4);
    fprintf(stderr, "consumer%d NGSnew push null cost %lf\n", my_rank, tsum5);
    fprintf(stderr, "consumer%d NGSnew release cost %lf\n", my_rank, tsum6);
    fprintf(stderr, "consumer%d NGSnew gz slave cost %lf\n", my_rank, t_slave_gz);
    fprintf(stderr, "consumer%d NGSnew gz slave sub cost [%lf %lf %lf %lf %lf %lf %lf]\n", my_rank, t_slave_gz1, t_slave_gz2, t_slave_gz3_1, t_slave_gz3_2, t_slave_gz3_3, t_slave_gz4, t_slave_gz5);
    fprintf(stderr, "consumer%d NGSnew push to queue cost %lf\n", my_rank, t_push_q);
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

//#ifdef WRITER_USE_LIBDEFLATE
//
//        tt0 = GetTime();
//        if (out_is_zip_) {
//            Para paras[64];
//            size_t out_size[64] = {0};
//            for(int i = 0; i < 64; i++) {
//                paras[i].in_buffer = Qitem.buffer[i];
//                paras[i].out_buffer = new char[BLOCK_SIZE];
//                paras[i].in_size = Qitem.buffer_len[i];
//                paras[i].out_size = &(out_size[i]);
//                paras[i].level = 1;
//            }
//            {
//                lock_guard<mutex> guard(globalMutex);
//                __real_athread_spawn((void *)slave_compressfunc, paras, 1);
//                athread_join();
//            }
//            for(int i = 0; i < 64; i++) {
//                off_idx << out_size[i] << endl;
//            }
//            for(int i = 0; i < 64; i++) {
//                delete[] Qitem.buffer[i];
//                Qitem.buffer[i] = paras[i].out_buffer;
//                Qitem.buffer_len[i] = out_size[i];
//
//                int now_size = Qitem.buffer_len[i];
//                int now_sizes[comm_size];
//                now_sizes[0] = now_size;
//                MPI_Barrier(MPI_COMM_WORLD);
//                if(my_rank) {
//                    MPI_Send(&now_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
//                } else {
//                    for(int ii = 1; ii < comm_size; ii++) {
//                        MPI_Recv(&(now_sizes[ii]), 1, MPI_INT, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//                    }
//                }
//                MPI_Barrier(MPI_COMM_WORLD);
//                if(my_rank == 0) {
//                    for(int ii = 1; ii < comm_size; ii++) {
//                        MPI_Send(now_sizes, comm_size, MPI_INT, ii, 0, MPI_COMM_WORLD);
//                    }
//                } else {
//                    MPI_Recv(now_sizes, comm_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//                }
//                MPI_Barrier(MPI_COMM_WORLD);
//                int pre_sizes[comm_size];
//                pre_sizes[0] = 0;
//                for(int ii = 1; ii < comm_size; ii++) {
//                    pre_sizes[ii] = pre_sizes[ii - 1] + now_sizes[ii - 1];
//                }
//                long long now_pos_base = zip_now_pos_ + pre_sizes[my_rank];
//                MPI_Barrier(MPI_COMM_WORLD);
//                for(int ii = 0; ii < comm_size; ii++) {
//                    zip_now_pos_ += now_sizes[ii];
//                }
//                Qitem.file_offset[i] = now_pos_base;
//
//            }
//        }
//        t_gz_slave += GetTime() - tt0;
//
//#endif

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
            off_idx.close();

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
    auto pre_state_tmp = State::MergeStates(pre_vec_state);
    auto aft_state_tmp = State::MergeStates(aft_vec_state);
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
    if(use_in_mem) {
        fprintf(stderr, "TOT TIME %lf\n", GetTime() - ttt);
    }
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
        athread_init();
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

