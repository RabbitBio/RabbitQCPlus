// 
// Created by ylf9811 on 2021/7/6.
//

#include "peqc.h"

using namespace std;


extern "C" {
#include <athread.h>
#include <pthread.h>
    void slave_ngspefunc();
    void slave_formatpefunc();
    void slave_peallfunc();
    void slave_decompressfunc();
    void slave_compressfunc();
}




void PrintRef(neoReference &ref) {
    printf("%.*s\n", ref.lname, (char *) ref.base + ref.pname);
    printf("%.*s\n", ref.lseq, (char *) ref.base + ref.pseq);
    printf("%.*s\n", ref.lstrand, (char *) ref.base + ref.pstrand);
    printf("%.*s\n", ref.lqual, (char *) ref.base + ref.pqual);
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

static bool use_in_mem = 0;
static bool use_out_mem = 0;

static bool use_swidx_file = 0;

#define use_mpi_file

char* OutMemData1;
char* OutMemData2;



/**
 * @brief Construct function
 * @param cmd_info1 : cmd information
 */
PeQc::PeQc(CmdInfo *cmd_info1, int my_rank_, int comm_size_) {
    //printf(" %d /// %d\n", my_rank_, comm_size_);
    my_rank = my_rank_;
    comm_size = comm_size_;
    part_sizes = new int64_t[comm_size];
    now_pos1_ = 0;
    now_pos2_ = 0;
    zip_now_pos1_ = 0;
    zip_now_pos2_ = 0;
    cmd_info_ = cmd_info1;
    filter_ = new Filter(cmd_info1);
    done_thread_number_ = 0;
    int out_block_nums = int(1.0 * cmd_info1->in_file_size1_ / cmd_info1->out_block_size_);
    //printf("out_block_nums %d\n", out_block_nums);
    out_queue_ = NULL;
    p_out_queue_ = NULL;
    //out_queue1_ = NULL;
    //out_queue2_ = NULL;
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

    gFile.open(cmd_info1->in_file_name2_.c_str());
    gFile.seekg(0, ios_base::end);
    long long real_file_size2 = gFile.tellg();
    gFile.close();

    if(!in_is_zip_ && real_file_size != real_file_size2 && comm_size > 1) {
        fprintf(stderr, "PE in files size error %lld %lld\n", real_file_size, real_file_size2);
        exit(0);
    }


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
        int res_size = fread(tmp_chunk, sizeof(char), 1 << 20, pre_fp);
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
    //MPI_Barrier(MPI_COMM_WORLD);
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

    p_out_queue_ = new vector<rabbit::fq::FastqDataPairChunk *>[1 << 20];
    p_queueP1 = 0;
    p_queueP2 = 0;
    p_queueNumNow = 0;
    p_queueSizeLim = Q_lim_pe;

 
    if (cmd_info1->write_data_) {

        out_queue_ = new pair<CIPair, CIPair>[1 << 20];
//        out_queue_ = new pair<QChunkItem, QChunkItem>[1 << 20];
        queueP1 = 0;
        queueP2 = 0;
        queueNumNow = 0;
        queueSizeLim = Q_lim_pe * 64;
        //queueSizeLim = Q_lim_pe;

        //out_queue1_ = new CIPair[1 << 20];
        //queue1P1 = 0;
        //queue1P2 = 0;
        //queueNumNow1 = 0;
        //queueSizeLim1 = Q_lim_pe;
        if (cmd_info_->interleaved_out_ == 0) {
            //out_queue2_ = new CIPair[1 << 20];
            //queue2P1 = 0;
            //queue2P2 = 0;
            //queueNumNow2 = 0;
            //queueSizeLim2 = Q_lim_pe;
        }
        if(use_out_mem) {
            OutMemData1 = new char[128 << 20];
            OutMemData2 = new char[128 << 20];
        }

        if (out_is_zip_) {
            if (cmd_info1->use_pigz_) {
                fprintf(stderr, "use pigz TODO...\n");
                exit(0);
            } else {
#ifdef Verbose
                printf("open gzip stream1 %s\n", cmd_info1->out_file_name1_.c_str());
                printf("open gzip stream2 %s\n", cmd_info1->out_file_name2_.c_str());
#endif

#ifdef USE_LIBDEFLATE
                //string idx_name = cmd_info1->out_file_name1_ + ".swidx";
                //off_idx1.open(idx_name);

                ofstream file1(cmd_info1->out_file_name1_.c_str());

                if (file1.is_open()) {
                    file1 << "Hello, file1!" << endl;
                    file1.close();
                } else {
                    cerr << "Failed to open file1: " << cmd_info1->out_file_name1_ << endl;
                }

                int err1 = MPI_File_open(MPI_COMM_WORLD, cmd_info1->out_file_name1_.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh1);

                if (err1 != MPI_SUCCESS) {
                    char error_string[BUFSIZ];
                    int length_of_error_string, error_class;

                    MPI_Error_class(err1, &error_class);
                    MPI_Error_string(error_class, error_string, &length_of_error_string);

                    fprintf(stderr, "%3d: Couldn't open file1, MPI Error: %s\n", my_rank, error_string);

                    MPI_Abort(MPI_COMM_WORLD, err1);
                    exit(0);
                }

                //string idx_name2 = cmd_info1->out_file_name2_ + ".swidx";
                //off_idx2.open(idx_name2);

                ofstream file2(cmd_info1->out_file_name2_.c_str());

                if (file2.is_open()) {
                    file2 << "Hello, file2!" << endl;
                    file2.close();
                } else {
                    cerr << "Failed to open file2: " << cmd_info1->out_file_name2_ << endl;
                }

                int err2 = MPI_File_open(MPI_COMM_WORLD, cmd_info1->out_file_name2_.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh2);

                if (err2 != MPI_SUCCESS) {
                    char error_string[BUFSIZ];
                    int length_of_error_string, error_class;

                    MPI_Error_class(err2, &error_class);
                    MPI_Error_string(error_class, error_string, &length_of_error_string);

                    fprintf(stderr, "%3d: Couldn't open file, MPI Error: %s\n", my_rank, error_string);

                    MPI_Abort(MPI_COMM_WORLD, err2);
                    exit(0);
                }
#else
                zip_out_stream1 = gzopen(cmd_info1->out_file_name1_.c_str(), "w");
                gzsetparams(zip_out_stream1, cmd_info1->compression_level_, Z_DEFAULT_STRATEGY);
                gzbuffer(zip_out_stream1, 1024 * 1024);
                zip_out_stream2 = gzopen(cmd_info1->out_file_name2_.c_str(), "w");
                gzsetparams(zip_out_stream2, cmd_info1->compression_level_, Z_DEFAULT_STRATEGY);
                gzbuffer(zip_out_stream2, 1024 * 1024);
#endif
            }
        } else {
#ifdef Verbose
            printf("open stream1 %s\n", cmd_info1->out_file_name1_.c_str());
            if (cmd_info_->interleaved_out_ == 0)
                printf("open stream2 %s\n", cmd_info1->out_file_name2_.c_str());
#endif



#ifdef use_mpi_file

            ofstream file1(cmd_info1->out_file_name1_.c_str());

            if (file1.is_open()) {
                file1 << "Hello, file!" << endl;
                file1.close();
            } else {
                cerr << "Failed to open file1: " << cmd_info1->out_file_name1_ << endl;
            }

            int err1 = MPI_File_open(MPI_COMM_WORLD, cmd_info1->out_file_name1_.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh1);

            if (err1 != MPI_SUCCESS) {
                char error_string[BUFSIZ];
                int length_of_error_string, error_class;

                MPI_Error_class(err1, &error_class);
                MPI_Error_string(error_class, error_string, &length_of_error_string);

                fprintf(stderr, "%3d: Couldn't open file, MPI Error: %s\n", my_rank, error_string);

                MPI_Abort(MPI_COMM_WORLD, err1);
                exit(0);
            }

            if (cmd_info_->interleaved_out_ == 0) {

                ofstream file2(cmd_info1->out_file_name2_.c_str());

                if (file2.is_open()) {
                    file2 << "Hello, file2!" << endl;
                    file2.close();
                } else {
                    cerr << "Failed to open file2: " << cmd_info1->out_file_name2_ << endl;
                }

                int err2 = MPI_File_open(MPI_COMM_WORLD, cmd_info1->out_file_name2_.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh2);

                if (err2 != MPI_SUCCESS) {
                    char error_string[BUFSIZ];
                    int length_of_error_string, error_class;

                    MPI_Error_class(err2, &error_class);
                    MPI_Error_string(error_class, error_string, &length_of_error_string);

                    fprintf(stderr, "%3d: Couldn't open file, MPI Error: %s\n", my_rank, error_string);

                    MPI_Abort(MPI_COMM_WORLD, err2);
                    exit(0);
                }               
            }


#else
            if(my_rank == 0) {
                printf("pre real file1 size %lld\n", real_file_size);
                int fd = open(cmd_info1->out_file_name1_.c_str(), O_CREAT | O_TRUNC | O_RDWR | O_EXCL, 0644);
                ftruncate(fd, sizeof(char) * real_file_size);
                out_stream1_ = fdopen(fd, "w");
            } else {
                int found = 0;
                do {
                    if(-1 == access(cmd_info1->out_file_name1_.c_str(), F_OK)) {
                        if(ENOENT == errno) {
                            //cerr << "waiting file be created..." << endl;
                            usleep(10000);
                        } else {
                            //cerr << "file open GG" << endl;
                            usleep(10000);
                            //exit(0);
                        }
                    } else {
                        out_stream1_ = fopen(cmd_info1->out_file_name1_.c_str(), "r+b");   
                        found = 1;
                    }
                } while(found == 0);
            }

            if (cmd_info_->interleaved_out_ == 0){

                if(my_rank == 0) {
                    printf("pre real file2 size %lld\n", real_file_size2);
                    int fd = open(cmd_info1->out_file_name2_.c_str(), O_CREAT | O_TRUNC | O_RDWR | O_EXCL, 0644);
                    ftruncate(fd, sizeof(char) * real_file_size2);
                    out_stream2_ = fdopen(fd, "w");
                } else {
                    int found = 0;
                    do {

                        if(-1 == access(cmd_info1->out_file_name2_.c_str(), F_OK)) {
                            if(ENOENT == errno) {
                                //cerr << "waiting file be created..." << endl;
                                usleep(10000);
                            } else {
                                //cerr << "file open GG" << endl;
                                usleep(10000);
                                //exit(0);
                            }
                        } else {
                            out_stream2_ = fopen(cmd_info1->out_file_name2_.c_str(), "r+b");   
                            found = 1;
                        }
                    } while(found == 0);
                }
                //MPI_Barrier(MPI_COMM_WORLD);
            }
#endif
        }
    }
    duplicate_ = NULL;
    if (cmd_info1->state_duplicate_) {
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

    pugzDone1 = 0;
    pugzDone2 = 0;
//    producerDone = 0;
    writerDone1 = 0;
    writerDone2 = 0;
//    consumerCommDone = 0;
//    producerStop = 0;
//    now_chunks = 0;
//    mx_chunks = 0;
}


PeQc::~PeQc() {
    delete[] part_sizes;
    delete filter_;
    delete[] p_out_queue_;
    if (cmd_info_->write_data_) {
        delete []out_queue_;
        //delete out_queue1_;
        //if (cmd_info_->interleaved_out_ == 0)
        //    delete out_queue2_;
    }
    if(in_is_zip_) {
#ifdef USE_CC_GZ
        for(int i = 0; i < 64; i++) {
            delete []cc_gz_in_buffer[i];
        }
#endif
    }
    if (cmd_info_->state_duplicate_) {
        delete duplicate_;
    }
    if (cmd_info_->add_umi_) {
        delete umier_;
    }
}


string PeQc::Read2String(neoReference &ref) {
    return string((char *) ref.base + ref.pname, ref.lname) + "\n" +
        string((char *) ref.base + ref.pseq, ref.lseq) + "\n" +
        string((char *) ref.base + ref.pstrand, ref.lstrand) + "\n" +
        string((char *) ref.base + ref.pqual, ref.lqual) + "\n";
}

void PeQc::Read2Chars(neoReference &ref, char *out_data, int &pos) {
    //PrintRef(ref);
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


static rabbit::fq::FastqFileReader *fqFileReader;

void PeQc::ProducerPeFastqTask64(string file, string file2, rabbit::fq::FastqDataPool *fastqPool) {
    double t0 = GetTime();
    rabbit::uint32 tmpSize = SWAP1_SIZE;
    if (cmd_info_->seq_len_ <= 200) tmpSize = SWAP2_SIZE;
    int64_t n_chunks = 0;

    double t_sum1 = 0;
    double t_sum2 = 0;
    bool is_first = 1;
    int64_t tot_size = 0;
    vector<rabbit::fq::FastqDataPairChunk *> tmp_chunks;
    while (true) {
        double tt0 = GetTime();
        rabbit::fq::FastqDataPairChunk *fqdatachunk;
        int64_t offset;
        if(is_first) {
            offset = 0;
            for(int i = 0; i < my_rank; i++)
                offset += part_sizes[i];
        } else {
            offset = -1;
        }
        is_first = 0;
        fqdatachunk = fqFileReader->readNextPairChunk();
        t_sum1 += GetTime() - tt0;

        tt0 = GetTime();
        if ((fqdatachunk == NULL) || (fqdatachunk != NULL && fqdatachunk->left_part->size == 1ll << 32)) {
            if(tmp_chunks.size()) {
//                fprintf(stderr, "producer %d put last %d\n", my_rank, tmp_chunks.size());
                //p_mylock.lock();
                while (p_queueNumNow >= p_queueSizeLim) {
                    usleep(1000);
                    //fprintf(stderr, "producer %d waiting2\n", my_rank);
                }
                p_out_queue_[p_queueP2++] = tmp_chunks;
                p_queueNumNow++;
                //p_mylock.unlock();
                tmp_chunks.clear();
                n_chunks++;
            }
//            producerStop = 1;
//            fprintf(stderr, "producer %d stop, tot chunk size %d\n", my_rank, n_chunks);
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
        tot_size += fqdatachunk->left_part->size;
        tmp_chunks.push_back(fqdatachunk); 
//        fprintf(stderr, "producer %d read a chunk\n", my_rank);
        if(tmp_chunks.size() == 64) {
//            fprintf(stderr, "producer %d put 64\n", my_rank);
            //p_mylock.lock();
            while (p_queueNumNow >= p_queueSizeLim) {
                usleep(1000);
                //fprintf(stderr, "producer %d waiting1\n", my_rank);
            }
            p_out_queue_[p_queueP2++] = tmp_chunks;
            p_queueNumNow++;
            //p_mylock.unlock();
            tmp_chunks.clear();
            n_chunks++;
        }
        t_sum2 += GetTime() - tt0;
    }
    //printf("totsize %lld\n", tot_size);

    fprintf(stderr, "producer%d sum1 cost %lf\n", my_rank, t_sum1);
    fprintf(stderr, "producer%d sum2 cost %lf\n", my_rank, t_sum2);
    fprintf(stderr, "producer%d cost %lf\n", my_rank, GetTime() - t0);
    fprintf(stderr, "producer %d in func done\n", my_rank);
    //dq.SetCompleted();
//    producerDone = 1;
}



struct formatpe_data {
    vector <neoReference> *data1;
    vector <neoReference> *data2;
    rabbit::fq::FastqDataPairChunk *fqdatachunk;
    uint64_t *res;
    int fqSize;
    int my_rank;
};

struct PeMerge_data {
    char* out_data1[64];
    char* out_data2[64];
    int out_lens1[64] = {0};
    int out_lens2[64] = {0};
    std::vector <neoReference> *pass_data1[64];
    std::vector <neoReference> *pass_data2[64];
};

struct PeAll_data{
    formatpe_data para1[64];
    qc_data *para2;
    PeMerge_data *para3;
    int out_len_slave1[64];
    int out_len_slave2[64];
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
double tsum4_3 = 0;
double tsum4_4 = 0;
double tsum5 = 0;
double tsum6 = 0;
double tsum7 = 0;
double t_wait_producer = 0;
double t_format = 0;
double t_resize = 0;
double t_decom = 0;
double t_ngsfunc = 0;
double t_slave_gz = 0;
double t_slave_gz2 = 0;
double t_push_q = 0;


struct ConsumerWriteInfoPE{
    char* buffer;
    char* buffer2;
    int buffer_len;
    int buffer_len2;
    long long file_offset;
    long long file_offset2;
};

vector<ConsumerWriteInfoPE> writeInfosPE;

static int pe_pre_out_len_slave1[64];
static int pe_pre_out_len_slave2[64];

void PeQc::ProcessFormatQCWrite(bool &allIsNull, vector <neoReference> *data1, vector <neoReference> *data2,
                          vector <neoReference> *pass_data1, vector <neoReference> *pass_data2,
                          vector <neoReference> *pre_pass_data1, vector <neoReference> *pre_pass_data2,
                          vector <rabbit::fq::FastqDataPairChunk *> fqdatachunks, vector <rabbit::fq::FastqDataPairChunk *> pre_fqdatachunks,
                          qc_data *para, rabbit::fq::FastqDataPool *fastq_data_pool) {
    double tt0 = GetTime();
    formatpe_data para2[64];
    uint64_t counts[64] = {0};
    int fq_size = fqdatachunks.size();
    for(int i = 0; i < 64; i++) {
        data1[i].clear();
        data2[i].clear();
        para2[i].fqSize = fq_size;
        para2[i].my_rank = my_rank;
    }
    for(int i = 0; i < fq_size; i++) {
        if(fqdatachunks[i] != NULL) {
            data1[i].resize(fqdatachunks[i]->left_part->size / (cmd_info_->seq_len_ * 2));
            data2[i].resize(fqdatachunks[i]->right_part->size / (cmd_info_->seq_len_ * 2));
            if(cmd_info_->write_data_) {
                pass_data1[i].resize(data1[i].size());
                pass_data2[i].resize(data2[i].size());
            }
        } else {
            data1[i].resize(0);
            data2[i].resize(0);
            if(cmd_info_->write_data_) {
                pass_data1[i].resize(0);
                pass_data2[i].resize(0);
            }
        }
        para2[i].data1 = &(data1[i]);
        para2[i].data2 = &(data2[i]);
        para2[i].fqdatachunk = fqdatachunks[i];
        para2[i].res = counts;
        para2[i].fqSize = fq_size;
    }
    for(int i = 0; i < 64; i++) {
        para->data1_[i] = &data1[i];
        para->data2_[i] = &data2[i];
        para->pass_data1_[i] = &pass_data1[i];
        para->pass_data2_[i] = &pass_data2[i];
    }
    tsum1 += GetTime() - tt0;

    tt0 = GetTime();
    PeMerge_data pe_merge_data;
    PeAll_data pe_all_data;
    pe_all_data.write_data = cmd_info_->write_data_;
    int out_lens1[64] = {0};
    int out_lens2[64] = {0};
    char *out_data1[64];
    char *out_data2[64];
    for(int i = 0; i < 64; i++) {
        pe_all_data.para1[i] = para2[i];
        pe_all_data.out_len_slave1[i] = 0;
        pe_all_data.out_len_slave2[i] = 0;
        if(cmd_info_->write_data_) {
            out_lens1[i] = pe_pre_out_len_slave1[i];
            if(out_lens1[i]) out_data1[i] = new char[out_lens1[i]];
            else out_data1[i] = (char *)(NULL);
            out_lens2[i] = pe_pre_out_len_slave2[i];
            if(out_lens2[i]) out_data2[i] = new char[out_lens2[i]];
            else out_data2[i] = (char *)(NULL);
        }
    }
    for(int i = 0; i < 64; i++) {
        pe_merge_data.out_data1[i] = out_data1[i];
        pe_merge_data.out_lens1[i] = out_lens1[i];
        pe_merge_data.pass_data1[i] = &(pre_pass_data1[i]);
        pe_merge_data.out_data2[i] = out_data2[i];
        pe_merge_data.out_lens2[i] = out_lens2[i];
        pe_merge_data.pass_data2[i] = &(pre_pass_data2[i]);
    }
    pe_all_data.para2 = para;
    pe_all_data.para3 = &pe_merge_data;
    tsum2 += GetTime() - tt0;

    tt0 = GetTime();
    __real_athread_spawn((void *)slave_peallfunc, &pe_all_data, 1);
    athread_join();
    tsum3 += GetTime() - tt0;


    tt0 = GetTime();
    if (cmd_info_->write_data_) {
        long long now_pos_base1 = 0;
        long long now_pos_base2 = 0;
        bool use_consumer_libdeflate = 0;

#ifdef CONSUMER_USE_LIBDEFLATE
        use_consumer_libdeflate = 1;
#endif
        if(!out_is_zip_ || use_consumer_libdeflate == 0) {
            int out_len1 = 0;
            for(int i = 0; i < 64; i++) {
                out_len1 += out_lens1[i];
            }

            int now_size1 = out_len1;
            int now_sizes1[comm_size];
            now_sizes1[0] = now_size1;
            MPI_Barrier(MPI_COMM_WORLD);
            if(my_rank) {
                MPI_Send(&now_size1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            } else {
                for(int ii = 1; ii < comm_size; ii++) {
                    MPI_Recv(&(now_sizes1[ii]), 1, MPI_INT, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            if(my_rank == 0) {
                for(int ii = 1; ii < comm_size; ii++) {
                    MPI_Send(now_sizes1, comm_size, MPI_INT, ii, 0, MPI_COMM_WORLD);
                }
            } else {
                MPI_Recv(now_sizes1, comm_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            int pre_sizes1[comm_size];
            pre_sizes1[0] = 0;
            for(int ii = 1; ii < comm_size; ii++) {
                pre_sizes1[ii] = pre_sizes1[ii - 1] + now_sizes1[ii - 1];
            }

            now_pos_base1 = now_pos1_ + pre_sizes1[my_rank];

            for(int ii = 0; ii < comm_size; ii++) {
                now_pos1_ += now_sizes1[ii];
            }


            int out_len2 = 0;
            for(int i = 0; i < 64; i++) {
                out_len2 += out_lens2[i];
            }

            int now_size2 = out_len2;
            int now_sizes2[comm_size];
            now_sizes2[0] = now_size2;
            MPI_Barrier(MPI_COMM_WORLD);
            if(my_rank) {
                MPI_Send(&now_size2, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            } else {
                for(int ii = 1; ii < comm_size; ii++) {
                    MPI_Recv(&(now_sizes2[ii]), 1, MPI_INT, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            if(my_rank == 0) {
                for(int ii = 1; ii < comm_size; ii++) {
                    MPI_Send(now_sizes2, comm_size, MPI_INT, ii, 0, MPI_COMM_WORLD);
                }
            } else {
                MPI_Recv(now_sizes2, comm_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            int pre_sizes2[comm_size];
            pre_sizes2[0] = 0;
            for(int ii = 1; ii < comm_size; ii++) {
                pre_sizes2[ii] = pre_sizes2[ii - 1] + now_sizes2[ii - 1];
            }

            now_pos_base2 = now_pos2_ + pre_sizes2[my_rank];

            for(int ii = 0; ii < comm_size; ii++) {
                now_pos2_ += now_sizes2[ii];
            }
        }
        writeInfosPE.clear();
        for(int i = 0; i < 64; i++) {
            writeInfosPE.push_back({out_data1[i], out_data2[i], out_lens1[i], out_lens2[i], now_pos_base1, now_pos_base2});
        }
    }
    for(int i = 0; i < 64; i++) {
        pe_pre_out_len_slave1[i] = pe_all_data.out_len_slave1[i];
        pe_pre_out_len_slave2[i] = pe_all_data.out_len_slave2[i];
    }
    tsum4 += GetTime() - tt0;

    tt0 = GetTime();
    int isNull = 1;
    int all_stop = 1;
    {
        for(int i = 0; i < 64; i++) {
            data1[i].resize(counts[i]);
            data2[i].resize(counts[i]);
            if(cmd_info_->write_data_) {
                pass_data1[i].resize(data1[i].size());
                pass_data2[i].resize(data2[i].size());
            }
            if(data1[i].size()) isNull = 0;
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
        vector<rabbit::fq::FastqDataPairChunk *> tmp_chunks;
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
        if(pre_fqdatachunks[i] != NULL) {
            fastq_data_pool->Release(pre_fqdatachunks[i]->left_part);
            fastq_data_pool->Release(pre_fqdatachunks[i]->right_part);
        }
    }
    tsum7 += GetTime() - tt0;
}

inline void gather_and_sort_vectors_pe(int rank, int comm_size, int out_round, vector<pair<int, size_t>>& local_vec, int root) {
    int local_size = local_vec.size();
    vector<int> send_data(local_size * 2);
    for (int i = 0; i < local_size; ++i) {
        send_data[2 * i] = local_vec[i].first;
        send_data[2 * i + 1] = static_cast<int>(local_vec[i].second);  // Assuming size_t can be safely cast to int
    }

    // Root process initializes a vector to hold all received data
    vector<int> recv_data;
    if (rank == root) {
        recv_data.reserve(local_size * 2 * comm_size);  // Reserve enough space
    }

    // Root process receives data from all other processes
    if (rank == root) {
        vector<int> buffer(local_size * 2);
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
        vector<pair<int, size_t>> all_data;
        for (int i = 0; i < recv_data.size() / 2; ++i) {
            all_data.emplace_back(recv_data[2 * i], static_cast<size_t>(recv_data[2 * i + 1]));
        }
        sort(all_data.begin(), all_data.end());

        //// Print sorted data for verification
        //for (const auto& p : all_data) {
        //    cout << "(" << p.first << ", " << p.second << ") ";
        //}
        //cout << endl;

        // Optionally clear and repopulate local_vec
        local_vec.clear();
        for (const auto& item : all_data) {
            local_vec.push_back(item);
        }
    }
}



inline void gather_and_sort_vectors_pe2(int rank, int comm_size, int out_round, vector<pair<int, size_t>>& local_vec, int root, MPI_Comm comm) {

    // Assuming each vector is of the same size, find the size
    int local_size = local_vec.size();
    //int expected_size = out_round * comm_size * 64; 

    //// Check if local vector size matches expected size
    //if (local_size != expected_size) {
    //    if (rank == root) {
    //        cerr << "Error: Vector sizes do not match the expected size on rank " << rank << endl;
    //    }
    //    MPI_Abort(comm, 1);  // Abort if sizes do not match
    //}

    // Each process sends its vector size to the root
    vector<int> sizes(comm_size, 0);
    MPI_Gather(&local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, root, comm);

    // Calculate displacements for receiving the data
    vector<int> displs(comm_size, 0);
    for (int i = 1; i < comm_size; ++i) {
        displs[i] = displs[i - 1] + sizes[i - 1];
    }

    // Flatten the vector for sending
    vector<int> send_data(local_size * 2);
    for (int i = 0; i < local_size; ++i) {
        send_data[2 * i] = local_vec[i].first;
        send_data[2 * i + 1] = local_vec[i].second;
    }

    // Root process prepares to receive flattened data
    vector<int> recv_data;
    if (rank == root) {
        recv_data.resize(2 * displs[comm_size - 1] + 2 * sizes[comm_size - 1]);
    }

    // Gather the flattened data at the root
    MPI_Gatherv(send_data.data(), 2 * local_size, MPI_INT,
                recv_data.data(), sizes.data(), displs.data(), MPI_INT, root, comm);

    // Root process unflattens and sorts the data
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


void PeQc::ConsumerPeFastqTask64(ThreadInfo **thread_infos, rabbit::fq::FastqDataPool *fastqPool) {
    rabbit::int64 id = 0;
    rabbit::fq::FastqDataPairChunk *fqdatachunk;
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

    vector<rabbit::fq::FastqDataPairChunk *> fqdatachunks;
    vector<rabbit::fq::FastqDataPairChunk *> pre_fqdatachunks;
    int out_round = 0;
    vector<neoReference> data1[64];
    vector<neoReference> data2[64];
    vector<neoReference> pass_data1[64];
    vector<neoReference> pass_data2[64];
    vector<neoReference> pre_pass_data1[64];
    vector<neoReference> pre_pass_data2[64];
    for(int i = 0; i < 64; i++) {
        data1[i].clear();
        data2[i].clear();
        pass_data1[i].clear();
        pass_data2[i].clear();
        pre_pass_data1[i].clear();
        pre_pass_data2[i].clear();
        pre_fqdatachunks.push_back(NULL);
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

            //left chunk
            {
                Para degz_paras[64];
                size_t out_size[64] = {0};
                for (int i = 0; i < 64; i++) {
                    if(fqdatachunks[i] == NULL) {
                        degz_paras[i].in_buffer = NULL;
                        degz_paras[i].in_size = 0;
                        continue;
                    }
                    memcpy(cc_gz_in_buffer[i], fqdatachunks[i]->left_part->data.Pointer(), fqdatachunks[i]->left_part->size);
                    degz_paras[i].in_buffer = cc_gz_in_buffer[i];
                    degz_paras[i].out_buffer = (char*)fqdatachunks[i]->left_part->data.Pointer();
                    degz_paras[i].in_size = fqdatachunks[i]->left_part->size;
                    degz_paras[i].out_size = &(out_size[i]);
                }


                {
                    lock_guard <mutex> guard(globalMutex);
                    __real_athread_spawn((void *) slave_decompressfunc, degz_paras, 1);
                    athread_join();
                }


                for (int i = 0; i < 64; i++) {
                    if(fqdatachunks[i] == NULL) {
                    } else {
                        if(out_size[i]) fqdatachunks[i]->left_part->size = out_size[i] - 1;
                        else fqdatachunks[i]->left_part->size = out_size[i];
                    }
                }

            }

            //right chunk
            {
                Para degz_paras[64];
                size_t out_size[64] = {0};
                for (int i = 0; i < 64; i++) {
                    if(fqdatachunks[i] == NULL) {
                        degz_paras[i].in_buffer = NULL;
                        degz_paras[i].in_size = 0;
                        continue;
                    }
                    memcpy(cc_gz_in_buffer[i], fqdatachunks[i]->right_part->data.Pointer(), fqdatachunks[i]->right_part->size);
                    degz_paras[i].in_buffer = cc_gz_in_buffer[i];
                    degz_paras[i].out_buffer = (char*)fqdatachunks[i]->right_part->data.Pointer();
                    degz_paras[i].in_size = fqdatachunks[i]->right_part->size;
                    degz_paras[i].out_size = &(out_size[i]);
                }


                {
                    lock_guard <mutex> guard(globalMutex);
                    __real_athread_spawn((void *) slave_decompressfunc, degz_paras, 1);
                    athread_join();
                }


                for (int i = 0; i < 64; i++) {
                    if(fqdatachunks[i] == NULL) {
                    } else {
                        if(out_size[i]) fqdatachunks[i]->right_part->size = out_size[i] - 1;
                        else fqdatachunks[i]->right_part->size = out_size[i];
                    }
                }

            }

        }
#endif
#endif

        t_decom += GetTime() - tt0;

        tt0 = GetTime();
        ProcessFormatQCWrite(allIsNull, data1, data2, pass_data1, pass_data2, pre_pass_data1, pre_pass_data2, fqdatachunks, pre_fqdatachunks, &para, fastqPool);
        t_ngsfunc += GetTime() - tt0;

        //fprintf(stderr, "rank%d pending write size %d\n", my_rank, writeInfosPE.size());
#ifdef CONSUMER_USE_LIBDEFLATE
        tt0 = GetTime();
        if (cmd_info_->write_data_ && out_is_zip_) {
            Para paras[64];
            size_t out_size[64] = {0};
            for(int i = 0; i < writeInfosPE.size(); i++) {
                paras[i].in_buffer = writeInfosPE[i].buffer;
                paras[i].out_buffer = new char[BLOCK_SIZE];
                paras[i].in_size = writeInfosPE[i].buffer_len;
                paras[i].out_size = &(out_size[i]);
                paras[i].level = 1;
            }
            for(int i = writeInfosPE.size(); i < 64; i++) {
                paras[i].in_size = 0;
            }
            double tt00 = GetTime();
            {
                lock_guard<mutex> guard(globalMutex);
                __real_athread_spawn((void *)slave_compressfunc, paras, 1);
                athread_join();
            }
            t_slave_gz2 += GetTime() - tt00;
            for(int i = 0; i < 64; i++) {
                out_gz_block_sizes1.push_back(make_pair(out_round * comm_size * 64 + my_rank * 64 + i, out_size[i]));
            }

            //TODO
//            for(int i = 0; i < 64; i++) {
//                off_idx1 << out_size[i] << endl;
//            }

            int sum_buffer_len1 = 0;
            for(int i = 0; i < writeInfosPE.size(); i++) {
                delete[] writeInfosPE[i].buffer;
                writeInfosPE[i].buffer = paras[i].out_buffer;
                writeInfosPE[i].buffer_len = out_size[i];
                sum_buffer_len1 += out_size[i];
            }

            int now_size = sum_buffer_len1;
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
            long long now_pos_base1 = zip_now_pos1_ + pre_sizes[my_rank];
            MPI_Barrier(MPI_COMM_WORLD);
            for(int ii = 0; ii < comm_size; ii++) {
                zip_now_pos1_ += now_sizes[ii];
            }
            for(int i = 0; i < writeInfosPE.size(); i++) {
                writeInfosPE[i].file_offset = now_pos_base1;
            }

            Para paras2[64];
            size_t out_size2[64] = {0};
            for(int i = 0; i < writeInfosPE.size(); i++) {
                paras2[i].in_buffer = writeInfosPE[i].buffer2;
                paras2[i].out_buffer = new char[BLOCK_SIZE];
                paras2[i].in_size = writeInfosPE[i].buffer_len2;
                paras2[i].out_size = &(out_size2[i]);
                paras2[i].level = 1;
            }
            for(int i = writeInfosPE.size(); i < 64; i++) {
                paras2[i].in_size = 0;
            }
            tt00 = GetTime();
            {
                lock_guard<mutex> guard(globalMutex);
                __real_athread_spawn((void *)slave_compressfunc, paras2, 1);
                athread_join();
            }
            t_slave_gz2 += GetTime() - tt00;
            for(int i = 0; i < 64; i++) {
                out_gz_block_sizes2.push_back(make_pair(out_round * comm_size * 64 + my_rank * 64 + i, out_size2[i]));
            }
            out_round++;


            //TODO
//            for(int i = 0; i < 64; i++) {
//                off_idx2 << out_size2[i] << endl;
//            }
            int sum_buffer_len2 = 0;
            for(int i = 0; i < writeInfosPE.size(); i++) {
                delete[] writeInfosPE[i].buffer2;
                writeInfosPE[i].buffer2 = paras2[i].out_buffer;
                writeInfosPE[i].buffer_len2 = out_size2[i];
                sum_buffer_len2 += out_size2[i];
            }

            now_size = sum_buffer_len2;
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
            pre_sizes[0] = 0;
            for(int ii = 1; ii < comm_size; ii++) {
                pre_sizes[ii] = pre_sizes[ii - 1] + now_sizes[ii - 1];
            }
            long long now_pos_base2 = zip_now_pos2_ + pre_sizes[my_rank];
            MPI_Barrier(MPI_COMM_WORLD);
            for(int ii = 0; ii < comm_size; ii++) {
                zip_now_pos2_ += now_sizes[ii];
            }
            for(int i = 0; i < writeInfosPE.size(); i++) {
                writeInfosPE[i].file_offset2 = now_pos_base2;
            }
        }
        
        t_slave_gz += GetTime() - tt0;
#endif

        if (cmd_info_->write_data_) {
            tt0 = GetTime();
            assert(writeInfosPE.size() == 64);
            long long real_pre_chunk_pos1 = writeInfosPE[0].file_offset;
            long long real_pre_chunk_pos2 = writeInfosPE[0].file_offset2;
            for(int i = 0; i < 64; i++) {
                //mylock.lock();
                while (queueNumNow >= queueSizeLim) {
#ifdef Verbose
                    //printf("waiting to push a chunk to out queue1\n");
#endif
                    usleep(1000);
                }
                out_queue_[queueP2++] = make_pair(make_pair(writeInfosPE[i].buffer, make_pair(writeInfosPE[i].buffer_len, real_pre_chunk_pos1)), make_pair(writeInfosPE[i].buffer2, make_pair(writeInfosPE[i].buffer_len2, real_pre_chunk_pos2)));
                //fprintf(stderr, "out info : %d %d\n", writeInfosPE[i].buffer_len, real_pre_chunk_pos1);
                //out_queue_[queueP2++] = make_pair(Qitem1, Qitem2);
                queueNumNow++;
                real_pre_chunk_pos1 += writeInfosPE[i].buffer_len;
                real_pre_chunk_pos2 += writeInfosPE[i].buffer_len2;
                //mylock.unlock();
            }

            t_push_q += GetTime() - tt0;
        }

        pre_fqdatachunks.clear();
        for(int i = 0; i < 64; i++) {
            pre_pass_data1[i].clear();
            pre_pass_data2[i].clear();
            pre_pass_data1[i] = pass_data1[i];
            pre_pass_data2[i] = pass_data2[i];
            pre_fqdatachunks.push_back(fqdatachunks[i]);
        }
        if(allIsNull) break;

    }


    double tt0 = GetTime();
    gather_and_sort_vectors_pe(my_rank, comm_size, out_round, out_gz_block_sizes1, 0);
    gather_and_sort_vectors_pe(my_rank, comm_size, out_round, out_gz_block_sizes2, 0);

    printf("comm gz block size cost %lf\n", GetTime() - tt0);


    fprintf(stderr, "consumer%d NGSnew tot cost %lf\n", my_rank, GetTime() - t0);
    fprintf(stderr, "consumer%d NGSnew wait producer cost %lf\n", my_rank, t_wait_producer);
    fprintf(stderr, "consumer%d NGSnew format cost %lf\n", my_rank, t_format);
    fprintf(stderr, "consumer%d NGSnew resize cost %lf\n", my_rank, t_resize);
    fprintf(stderr, "consumer%d NGSnew ngsfunc cost %lf\n", my_rank, t_ngsfunc);
    fprintf(stderr, "consumer%d NGSnew      pre_qc1 cost %lf\n", my_rank, tsum1);
    fprintf(stderr, "consumer%d NGSnew      pre_qc1 cost %lf\n", my_rank, tsum2);
    fprintf(stderr, "consumer%d NGSnew      slave cost %lf\n", my_rank, tsum3);
    fprintf(stderr, "consumer%d NGSnew      write cost %lf (%lf [%lf %lf %lf], %lf, %lf, %lf)\n", my_rank, tsum4, tsum4_1, tsum4_1_1, tsum4_1_2, tsum4_1_3, tsum4_2, tsum4_3, tsum4_4);
    fprintf(stderr, "consumer%d NGSnew      isnull comm cost %lf\n", my_rank, tsum5);
    fprintf(stderr, "consumer%d NGSnew      push null comm cost %lf\n", my_rank, tsum6);
    fprintf(stderr, "consumer%d NGSnew      release cost %lf\n", my_rank, tsum7);
    fprintf(stderr, "consumer%d NGSnew gz slave cost %lf\n", my_rank, t_slave_gz);
    fprintf(stderr, "consumer%d NGSnew gz slave2 cost %lf\n", my_rank, t_slave_gz2);
    fprintf(stderr, "consumer%d NGSnew push to queue cost %lf\n", my_rank, t_push_q);
    done_thread_number_++;
    athread_halt();
    fprintf(stderr, "consumer %d in func done\n", my_rank);
}


/**
 * @brief a function to write pe data from out_data1 queue to file1
 */
void PeQc::WriteSeFastqTask12() {
#ifdef Verbose
    double t0 = GetTime();
#endif
    int cnt = 0;
    long long tot_size1 = 0;
    long long tot_size2 = 0;
    bool overWhile = 0;
    //pair<QChunkItem, QChunkItem> QitemPair;
    //QChunkItem Qitem1;
    //QChunkItem Qitem2;
    pair<CIPair, CIPair> now0;
    CIPair now1;
    CIPair now2;
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
                printf("writer %d write done\n", my_rank);
                break;
            }

            //fprintf(stderr, "writer %d waiting\n", my_rank);
            usleep(1000);
        }
        if (overWhile) break;
        //QitemPair = out_queue_[queueP1++];
        //Qitem1 = QitemPair.first;
        //Qitem2 = QitemPair.second;
        now0 = out_queue_[queueP1++];
        now1 = now0.first;
        now2 = now0.second;
        queueNumNow--;
        t_wait += GetTime() - tt0;



        if (out_is_zip_) {
            if (cmd_info_->use_pigz_) {
                fprintf(stderr, "use pigz TODO...\n");
                exit(0);
            } else {

                tt0 = GetTime();
                if(use_out_mem) {
                    memcpy(OutMemData1, now1.first, now1.second.first);
                    memcpy(OutMemData2, now2.first, now2.second.first);

                } else {
#ifdef USE_LIBDEFLATE

#ifdef use_mpi_file

                    MPI_File_seek(fh1, now1.second.second, MPI_SEEK_SET);
                    MPI_File_write(fh1, now1.first, now1.second.first, MPI_CHAR, &status1);

                    MPI_File_seek(fh2, now2.second.second, MPI_SEEK_SET);
                    MPI_File_write(fh2, now2.first, now2.second.first, MPI_CHAR, &status2);
#else

                    fseek(out_stream1_, now1.second.second, SEEK_SET);
                    fwrite(now1.first, sizeof(char), now1.second.first, out_stream1_);

                    fseek(out_stream2_, now2.second.second, SEEK_SET);
                    fwrite(now2.first, sizeof(char), now2.second.first, out_stream2_);
#endif


#else
                    int written1 = gzwrite(zip_out_stream1, now1.first, now1.second.first);
                    if (written1 != now1.second.first) {
                        printf("gzwrite error\n");
                        exit(0);
                    }

                    int written2 = gzwrite(zip_out_stream2, now2.first, now2.second.first);
                    if (written2 != now2.second.first) {
                        printf("gzwrite error\n");
                        exit(0);
                    }
#endif
                }

                t_write += GetTime() - tt0;

                tt0 = GetTime();
                delete[] now1.first;
                delete[] now2.first;
                t_del += GetTime() - tt0;
            }
        } else {

            tot_size1 += now1.second.first;
            if(now1.second.first) {
                tt0 = GetTime();
                if(use_out_mem) {
                    memcpy(OutMemData1, now1.first, now1.second.first);
                } else {
#ifdef use_mpi_file
                    //fprintf(stderr, "writer %d write seek\n", my_rank);
                    MPI_File_seek(fh1, now1.second.second, MPI_SEEK_SET);
                    //fprintf(stderr, "writer %d write ww\n", my_rank);
                    MPI_File_write(fh1, now1.first, now1.second.first, MPI_CHAR, &status1);
                    //fprintf(stderr, "writer %d write ww done\n", my_rank);
#else
                    fseek(out_stream1_, now1.second.second, SEEK_SET);
                    //fprintf(stderr, "rank %d write ww %lld\n", my_rank, now.second.first);
                    fwrite(now1.first, sizeof(char), now1.second.first, out_stream1_);
                    //fprintf(stderr, "rank %d write ww done\n", my_rank);
#endif
                }
                t_write += GetTime() - tt0;

                tt0 = GetTime();
                delete[] now1.first;
                t_del += GetTime() - tt0;
                //fprintf(stderr, "writer %d write del done\n", my_rank);
            }


            tot_size2 += now2.second.first;
            if(now2.second.first) {
                tt0 = GetTime();
                if(use_out_mem) {
                    memcpy(OutMemData2, now2.first, now2.second.first);
                } else {

#ifdef use_mpi_file
                    //fprintf(stderr, "writer %d write seek\n", my_rank);
                    MPI_File_seek(fh2, now2.second.second, MPI_SEEK_SET);
                    //fprintf(stderr, "writer %d write ww\n", my_rank);
                    MPI_File_write(fh2, now2.first, now2.second.first, MPI_CHAR, &status2);
                    //fprintf(stderr, "writer %d write ww done\n", my_rank);
#else
                    fseek(out_stream2_, now2.second.second, SEEK_SET);
                    //fprintf(stderr, "rank %d write ww %lld\n", my_rank, now.second.first);
                    fwrite(now2.first, sizeof(char), now2.second.first, out_stream2_);
                    //fprintf(stderr, "rank %d write ww done\n", my_rank);
#endif
                }
                t_write += GetTime() - tt0;

                tt0 = GetTime();
                delete[] now2.first;
                t_del += GetTime() - tt0;
                //fprintf(stderr, "writer %d write del done\n", my_rank);
            }
            //fprintf(stderr, "writer %d %d done\n", my_rank, cnt++);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    //fprintf(stderr, "writer rank%d write donedone\n", my_rank);
    if (out_is_zip_) {
        if (cmd_info_->use_pigz_) {
            fprintf(stderr, "use pigz TODO...\n");
            exit(0);

        } else {
#ifdef USE_LIBDEFLATE
#ifdef use_mpi_file
            MPI_File_close(&fh1);
            MPI_File_close(&fh2);
#else
            fclose(out_stream1_);
            if(my_rank == 0) {
                truncate(cmd_info_->out_file_name1_.c_str(), sizeof(char) * zip_now_pos1_);
            }
            fclose(out_stream2_);
            if(my_rank == 0) {
                truncate(cmd_info_->out_file_name2_.c_str(), sizeof(char) * zip_now_pos2_);
            }
#endif
            //off_idx1.close();
            //off_idx2.close();

#else
            if (zip_out_stream1) {
                gzflush(zip_out_stream1, Z_FINISH);
                gzclose(zip_out_stream1);
                zip_out_stream1 = NULL;
            }
            if (zip_out_stream2) {
                gzflush(zip_out_stream2, Z_FINISH);
                gzclose(zip_out_stream2);
                zip_out_stream2 = NULL;
            }
#endif
        }

    } else {
#ifdef use_mpi_file
        MPI_File_close(&fh1);
        MPI_File_close(&fh2);
#else
        fclose(out_stream1_);
        if(my_rank == 0) {
            truncate(cmd_info_->out_file_name1_.c_str(), sizeof(char) * now_pos1_);
        }
        fclose(out_stream2_);
        if(my_rank == 0) {
            truncate(cmd_info_->out_file_name2_.c_str(), sizeof(char) * now_pos2_);
        }
#endif
    }

    if(my_rank == 0 && out_is_zip_) {
        double tt0 = GetTime();
        ofstream ofs(cmd_info_->out_file_name1_, ios::binary | ios::app);
        for (const auto& pair : out_gz_block_sizes1) {
            ofs.write(reinterpret_cast<const char*>(&pair.second), sizeof(size_t));
        }
        size_t vector_size = out_gz_block_sizes1.size();
        ofs.write(reinterpret_cast<const char*>(&vector_size), sizeof(size_t));
        ofs.close();

        ofstream ofs2(cmd_info_->out_file_name2_, ios::binary | ios::app);
        for (const auto& pair : out_gz_block_sizes2) {
            ofs2.write(reinterpret_cast<const char*>(&pair.second), sizeof(size_t));
        }
        size_t vector_size2 = out_gz_block_sizes2.size();
        ofs2.write(reinterpret_cast<const char*>(&vector_size2), sizeof(size_t));
        ofs2.close();
        printf("writer final cost %lf\n", GetTime() - tt0);
    }

#ifdef Verbose
    printf("writer wait queue cost %.4f\n", t_wait);
    printf("writer gz slave cost %.4f\n", t_gz_slave);
    printf("writer write cost %.4f\n", t_write);
    printf("writer del cost %.4f\n", t_del);
    printf("writer free cost %.4f\n", t_free);
    printf("writer %d cost %lf, tot size %lld %lld\n", my_rank, GetTime() - t0, tot_size1, tot_size2);
    fprintf(stderr, "writer %d in func done\n", my_rank);
    //printf("write %d cost %.5f --- %lld\n", my_rank, GetTime() - t0, tot_size);
    //
#endif
}



///**
// * @brief a function to write pe data from out_data1 queue to file1
// */
//void PeQc::WriteSeFastqTask1() {
//#ifdef Verbose
//    double t0 = GetTime();
//#endif
//    int cnt = 0;
//    long long tot_size = 0;
//    bool overWhile = 0;
//    CIPair now;
//    while (true) {
//        while (queueNumNow1 == 0) {
//            if (done_thread_number_ == cmd_info_->thread_number_) {
//                overWhile = 1;
//                printf("rank %d write done\n", my_rank);
//                break;
//            }
//
//            //fprintf(stderr, "rank %d write wait\n", my_rank);
//            usleep(1000);
//        }
//        if (overWhile) break;
//        now = out_queue1_[queue1P1++];
//        queueNumNow1--;
//        //fprintf(stderr, "rank %d write get %p\n", my_rank, now.first);
//        if (out_is_zip_) {
//            if (cmd_info_->use_pigz_) {
//                fprintf(stderr, "use pigz TODO...\n");
//                exit(0);
//            } else {
//                int written = gzwrite(zip_out_stream1, now.first, now.second.first);
//                if (written != now.second.first) {
//                    printf("gzwrite error\n");
//                    exit(0);
//                }
//                delete[] now.first;
//            }
//        } else {
//            tot_size += now.second.first;
//            if(now.second.first) {
//
//#ifdef use_mpi_file
//                //fprintf(stderr, "rank %d write seek %lld\n", my_rank, now.second.first);
//                MPI_File_seek(fh1, now.second.second, MPI_SEEK_SET);
//                //fprintf(stderr, "rank %d write ww %lld\n", my_rank, now.second.first);
//                MPI_File_write(fh1, now.first, now.second.first, MPI_CHAR, &status1);
//                //fprintf(stderr, "rank %d write ww done\n", my_rank);
//#else
//                fseek(out_stream1_, now.second.second, SEEK_SET);
//                //fprintf(stderr, "rank %d write ww %lld\n", my_rank, now.second.first);
//                fwrite(now.first, sizeof(char), now.second.first, out_stream1_);
//                //fprintf(stderr, "rank %d write ww done\n", my_rank);
//#endif
//                delete[] now.first;
//                //fprintf(stderr, "rank %d write del done\n", my_rank);
//            }
//            //fprintf(stderr, "write%d %d done\n", my_rank, cnt++);
//        }
//    }
//    //fprintf(stderr, "111111 rank%d donedone1\n", my_rank);
//    MPI_Barrier(MPI_COMM_WORLD);
//    //fprintf(stderr, "111111 rank%d donedone2\n", my_rank);
//    if (out_is_zip_) {
//        if (cmd_info_->use_pigz_) {
//            fprintf(stderr, "use pigz TODO...\n");
//            exit(0);
//        } else {
//            if (zip_out_stream1) {
//                gzflush(zip_out_stream1, Z_FINISH);
//                gzclose(zip_out_stream1);
//                zip_out_stream1 = NULL;
//            }
//        }
//
//    } else {
//#ifdef use_mpi_file
//        //fprintf(stderr, "111111 rank%d donedone3\n", my_rank);
//        MPI_File_close(&fh1);
//        //fprintf(stderr, "111111 rank%d donedone4\n", my_rank);
//#else
//        fclose(out_stream1_);
//        if(my_rank == 0) {
//            truncate(cmd_info_->out_file_name1_.c_str(), sizeof(char) * now_pos1_);
//        }
//#endif
//    }
//#ifdef Verbose
//    cerr << "write" << my_rank << " cost " << GetTime() - t0 << ", tot size " << tot_size << endl;
//    //printf("write %d cost %.5f --- %lld\n", my_rank, GetTime() - t0, tot_size);
//    //
//#endif
//}
//
///**
// * @brief a function to write pe data from out_data2 queue to file2
// */
//void PeQc::WriteSeFastqTask2() {
//#ifdef Verbose
//    double t0 = GetTime();
//#endif
//    int cnt = 0;
//    long long tot_size = 0;
//    bool overWhile = 0;
//    CIPair now;
//    while (true) {
//        while (queueNumNow2 == 0) {
//            if (done_thread_number_ == cmd_info_->thread_number_) {
//                overWhile = 1;
//                printf("rank %d write done\n", my_rank);
//                break;
//            }
//
//            //fprintf(stderr, "rank %d write wait\n", my_rank);
//            usleep(1000);
//        }
//        if (overWhile) break;
//        now = out_queue2_[queue2P1++];
//        queueNumNow2--;
//        //fprintf(stderr, "rank %d write get %p\n", my_rank, now.first);
//        if (out_is_zip_) {
//            if (cmd_info_->use_pigz_) {
//                fprintf(stderr, "use pigz TODO...\n");
//                exit(0);
//            } else {
//                int written = gzwrite(zip_out_stream2, now.first, now.second.first);
//                if (written != now.second.first) {
//                    printf("gzwrite error\n");
//                    exit(0);
//                }
//                delete[] now.first;
//            }
//        } else {
//            tot_size += now.second.first;
//            if(now.second.first) {
//
//#ifdef use_mpi_file
//                //fprintf(stderr, "rank %d write seek %lld\n", my_rank, now.second.first);
//                MPI_File_seek(fh2, now.second.second, MPI_SEEK_SET);
//                //fprintf(stderr, "rank %d write ww %lld\n", my_rank, now.second.first);
//                MPI_File_write(fh2, now.first, now.second.first, MPI_CHAR, &status2);
//                //fprintf(stderr, "rank %d write ww done\n", my_rank);
//#else
//                fseek(out_stream2_, now.second.second, SEEK_SET);
//                //fprintf(stderr, "rank %d write ww %lld\n", my_rank, now.second.first);
//                fwrite(now.first, sizeof(char), now.second.first, out_stream2_);
//                //fprintf(stderr, "rank %d write ww done\n", my_rank);
//#endif
//                delete[] now.first;
//                //fprintf(stderr, "rank %d write del done\n", my_rank);
//            }
//            //fprintf(stderr, "write%d %d done\n", my_rank, cnt++);
//        }
//    }
//    //fprintf(stderr, "222222 rank%d donedone1\n", my_rank);
//    MPI_Barrier(MPI_COMM_WORLD);
//    //fprintf(stderr, "222222 rank%d donedone2\n", my_rank);
//    if (out_is_zip_) {
//        if (cmd_info_->use_pigz_) {
//            fprintf(stderr, "use pigz TODO...\n");
//            exit(0);
//        } else {
//            if (zip_out_stream2) {
//                gzflush(zip_out_stream2, Z_FINISH);
//                gzclose(zip_out_stream2);
//                zip_out_stream2 = NULL;
//            }
//        }
//
//    } else {
//#ifdef use_mpi_file
//        //fprintf(stderr, "222222 rank%d donedone3\n", my_rank);
//        MPI_File_close(&fh2);
//        //fprintf(stderr, "222222 rank%d donedone4\n", my_rank);
//#else
//        fclose(out_stream2_);
//        if(my_rank == 0) {
//            truncate(cmd_info_->out_file_name2_.c_str(), sizeof(char) * now_pos2_);
//        }
//#endif
//    }
//#ifdef Verbose
//    cerr << "write" << my_rank << " cost " << GetTime() - t0 << ", tot size " << tot_size << endl;
//    //printf("write %d cost %.5f --- %lld\n", my_rank, GetTime() - t0, tot_size);
//    //
//#endif
//}


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


void PeQc::ProcessPeFastq() {
    //if(my_rank == 0) block_size = 6 * (1 << 20);
    //else block_size = 4 * (1 << 20);
    auto *fastqPool = new rabbit::fq::FastqDataPool(Q_lim_pe * 64 * 2, BLOCK_SIZE);
    //rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> queue1(Q_lim_pe, 1);
    auto **p_thread_info = new ThreadInfo *[slave_num];
    for (int t = 0; t < slave_num; t++) {
        p_thread_info[t] = new ThreadInfo(cmd_info_, true);
    }
    rabbit::uint32 tmpSize = SWAP1_SIZE;
    if (cmd_info_->seq_len_ <= 200) tmpSize = SWAP2_SIZE;
    
    if(in_is_zip_) {
        fqFileReader = new rabbit::fq::FastqFileReader(cmd_info_->in_file_name1_, *fastqPool, cmd_info_->in_file_name2_, in_is_zip_, tmpSize, start_line_, end_line_, use_in_mem);
    } else {
        fqFileReader = new rabbit::fq::FastqFileReader(cmd_info_->in_file_name1_, *fastqPool, cmd_info_->in_file_name2_, in_is_zip_, tmpSize, start_pos_, end_pos_, use_in_mem);
    }
    if(use_in_mem) {
        MPI_Barrier(MPI_COMM_WORLD);
    }
    double ttt = GetTime();


    //tagg

    thread producer(bind(&PeQc::ProducerPeFastqTask64, this, cmd_info_->in_file_name1_, cmd_info_->in_file_name2_, fastqPool));
    thread consumer(bind(&PeQc::ConsumerPeFastqTask64, this, p_thread_info, fastqPool));

    if (cmd_info_->write_data_) {
        WriteSeFastqTask12();
        printf("write%d done\n", my_rank);
    }


    consumer.join();
    printf("consumer%d done\n", my_rank);

    producer.join();
    printf("producer%d done\n", my_rank);


    printf("rank%d all pro done1\n", my_rank);
    MPI_Barrier(MPI_COMM_WORLD);
    printf("rank%d all pro done2\n", my_rank);


#ifdef Verbose
    printf("all thread done\n");
    printf("now merge thread info\n");
#endif

    printf("TOT TIME1 %lf\n", GetTime() - ttt);

    double tt00 = GetTime();
    vector<State *> pre_vec_state1;
    vector<State *> pre_vec_state2;
    vector<State *> aft_vec_state1;
    vector<State *> aft_vec_state2;

    for (int t = 0; t < slave_num; t++) {
        pre_vec_state1.push_back(p_thread_info[t]->pre_state1_);
        pre_vec_state2.push_back(p_thread_info[t]->pre_state2_);
        aft_vec_state1.push_back(p_thread_info[t]->aft_state1_);
        aft_vec_state2.push_back(p_thread_info[t]->aft_state2_);
    }
    State* pre_state_tmp1;
    State* pre_state_tmp2;
    State* aft_state_tmp1;
    State* aft_state_tmp2;
    
    if(cmd_info_->do_overrepresentation_) {
    //if(0) {
        pre_state_tmp1 = State::MergeStatesSlave(pre_vec_state1);
        pre_state_tmp2 = State::MergeStatesSlave(pre_vec_state2);
        aft_state_tmp1 = State::MergeStatesSlave(aft_vec_state1);
        aft_state_tmp2 = State::MergeStatesSlave(aft_vec_state2);
    } else {
        pre_state_tmp1 = State::MergeStates(pre_vec_state1);
        pre_state_tmp2 = State::MergeStates(pre_vec_state2);
        aft_state_tmp1 = State::MergeStates(aft_vec_state1);
        aft_state_tmp2 = State::MergeStates(aft_vec_state2);
    }
    printf("merge1 done\n");
    printf("merge cost %lf\n", GetTime() - tt00);

    tt00 = GetTime();
    vector<State *> pre_state_mpis1;
    pre_state_mpis1.push_back(pre_state_tmp1);
    vector<State *> aft_state_mpis1;
    aft_state_mpis1.push_back(aft_state_tmp1);
    vector<State *> pre_state_mpis2;
    pre_state_mpis2.push_back(pre_state_tmp2);
    vector<State *> aft_state_mpis2;
    aft_state_mpis2.push_back(aft_state_tmp2);
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
            pre_state_mpis1.push_back(tmp_state);
        }
    } else {
        string pre_state_is = pre_state_tmp1->ParseString();
        State* tmp_state = new State(pre_state_is.c_str(), pre_state_is.length(), cmd_info_, cmd_info_->seq_len_, cmd_info_->qul_range_, false);
        //bool res = checkStates(pre_state_tmp1, tmp_state);
        //if(!res) cerr << "GGGGGGGG" << endl;
        int now_size = pre_state_is.length();
        MPI_Send(&now_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(pre_state_is.c_str(), now_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0) {
        for(int i = 1; i < comm_size; i++) {
            int now_size = 0;
            MPI_Recv(&now_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            char* info = new char[now_size];
            MPI_Recv(info, now_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            State* tmp_state = new State(info, now_size, cmd_info_, cmd_info_->seq_len_, cmd_info_->qul_range_, false);
            delete[] info;
            pre_state_mpis2.push_back(tmp_state);
        }
    } else {
        string pre_state_is = pre_state_tmp2->ParseString();
        State* tmp_state = new State(pre_state_is.c_str(), pre_state_is.length(), cmd_info_, cmd_info_->seq_len_, cmd_info_->qul_range_, false);
        //bool res = checkStates(pre_state_tmp2, tmp_state);
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
            aft_state_mpis1.push_back(tmp_state);
        }
    } else {
        string aft_state_is = aft_state_tmp1->ParseString();
        int now_size = aft_state_is.length();
        MPI_Send(&now_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(aft_state_is.c_str(), now_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank == 0) {
        for(int i = 1; i < comm_size; i++) {
            int now_size = 0;
            MPI_Recv(&now_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            char* info = new char[now_size];
            MPI_Recv(info, now_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            State* tmp_state = new State(info, now_size, cmd_info_, cmd_info_->seq_len_, cmd_info_->qul_range_, false);
            delete[] info;
            aft_state_mpis2.push_back(tmp_state);
        }
    } else {
        string aft_state_is = aft_state_tmp2->ParseString();
        int now_size = aft_state_is.length();
        MPI_Send(&now_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(aft_state_is.c_str(), now_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    printf("mpi comm cost %lf\n", GetTime() - tt00);
    printf("aft merge done\n");

    tt00 = GetTime();
    if(my_rank) {
        delete pre_state_tmp1;
        delete pre_state_tmp2;
        delete aft_state_tmp1;
        delete aft_state_tmp2;
        delete fastqPool;
        for (int t = 0; t < slave_num; t++) {
            delete p_thread_info[t];
        }
        delete[] p_thread_info;
        delete fqFileReader;
        if(use_out_mem) {
            delete[] OutMemData1;
            delete[] OutMemData2;
        }
        return;
    }



    auto pre_state1 = State::MergeStates(pre_state_mpis1);
    auto pre_state2 = State::MergeStates(pre_state_mpis2);
    auto aft_state1 = State::MergeStates(aft_state_mpis1);
    auto aft_state2 = State::MergeStates(aft_state_mpis2);

    printf("merge2 cost %lf\n", GetTime() - tt00);

    tt00 = GetTime();
#ifdef Verbose
    if (cmd_info_->do_overrepresentation_) {
        printf("orp cost %f\n", pre_state1->GetOrpCost() + pre_state2->GetOrpCost() + aft_state1->GetOrpCost() + aft_state2->GetOrpCost());
    }
    printf("merge done\n");
#endif
    printf("\nprint read1 (before filter) info :\n");
    State::PrintStates(pre_state1);
    printf("\nprint read1 (after filter) info :\n");
    State::PrintStates(aft_state1);
    printf("\n");

    printf("\nprint read2 (before filter) info :\n");
    State::PrintStates(pre_state2);
    printf("\nprint read2 (after filter) info :\n");
    State::PrintStates(aft_state2);
    printf("\n");
    if (cmd_info_->print_what_trimmed_) {
        State::PrintAdapterToFile(aft_state1);
        State::PrintAdapterToFile(aft_state2);
    }
    State::PrintFilterResults(aft_state1);
    printf("\n");

    if (cmd_info_->do_overrepresentation_ && cmd_info_->print_ORP_seqs_) {

        auto pre_hash_graph1 = pre_state1->GetHashGraph();
        int pre_hash_num1 = pre_state1->GetHashNum();
        auto pre_hash_graph2 = pre_state2->GetHashGraph();
        int pre_hash_num2 = pre_state2->GetHashNum();

        auto aft_hash_graph1 = aft_state1->GetHashGraph();
        int aft_hash_num1 = aft_state1->GetHashNum();
        auto aft_hash_graph2 = aft_state2->GetHashGraph();
        int aft_hash_num2 = aft_state2->GetHashNum();


        int spg = cmd_info_->overrepresentation_sampling_;
        ofstream ofs;

        string srr_name1 = cmd_info_->in_file_name1_;
        srr_name1 = PaseFileName(srr_name1);

        string srr_name2 = cmd_info_->in_file_name2_;
        srr_name2 = PaseFileName(srr_name2);

        string out_name1 = "pe_" + srr_name1 + "_before_ORP_sequences.txt";
        ofs.open(out_name1, ifstream::out);

        int cnt1 = 0;
        ofs << "sequence count"
            << "\n";
        for (int i = 0; i < pre_hash_num1; i++) {
            if (!overRepPassed(pre_hash_graph1[i].seq, pre_hash_graph1[i].cnt, spg)) continue;
            ofs << pre_hash_graph1[i].seq << " " << pre_hash_graph1[i].cnt << "\n";
            cnt1++;
        }
        ofs.close();
        printf("in %s (before filter) find %d possible overrepresented sequences (store in %s)\n",
                srr_name1.c_str(), cnt1, out_name1.c_str());

        string out_name2 = "pe_" + srr_name2 + "_before_ORP_sequences.txt";
        ofs.open(out_name2, ifstream::out);
        int cnt2 = 0;
        ofs << "sequence count"
            << "\n";
        for (int i = 0; i < pre_hash_num2; i++) {
            if (!overRepPassed(pre_hash_graph2[i].seq, pre_hash_graph2[i].cnt, spg)) continue;
            ofs << pre_hash_graph2[i].seq << " " << pre_hash_graph2[i].cnt << "\n";
            cnt2++;
        }
        ofs.close();
        printf("in %s (before filter) find %d possible overrepresented sequences (store in %s)\n",
                srr_name2.c_str(), cnt2, out_name2.c_str());


        out_name1 = "pe_" + srr_name1 + "_after_ORP_sequences.txt";
        ofs.open(out_name1, ifstream::out);
        cnt1 = 0;
        ofs << "sequence count"
            << "\n";
        for (int i = 0; i < aft_hash_num1; i++) {
            if (!overRepPassed(aft_hash_graph1[i].seq, aft_hash_graph1[i].cnt, spg)) continue;
            ofs << aft_hash_graph1[i].seq << " " << aft_hash_graph1[i].cnt << "\n";
            cnt1++;
        }
        ofs.close();
        printf("in %s (after filter) find %d possible overrepresented sequences (store in %s)\n", srr_name1.c_str(),
                cnt1, out_name1.c_str());

        out_name2 = "pe_" + srr_name2 + "_after_ORP_sequences.txt";
        ofs.open(out_name2, ifstream::out);
        cnt2 = 0;
        ofs << "sequence count"
            << "\n";
        for (int i = 0; i < aft_hash_num2; i++) {
            if (!overRepPassed(aft_hash_graph2[i].seq, aft_hash_graph2[i].cnt, spg)) continue;
            ofs << aft_hash_graph2[i].seq << " " << aft_hash_graph2[i].cnt << "\n";
            cnt2++;
        }
        ofs.close();
        printf("in %s (after filter) find %d possible overrepresented sequences (store in %s)\n", srr_name2.c_str(),
                cnt2, out_name2.c_str());

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
        printf("Duplication rate : %.5f %%\n", dupRate * 100.0);
        delete[] dupHist;
        delete[] dupMeanGC;
    }

    int *merge_insert_size;
    if (!cmd_info_->no_insert_size_) {
        merge_insert_size = new int[cmd_info_->max_insert_size_ + 1];
        memset(merge_insert_size, 0, sizeof(int) * (cmd_info_->max_insert_size_ + 1));

        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            for (int i = 0; i <= cmd_info_->max_insert_size_; i++) {
                merge_insert_size[i] += p_thread_info[t]->insert_size_dist_[i];
            }
        }
        int mx_id = 0;
        for (int i = 0; i < cmd_info_->max_insert_size_; i++) {
            if (merge_insert_size[i] > merge_insert_size[mx_id]) mx_id = i;
        }
        //printf("Insert size peak (evaluated by paired-end reads): %d\n", mx_id);
        printf("Insert size peak (based on PE overlap analyze): %d\n", mx_id);
    }
    printf("print and other cost %lf\n", GetTime() - tt00);

    tt00 = GetTime();
    string srr_name1 = cmd_info_->in_file_name1_;
    srr_name1 = PaseFileName(srr_name1);
    string srr_name2 = cmd_info_->in_file_name2_;
    srr_name2 = PaseFileName(srr_name2);
    Repoter::ReportHtmlPe(srr_name1 + "_" + srr_name2 + "_RabbitQCPlus.html", pre_state1, pre_state2, aft_state1,
            aft_state2, cmd_info_->in_file_name1_,
            cmd_info_->in_file_name2_, dupRate * 100.0, merge_insert_size);
#ifdef Verbose
    printf("report done\n");
#endif
    printf("report cost %lf\n", GetTime() - tt00);

    if(use_in_mem) {
        fprintf(stderr, "TOT TIME %lf\n", GetTime() - ttt);
    }
    tt00 = GetTime();
    delete pre_state1;
    delete pre_state2;
    delete aft_state1;
    delete aft_state2;
    if (!cmd_info_->no_insert_size_)
        delete[] merge_insert_size;

    delete fastqPool;
    for (int t = 0; t < slave_num; t++) {
        delete p_thread_info[t];
    }

    delete[] p_thread_info;
    
    delete fqFileReader;
    if(use_out_mem) {
        delete[] OutMemData1;
        delete[] OutMemData2;
    }
    printf("delete cost %lf\n", GetTime() - tt00);

}
