// O
// Created by ylf9811 on 2021/7/6.
//
//

#include "peqc.h"
using namespace std;

int block_size;

extern "C" {
#include <athread.h>
#include <pthread.h>
    void slave_ngspefunc();
    void slave_formatpefunc();
}



int Q_lim = 4;

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
    if (data_[pos_] != '+') std::cout << "core dump is pos: " << pos_ << " char: " << data_[pos_] << std::endl;
    ASSERT(data_[pos_] == '+');// pos0 was the start of tag
    return pos0;
}


//#define use_in_mem
//#define use_out_mem
//
#define use_mpi_file

#ifdef use_out_mem
char* OutMemData1;
char* OutMemData2;
#endif



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

    ifstream gFile;
    gFile.open(cmd_info1->in_file_name1_.c_str());
    gFile.seekg(0, ios_base::end);
    long long real_file_size = gFile.tellg();
    gFile.close();

    gFile.open(cmd_info1->in_file_name2_.c_str());
    gFile.seekg(0, ios_base::end);
    long long real_file_size2 = gFile.tellg();
    gFile.close();

    if(real_file_size != real_file_size2) {
        fprintf(stderr, "PE in files size error %lld %lld\n", real_file_size, real_file_size2);
    }

    int64_t pre_size = ceil(1.0 * real_file_size / (2 * comm_size + 1));
    int64_t start_pos, end_pos;
    if(my_rank == 0) {
        start_pos = 0;
        end_pos = start_pos + pre_size * 3;
    } else {
        start_pos = pre_size * (my_rank * 2 + 1);
        end_pos = start_pos + pre_size * 2;
    }
    if(end_pos > real_file_size) end_pos = real_file_size;

    //cerr << "startpos" << my_rank << " " << start_pos << "-" << end_pos << endl;
    FILE *pre_fp;
    pre_fp = fopen(cmd_info1->in_file_name1_.c_str(), "rb");
    fseek(pre_fp, start_pos, SEEK_SET);
    char *tmp_chunk = new char[1 << 10];
    int res_size = fread(tmp_chunk, sizeof(char), 1 << 10, pre_fp);
    MPI_Barrier(MPI_COMM_WORLD);
    //cerr << "res_size" << my_rank << " " << res_size << endl;
    int64_t right_pos;
    if(my_rank == 0) right_pos = 0;
    else right_pos = GetNextFastq(tmp_chunk, 0, res_size);
    //cerr << "right_pos" << my_rank << " " << right_pos << endl;
    fclose(pre_fp);

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


    p_out_queue_ = new vector<rabbit::fq::FastqDataPairChunk *>[1 << 20];
    p_queueP1 = 0;
    p_queueP2 = 0;
    p_queueNumNow = 0;
    p_queueSizeLim = Q_lim;

 
    if (cmd_info1->write_data_) {

        out_queue_ = new pair<CIPair, CIPair>[1 << 20];
        queueP1 = 0;
        queueP2 = 0;
        queueNumNow = 0;
        queueSizeLim = Q_lim;

        //out_queue1_ = new CIPair[1 << 20];
        //queue1P1 = 0;
        //queue1P2 = 0;
        //queueNumNow1 = 0;
        //queueSizeLim1 = Q_lim;
        if (cmd_info_->interleaved_out_ == 0) {
            //out_queue2_ = new CIPair[1 << 20];
            //queue2P1 = 0;
            //queue2P2 = 0;
            //queueNumNow2 = 0;
            //queueSizeLim2 = Q_lim;
        }
        if (out_is_zip_) {
            if (cmd_info1->use_pigz_) {
                pigzQueueNumNow1 = 0;
                pigzQueueSizeLim1 = 1 << 5;
                string out_name1 = cmd_info1->out_file_name1_;
                out_name1 = out_name1.substr(0, out_name1.find(".gz"));
                if((out_stream1_ = fopen(out_name1.c_str(),"w")) == NULL) {
                    printf("File cannot be opened\n");
                }
                fclose(out_stream1_);

                //out_stream1_.open(out_name1);
                //out_stream1_.close();

                pigzQueueNumNow2 = 0;
                pigzQueueSizeLim2 = 1 << 5;
                string out_name2 = cmd_info1->out_file_name2_;
                out_name2 = out_name2.substr(0, out_name2.find(".gz"));
                if((out_stream2_ = fopen(out_name2.c_str(),"w")) == NULL) {
                    printf("File cannot be opened\n");
                }
                fclose(out_stream2_);

                //out_stream2_.open(out_name2);
                //out_stream2_.close();
#ifdef Verbose
                printf("now use pigz to compress output data\n");
#endif

            } else {
#ifdef Verbose
                printf("open gzip stream1 %s\n", cmd_info1->out_file_name1_.c_str());
                printf("open gzip stream2 %s\n", cmd_info1->out_file_name2_.c_str());
#endif
                zip_out_stream1 = gzopen(cmd_info1->out_file_name1_.c_str(), "w");
                gzsetparams(zip_out_stream1, cmd_info1->compression_level_, Z_DEFAULT_STRATEGY);
                gzbuffer(zip_out_stream1, 1024 * 1024);
                zip_out_stream2 = gzopen(cmd_info1->out_file_name2_.c_str(), "w");
                gzsetparams(zip_out_stream2, cmd_info1->compression_level_, Z_DEFAULT_STRATEGY);
                gzbuffer(zip_out_stream2, 1024 * 1024);
            }
        } else {
#ifdef Verbose
            printf("open stream1 %s\n", cmd_info1->out_file_name1_.c_str());
            if (cmd_info_->interleaved_out_ == 0)
                printf("open stream2 %s\n", cmd_info1->out_file_name2_.c_str());
#endif

#ifdef use_out_mem
            //OutMemData1 = new char[real_file_size];
            OutMemData1 = new char[128 << 20];
            //OutMemData2 = new char[real_file_size];
            OutMemData2 = new char[128 << 20];
#endif
 

#ifdef use_mpi_file

            std::ofstream file1(cmd_info1->out_file_name1_.c_str());

            if (file1.is_open()) {
                file1 << "Hello, file!" << std::endl;
                file1.close();
            } else {
                std::cerr << "Failed to open file1: " << cmd_info1->out_file_name1_ << std::endl;
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

                std::ofstream file2(cmd_info1->out_file_name2_.c_str());

                if (file2.is_open()) {
                    file2 << "Hello, file!" << std::endl;
                    file2.close();
                } else {
                    std::cerr << "Failed to open file1: " << cmd_info1->out_file_name2_ << std::endl;
                }

                int err2 = MPI_File_open(MPI_COMM_WORLD, cmd_info1->out_file_name2_.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh2);

                if (err2 != MPI_SUCCESS) {
                    char error_string[BUFSIZ];
                    int length_of_error_string, error_class;

                    MPI_Error_class(err2, &error_class);
                    MPI_Error_string(error_class, error_string, &length_of_error_string);

                    fprintf(stderr, "%3d: Couldn't open file, MPI Error: %s\n", my_rank, error_string);

                    MPI_Abort(MPI_COMM_WORLD, err1);
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
        pugzQueue1 = new moodycamel::ReaderWriterQueue<pair<char *, int>>(1 << 10);
        pugzQueue2 = new moodycamel::ReaderWriterQueue<pair<char *, int>>(1 << 10);
    }
    if (cmd_info1->use_pigz_) {
        pigzQueue1 = new moodycamel::ReaderWriterQueue<pair<char *, int>>;
        pigzLast1.first = new char[1 << 24];
        pigzLast1.second = 0;
        pigzQueue2 = new moodycamel::ReaderWriterQueue<pair<char *, int>>;
        pigzLast2.first = new char[1 << 24];
        pigzLast2.second = 0;
    }

    pugzDone1 = 0;
    pugzDone2 = 0;
    producerDone = 0;
    writerDone1 = 0;
    writerDone2 = 0;
    consumerCommDone = 0;
    producerStop = 0;
    now_chunks = 0;
    mx_chunks = 0;
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

void PeQc::ProducerPeInterFastqTask(string file, rabbit::fq::FastqDataPool *fastq_data_pool,
        rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {
#ifdef Verbose
    double t0 = GetTime();
#endif
    rabbit::fq::FastqFileReader *fqFileReader;
    rabbit::uint32 tmpSize = 1 << 20;
    if (cmd_info_->seq_len_ <= 200) tmpSize = 1 << 14;
    fqFileReader = new rabbit::fq::FastqFileReader(file, *fastq_data_pool, "", in_is_zip_);
    int64_t n_chunks = 0;
    while (true) {
        rabbit::fq::FastqDataChunk *fqdatachunk;
        fqdatachunk = fqFileReader->readNextInterChunk();
        if (fqdatachunk == NULL) break;
        n_chunks++;
        dq.Push(n_chunks, fqdatachunk);
    }

    dq.SetCompleted();
    delete fqFileReader;
#ifdef Verbose
    cout << "file " << file << " has " << n_chunks << " chunks" << endl;
    printf("producer cost %.3f\n", GetTime() - t0);
#endif
}

rabbit::fq::FastqFileReader *fqFileReader;

void PeQc::ProducerPeFastqTask64(string file, string file2, rabbit::fq::FastqDataPool *fastqPool) {
    double t0 = GetTime();
    rabbit::uint32 tmpSize = 1 << 20;
    if (cmd_info_->seq_len_ <= 200) tmpSize = 1 << 14;
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
#ifdef use_in_mem
        fqdatachunk = fqFileReader->readNextPairChunkFromMem(offset, part_sizes[my_rank]);
#else 
        fqdatachunk = fqFileReader->readNextPairChunk(offset, part_sizes[my_rank]);
#endif
        t_sum1 += GetTime() - tt0;

        tt0 = GetTime();
        if (fqdatachunk == NULL) {
            if(tmp_chunks.size()) {
                //fprintf(stderr, "producer %d put last %d\n", my_rank, tmp_chunks.size());
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
            producerStop = 1;
            //fprintf(stderr, "producer %d stop, tot chunk size %d\n", my_rank, n_chunks);
            if(cmd_info_->write_data_) {
                now_chunks = n_chunks;
                while(consumerCommDone == 0) {
                    usleep(1000);
                    //fprintf(stderr, "producer %d waiting consumerCommDone\n", my_rank);
                }
                int bu_chunks = mx_chunks - now_chunks;
                cerr << "producer bu rank" << my_rank << " : " << bu_chunks << " of (" << now_chunks << ", " << mx_chunks << ")" << endl;
                if(bu_chunks) {
                    for(int i = 0; i < bu_chunks; i++) {
                        if(tmp_chunks.size() == 64) {
                            //p_mylock.lock();
                            while (p_queueNumNow >= p_queueSizeLim) {
                                usleep(1000);
                            }
                            p_out_queue_[p_queueP2++] = tmp_chunks;
                            p_queueNumNow++;
                            //p_mylock.unlock();
                            tmp_chunks.clear();
                            n_chunks++;
                        } else {
                            tmp_chunks.push_back(fqdatachunk); 
                        }
                    }
                    if(tmp_chunks.size()) {
                        //p_mylock.lock();
                        while (p_queueNumNow >= p_queueSizeLim) {
                            usleep(1000);
                        }
                        p_out_queue_[p_queueP2++] = tmp_chunks;
                        p_queueNumNow++;
                        //p_mylock.unlock();
                        tmp_chunks.clear();
                        n_chunks++;
                    }
                }
            }
            break;
        }
        tot_size += fqdatachunk->left_part->size;
        if(fqdatachunk->left_part->size <= 0 || tot_size <= 0) {
            if(tmp_chunks.size()) {
                //p_mylock.lock();
                while (p_queueNumNow >= p_queueSizeLim) {
                    usleep(1000);
                    //fprintf(stderr, "producer %d waiting0\n", my_rank);
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
        //fprintf(stderr, "producer %d read a chunk\n", my_rank);
        if(tmp_chunks.size() == 64) {
            //fprintf(stderr, "producer %d put 64\n", my_rank);
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

    printf("producer sum1 cost %lf\n", t_sum1);
    printf("producer sum2 cost %lf\n", t_sum2);
    printf("producer%d cost %lf\n", my_rank, GetTime() - t0);
    fprintf(stderr, "producer %d in func done\n", my_rank);
    //dq.SetCompleted();
    producerDone = 1;
}



void PeQc::ProducerPeFastqTask(string file, string file2, rabbit::fq::FastqDataPool *fastqPool,
        rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> &dq) {
#ifdef Verbose
    double t0 = GetTime();
#endif
    
    rabbit::uint32 tmpSize = 1 << 20;
    if (cmd_info_->seq_len_ <= 200) tmpSize = 1 << 14;
    int64_t n_chunks = 0;


    if (cmd_info_->use_pugz_) {
        pair<char *, int> last1;
        pair<char *, int> last2;
        last1.first = new char[1 << 20];
        last1.second = 0;
        last2.first = new char[1 << 20];
        last2.second = 0;
        while (true) {
            rabbit::fq::FastqDataPairChunk *fqdatachunk;
            fqdatachunk = fqFileReader->readNextPairChunkParallel(pugzQueue1, pugzQueue2, &pugzDone1, &pugzDone2, last1,
                    last2);
            if (fqdatachunk == NULL) break;
            n_chunks++;
            dq.Push(n_chunks, fqdatachunk);
        }
        delete[] last1.first;
        delete[] last2.first;
    } else {
        bool is_first = 1;
        int64_t tot_size = 0;
        while (true) {
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
     
#ifdef use_in_mem
            fqdatachunk = fqFileReader->readNextPairChunkFromMem(offset, part_sizes[my_rank]);
#else 
            fqdatachunk = fqFileReader->readNextPairChunk(offset, part_sizes[my_rank]);
#endif
            if (fqdatachunk == NULL) {

                //printf("null\n");
                if(cmd_info_->write_data_) {
                //if(0) {
                    //printf("producer%d stop\n", my_rank);
                    //cerr << "producer" << my_rank << "stop" << endl;
                    now_chunks = n_chunks;
                    producerStop = 1;
                    while(consumerCommDone == 0) {
                        usleep(1000);
                    }
                    //printf("producer%d get val done %lld %lld\n", my_rank, now_chunks, mx_chunks);
                    int bu_chunks = mx_chunks - now_chunks;
                    //cerr << "bu rank" << my_rank << " : " << bu_chunks << " of (" << now_chunks << ", " << mx_chunks << ")" << endl;
                    if(bu_chunks) {
                        for(int i = 0; i < bu_chunks; i++) {
                            dq.Push(n_chunks, fqdatachunk);
                            n_chunks++;
                        }
                    }
                }
                //fprintf(stderr, "rank%d producer done === %lld\n", my_rank, n_chunks);
                break;
            }
            tot_size += fqdatachunk->left_part->size;
            if(fqdatachunk->left_part->size <= 0 || tot_size <= 0) break;
            n_chunks++;
            dq.Push(n_chunks, fqdatachunk);
            //printf("producer push chunk %lld\n", n_chunks);
            //cerr << "producer" << my_rank << " push chunk " << n_chunks << endl;
        }
    }
    dq.SetCompleted();
    producerDone = 1;
    printf("rank%d producer baba cost %lf\n", my_rank, GetTime() - t0);
}

struct formatpe_data {
    vector <neoReference> *data1;
    vector <neoReference> *data2;
    rabbit::fq::FastqDataPairChunk *fqdatachunk;
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
double tsum4_3 = 0;
double tsum4_4 = 0;
double tsum5 = 0;



void PeQc::ProcessNgsData(bool &proDone, std::vector <neoReference> &data1, std::vector <neoReference> &data2, rabbit::fq::FastqDataPairChunk *fqdatachunk, qc_data *para, rabbit::fq::FastqDataPool *fastqPool) {

    //fprintf(stderr, "consumer %d data size %d %d\n", my_rank, data1.size(), data2.size());
    double tt0 = GetTime();
    vector <neoReference> pass_data1(data1);
    vector <neoReference> pass_data2(data2);
    vector <dupInfo> dups;
    //if(cmd_info_->write_data_) {
    //    pass_data1.resize(data1.size());
    //    pass_data2.resize(data2.size());
    //}
    if(cmd_info_->state_duplicate_) {
        dups.resize(data1.size());
    }
    tsum1 += GetTime() - tt0;

    tt0 = GetTime();
    if(fqdatachunk != NULL) {
        para->data1_ = &data1;
        para->data2_ = &data2;
        para->pass_data1_ = &pass_data1;
        para->pass_data2_ = &pass_data2;
        para->dups = &dups;
        __real_athread_spawn((void *)slave_ngspefunc, para, 1);
        athread_join();
    }
    tsum2 += GetTime() - tt0;
    //fprintf(stderr, "consumer %d data_pass size %d %d\n", my_rank, pass_data1.size(), pass_data2.size());

    tt0 = GetTime(); 
    if(cmd_info_->state_duplicate_) {
        for(auto item : dups) {
            auto key = item.key;
            auto kmer32 = item.kmer32;
            auto gc = item.gc;
            if (duplicate_->counts_[key] == 0) {
                duplicate_->counts_[key] = 1;
                duplicate_->dups_[key] = kmer32;
                duplicate_->gcs_[key] = gc;
            } else {
                if (duplicate_->dups_[key] == kmer32) {
                    duplicate_->counts_[key]++;
                    if (duplicate_->gcs_[key] > gc) duplicate_->gcs_[key] = gc;
                } else if (duplicate_->dups_[key] > kmer32) {
                    duplicate_->dups_[key] = kmer32;
                    duplicate_->counts_[key] = 1;
                    duplicate_->gcs_[key] = gc;
                }
            }
        }
    }
    tsum3 += GetTime() - tt0;

    tt0 = GetTime();
    if (cmd_info_->write_data_) {
        double tt1 = GetTime();

        double tt2 = GetTime();
        int out_len1 = 0;
        for(auto item : pass_data1){
            if(item.lname == 0) continue;
            out_len1 += item.lname + item.lseq + item.lstrand + item.lqual + 4;
        }
        int out_len2 = 0;
        for(auto item : pass_data2){
            if(item.lname == 0) continue;
            out_len2 += item.lname + item.lseq + item.lstrand + item.lqual + 4;
        }
        tsum4_1_1 += GetTime() - tt2;

        if (cmd_info_->interleaved_out_) {
            printf("interleaved out TODO...\n");
        } else {
            tt2 = GetTime();
            char *out_data1;
            if(out_len1) out_data1 = new char[out_len1];
            else out_data1 = (char *)(NULL);
            if(out_len1) memset(out_data1, 0, sizeof(char) * out_len1);
            char *out_data2;
            if(out_len2) out_data2 = new char[out_len2];
            else out_data2 = (char *)(NULL);
            if(out_len2) memset(out_data2, 0, sizeof(char) * out_len2);
            tsum4_1_2 += GetTime() - tt2;

            //fprintf(stderr, "consumer %d out_len %d %d\n", my_rank, out_len1, out_len2);

            tt2 = GetTime();
            //OPT: reduce memcpy times, merge contiguous datas.
            int pos = 0;
            char* start_pos = NULL;
            int tmp_len = 0;
            for(int i = 0; i < pass_data1.size(); i++) {
                if(pass_data1[i].lname == 0 || pass_data1[i].pname + pass_data1[i].lname + pass_data1[i].lseq + pass_data1[i].lstrand + 3 != pass_data1[i].pqual) { 
                    //this record has been trimmed
                    if(tmp_len) memcpy(out_data1 + pos, start_pos, tmp_len);
                    pos += tmp_len;
                    if(pass_data1[i].lname) Read2Chars(pass_data1[i], out_data1, pos);
                    tmp_len = 0;
                    if(i < pass_data1.size() - 1) {
                        start_pos = (char*)pass_data1[i + 1].base + pass_data1[i + 1].pname;
                    }
                    continue;
                }
                if((char*)pass_data1[i].base + pass_data1[i].pname != start_pos + tmp_len) {
                    //record is complete, but can extend to last record
                    if(tmp_len) memcpy(out_data1 + pos, start_pos, tmp_len);
                    pos += tmp_len;
                    tmp_len = 0;
                    start_pos = (char*)pass_data1[i].base + pass_data1[i].pname;
                } 
                tmp_len += pass_data1[i].lname + pass_data1[i].lseq + pass_data1[i].lstrand + pass_data1[i].lqual + 4;
            }
            if(tmp_len) {
                memcpy(out_data1 + pos, start_pos, tmp_len);
                pos += tmp_len;
            }


            //int pos = 0;
            //for (auto item: pass_data1) {
            //    if(item.lname == 0) continue;
            //    Read2Chars(item, out_data1, pos);
            //}
            ASSERT(pos == out_len1);

            //OPT: reduce memcpy times, merge contiguous datas.
            pos = 0;
            start_pos = NULL;
            tmp_len = 0;
            for(int i = 0; i < pass_data2.size(); i++) {
                if(pass_data2[i].lname == 0 || pass_data2[i].pname + pass_data2[i].lname + pass_data2[i].lseq + pass_data2[i].lstrand + 3 != pass_data2[i].pqual) { 
                    //this record has been trimmed
                    if(tmp_len) memcpy(out_data2 + pos, start_pos, tmp_len);
                    pos += tmp_len;
                    if(pass_data2[i].lname) Read2Chars(pass_data2[i], out_data2, pos);
                    tmp_len = 0;
                    if(i < pass_data2.size() - 1) {
                        start_pos = (char*)pass_data2[i + 1].base + pass_data2[i + 1].pname;
                    }
                    continue;
                }
                if((char*)pass_data2[i].base + pass_data2[i].pname != start_pos + tmp_len) {
                    //record is complete, but can extend to last record
                    if(tmp_len) memcpy(out_data2 + pos, start_pos, tmp_len);
                    pos += tmp_len;
                    tmp_len = 0;
                    start_pos = (char*)pass_data2[i].base + pass_data2[i].pname;
                } 
                tmp_len += pass_data2[i].lname + pass_data2[i].lseq + pass_data2[i].lstrand + pass_data2[i].lqual + 4;
            }
            if(tmp_len) {
                memcpy(out_data2 + pos, start_pos, tmp_len);
                pos += tmp_len;
            }


            //pos = 0;
            //for (auto item: pass_data2) {
            //    if(item.lname == 0) continue;
            //    Read2Chars(item, out_data2, pos);
            //}
            ASSERT(pos == out_len2);
            tsum4_1_3 += GetTime() - tt2;

            tsum4_1 += GetTime() - tt1;
        
            tt1 = GetTime();
            int now_size1 = out_len1;
            char *now_pos1 = out_data1;
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

            long long now_pos_base1 = now_pos1_ + pre_sizes1[my_rank];

            for(int ii = 0; ii < comm_size; ii++) {
                now_pos1_ += now_sizes1[ii];
            }



            int now_size2 = out_len2;
            char *now_pos2 = out_data2;
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

            long long now_pos_base2 = now_pos2_ + pre_sizes2[my_rank];

            for(int ii = 0; ii < comm_size; ii++) {
                now_pos2_ += now_sizes2[ii];
            }
            tsum4_2 += GetTime() - tt1;

            tt1 = GetTime();
            if(proDone == 0) {
                int proStop = producerStop.load();
                int stops[comm_size];
                stops[0] = proStop;
                MPI_Barrier(MPI_COMM_WORLD);
                if(my_rank) {
                    MPI_Send(&proStop, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                } else {
                    for(int ii = 1; ii < comm_size; ii++) {
                        MPI_Recv(&(stops[ii]), 1, MPI_INT, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    } 
                }
                MPI_Barrier(MPI_COMM_WORLD);
                if(my_rank == 0) {
                    for(int ii = 1; ii < comm_size; ii++) {
                        MPI_Send(stops, comm_size, MPI_INT, ii, 0, MPI_COMM_WORLD);
                    }
                } else {
                    MPI_Recv(stops, comm_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                int all_stop = 1;
                int some_stop = 0;
                for(int ii = 0; ii < comm_size; ii++) {
                    if(stops[ii] == 0) all_stop = 0;
                    if(stops[ii] == 1) some_stop = 1;
                }
                if(all_stop) {
                    int64_t chunk_sizes[comm_size];
                    chunk_sizes[0] = now_chunks;
                    MPI_Barrier(MPI_COMM_WORLD);
                    if(my_rank) {
                        MPI_Send(&now_chunks, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
                    } else {
                        for(int ii = 1; ii < comm_size; ii++) {
                            MPI_Recv(&(chunk_sizes[ii]), 1, MPI_LONG, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        } 
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                    if(my_rank == 0) {
                        for(int ii = 1; ii < comm_size; ii++) {
                            MPI_Send(chunk_sizes, comm_size, MPI_LONG, ii, 0, MPI_COMM_WORLD);
                        }
                    } else {
                        MPI_Recv(chunk_sizes, comm_size, MPI_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                    for(int ii = 0; ii < comm_size; ii++) {
                        mx_chunks = max(mx_chunks, chunk_sizes[ii]);
                    }
                    consumerCommDone = 1;
                    proDone = 1;
                } else if(some_stop){
                    int goonn = 1;
                    while(goonn) {
                        usleep(1000);
                        int proStop = producerStop.load();
                        int stops[comm_size];
                        stops[0] = proStop;
                        MPI_Barrier(MPI_COMM_WORLD);
                        if(my_rank) {
                            MPI_Send(&proStop, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                        } else {
                            for(int ii = 1; ii < comm_size; ii++) {
                                MPI_Recv(&(stops[ii]), 1, MPI_INT, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            } 
                        }
                        MPI_Barrier(MPI_COMM_WORLD);
                        if(my_rank == 0) {
                            for(int ii = 1; ii < comm_size; ii++) {
                                MPI_Send(stops, comm_size, MPI_INT, ii, 0, MPI_COMM_WORLD);
                            }
                        } else {
                            MPI_Recv(stops, comm_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        }
                        int all_stop = 1;
                        int some_stop = 0;
                        for(int ii = 0; ii < comm_size; ii++) {
                            if(stops[ii] == 0) all_stop = 0;
                            if(stops[ii] == 1) some_stop = 1;
                        }
                        //cerr << "===waiting stop" << my_rank << " " << proStop << endl;
                        //printf("consumer %d all stop: %d\n", my_rank, all_stop);
                        if(all_stop) {
                            int64_t chunk_sizes[comm_size];
                            chunk_sizes[0] = now_chunks;
                            MPI_Barrier(MPI_COMM_WORLD);
                            if(my_rank) {
                                MPI_Send(&now_chunks, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
                            } else {
                                for(int ii = 1; ii < comm_size; ii++) {
                                    MPI_Recv(&(chunk_sizes[ii]), 1, MPI_LONG, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                } 
                            }
                            MPI_Barrier(MPI_COMM_WORLD);
                            if(my_rank == 0) {
                                for(int ii = 1; ii < comm_size; ii++) {
                                    MPI_Send(chunk_sizes, comm_size, MPI_LONG, ii, 0, MPI_COMM_WORLD);
                                }
                            } else {
                                MPI_Recv(chunk_sizes, comm_size, MPI_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            }
                            for(int ii = 0; ii < comm_size; ii++) {
                                mx_chunks = max(mx_chunks, chunk_sizes[ii]);
                            }
                            consumerCommDone = 1;
                            proDone = 1;
                            goonn = 0;
                        }

                    }
                }
            }
            tsum4_3 += GetTime() - tt1;


            tt1 = GetTime();

            //mylock.lock();
            while (queueNumNow >= queueSizeLim) {
#ifdef Verbose
                //printf("waiting to push a chunk to out queue1\n");
#endif
                usleep(1000);
            }
            out_queue_[queueP2++] = make_pair(make_pair(out_data1, make_pair(out_len1, now_pos_base1)), make_pair(out_data2, make_pair(out_len2, now_pos_base2)));
            queueNumNow++;
            //fprintf(stderr, "consumer %d push to queue\n", my_rank);
         

//            //mylock.lock();
//            while (queueNumNow1 >= queueSizeLim1 || queueNumNow2 >= queueSizeLim2) {
//#ifdef Verbose
//                //printf("waiting to push a chunk to out queue1\n");
//#endif
//                usleep(10000);
//            }
//            out_queue1_[queue1P2++] = make_pair(out_data1, make_pair(out_len1, now_pos_base1));
//            queueNumNow1++;
//            out_queue2_[queue2P2++] = make_pair(out_data2, make_pair(out_len2, now_pos_base2));
//            queueNumNow2++;
//            //mylock.unlock();
            tsum4_4 += GetTime() - tt1;
        }
    }
    tsum4 += GetTime() - tt0;

    //fprintf(stderr, "consumer %d release1\n", my_rank);
    tt0 = GetTime();
    if(fqdatachunk != NULL) {
        fastqPool->Release(fqdatachunk->left_part);
        fastqPool->Release(fqdatachunk->right_part);
    }
    //fprintf(stderr, "consumer %d release2\n", my_rank);
    tsum5 += GetTime() - tt0;
}

void PeQc::ConsumerPeFastqTask64(ThreadInfo **thread_infos, rabbit::fq::FastqDataPool *fastqPool) {
    rabbit::int64 id = 0;
    rabbit::fq::FastqDataPairChunk *fqdatachunk;
    qc_data para;
    para.cmd_info_ = cmd_info_;
    para.thread_info_ = thread_infos;
    para.bit_len = 0;

    bool allProDone = 0;
    if(cmd_info_->state_duplicate_) {
        para.bit_len = duplicate_->key_len_base_;
    }
    athread_init();

    double t0 = GetTime();
    double t_format = 0;
    bool p_overWhile = 0;
    vector<rabbit::fq::FastqDataPairChunk *> fqdatachunks;
    while(true) {
        bool pushNull = 0;
        while (p_queueNumNow == 0) {
            int proDone = producerDone.load();
            int proStop = producerStop.load();
            if (proStop && proDone) {
                p_overWhile = 1;
                printf("consumer %d process done\n", my_rank);
                break;
            } else if(proStop && !proDone) {
                pushNull = 1;
                break;
            }
            //fprintf(stderr, "consumer %d waiting\n", my_rank);
            usleep(1000);
        }
        //fprintf(stderr, "consumer %d pushNull: %d\n", my_rank, pushNull);
        if (p_overWhile) {
            consumerCommDone = 1;
            break;
        }
        if(pushNull) {
            fqdatachunks.clear();
            fqdatachunks.push_back(NULL);
        } else {
            fqdatachunks = p_out_queue_[p_queueP1++];
            p_queueNumNow--;
        } 
        
        double tt0 = GetTime();
        vector<neoReference> data1[64];
        vector<neoReference> data2[64];
        formatpe_data para2[64];
        uint64_t counts[64] = {0};
        for(int i = fqdatachunks.size(); i < 64; i++) {
            fqdatachunks.push_back(NULL);
        }
        int fq_size = fqdatachunks.size();
        for(int i = 0; i < 64; i++) {
            para2[i].fqSize = fq_size;
            para2[i].my_rank = my_rank;
        }
        //fprintf(stderr, "consumer %d fq_size %d\n", my_rank, fq_size);
        for(int i = 0; i < fq_size; i++) {
            //if(fqdatachunks[i] != NULL) fprintf(stderr, "rank%d %p %p - %d %d\n", my_rank, fqdatachunks[i]->left_part, fqdatachunks[i]->right_part, fqdatachunks[i]->left_part->size, fqdatachunks[i]->right_part->size); 
            if(fqdatachunks[i] != NULL) data1[i].resize(fqdatachunks[i]->left_part->size / (cmd_info_->seq_len_ * 2));
            if(fqdatachunks[i] != NULL) data2[i].resize(fqdatachunks[i]->right_part->size / (cmd_info_->seq_len_ * 2));
            para2[i].data1 = &(data1[i]);
            para2[i].data2 = &(data2[i]);
            para2[i].fqdatachunk = fqdatachunks[i];
            para2[i].res = counts;
            para2[i].fqSize = fq_size;
        }

        __real_athread_spawn((void *)slave_formatpefunc, para2, 1);
        athread_join();
        t_format += GetTime() - tt0;

        for(int i = 0; i < fq_size; i++) {
            //fprintf(stderr, "consumer rank%d count %d\n", my_rank, counts[i]);
            data1[i].resize(counts[i]);
            data2[i].resize(counts[i]);
            ProcessNgsData(allProDone, data1[i], data2[i], fqdatachunks[i], &para, fastqPool);
        }
        fqdatachunks.clear();
    }
    printf("NGSnew tot cost %lf\n", GetTime() - t0);
    printf("NGSnew format cost %lf\n", t_format);
    printf("NGSnew resize cost %lf\n", tsum1);
    printf("NGSnew slave cost %lf\n", tsum2);
    printf("NGSnew dup cost %lf\n", tsum3);
    printf("NGSnew write cost %lf (%lf [%lf %lf %lf], %lf, %lf, %lf)\n", tsum4, tsum4_1, tsum4_1_1, tsum4_1_2, tsum4_1_3, tsum4_2, tsum4_3, tsum4_4);
    printf("NGSnew release cost %lf\n", tsum5);
    done_thread_number_++;
    athread_halt();
    fprintf(stderr, "consumer %d in func done\n", my_rank);
}


/**
 * @brief get fastq data chunks from the data queue and do QC for them
 * @param thread_info : thread information
 * @param fastq_data_pool :a fastq data pool, it will be used to release data chunk
 * @param dq : data queue
 */
void PeQc::ConsumerPeFastqTask(ThreadInfo **thread_infos, rabbit::fq::FastqDataPool *fastqPool,
        rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> &dq) {
    rabbit::int64 id = 0;
    rabbit::fq::FastqDataPairChunk *fqdatachunk;
    qc_data para;
    para.cmd_info_ = cmd_info_;
    para.thread_info_ = thread_infos;
    para.bit_len = 0;
    bool proDone = 0;
    if(cmd_info_->state_duplicate_) {
        para.bit_len = duplicate_->key_len_base_;
    }
    athread_init();
    double t0 = GetTime();
    double tsum1 = 0;
    double tsum2 = 0;
    double tsum3 = 0;
    double tsum4 = 0;
    double tsum4_1 = 0;
    double tsum4_2 = 0;
    while (dq.Pop(id, fqdatachunk)) {
        double tt0 = GetTime();
        tsum1 += GetTime() - tt0;
        tt0 = GetTime();
        vector <neoReference> data1;
        vector <neoReference> data2;
        if(fqdatachunk != NULL) rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk *) (fqdatachunk->left_part), data1, true);
        if(fqdatachunk != NULL) rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk *) (fqdatachunk->right_part), data2, true);
        if(data1.size() != data2.size()) printf("GG size pe\n");
        vector <neoReference> pass_data1(data1);
        vector <neoReference> pass_data2(data2);
        vector <dupInfo> dups;
        if(cmd_info_->state_duplicate_) {
            dups.resize(data1.size());
        }
        tsum2 += GetTime() - tt0;
        tt0 = GetTime();
        if(fqdatachunk != NULL) {
            para.data1_ = &data1;
            para.data2_ = &data2;
            para.pass_data1_ = &pass_data1;
            para.pass_data2_ = &pass_data2;
            para.dups = &dups;
            __real_athread_spawn((void *)slave_ngspefunc, &para, 1);
            athread_join();
        }
        //printf("pass size1 %d\n", pass_data1.size());
        //printf("pass size2 %d\n", pass_data2.size());
        tsum3 += GetTime() - tt0;
        tt0 = GetTime(); 
        if(cmd_info_->state_duplicate_) {
            for(auto item : dups) {
                auto key = item.key;
                auto kmer32 = item.kmer32;
                auto gc = item.gc;
                if (duplicate_->counts_[key] == 0) {
                    duplicate_->counts_[key] = 1;
                    duplicate_->dups_[key] = kmer32;
                    duplicate_->gcs_[key] = gc;
                } else {
                    if (duplicate_->dups_[key] == kmer32) {
                        duplicate_->counts_[key]++;
                        if (duplicate_->gcs_[key] > gc) duplicate_->gcs_[key] = gc;
                    } else if (duplicate_->dups_[key] > kmer32) {
                        duplicate_->dups_[key] = kmer32;
                        duplicate_->counts_[key] = 1;
                        duplicate_->gcs_[key] = gc;
                    }
                }
            }
        }
        if (cmd_info_->write_data_) {
            int out_len1 = 0;
            for(auto item : pass_data1){
                if(item.lname == 0) continue;
                out_len1 += item.lname + item.lseq + item.lstrand + item.lqual + 4;
            }
            int out_len2 = 0;
            for(auto item : pass_data2){
                if(item.lname == 0) continue;
                out_len2 += item.lname + item.lseq + item.lstrand + item.lqual + 4;
            }
            if (cmd_info_->interleaved_out_) {
                printf("interleaved out TODO...\n");
            } else {
                char *out_data1;
                if(out_len1) out_data1 = new char[out_len1];
                else out_data1 = (char *)(NULL);
                int pos = 0;
                for (auto item: pass_data1) {
                    if(item.lname == 0) continue;
                    Read2Chars(item, out_data1, pos);
                }
                ASSERT(pos == out_len1);

                char *out_data2;
                if(out_len2) out_data2 = new char[out_len2];
                else out_data2 = (char *)(NULL);
                pos = 0;
                for (auto item: pass_data2) {
                    if(item.lname == 0) continue;
                    Read2Chars(item, out_data2, pos);
                }
                ASSERT(pos == out_len2);

                int now_size1 = out_len1;
                char *now_pos1 = out_data1;
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

                long long now_pos_base1 = now_pos1_ + pre_sizes1[my_rank];

                for(int ii = 0; ii < comm_size; ii++) {
                    now_pos1_ += now_sizes1[ii];
                }



                int now_size2 = out_len2;
                char *now_pos2 = out_data2;
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

                long long now_pos_base2 = now_pos2_ + pre_sizes2[my_rank];

                for(int ii = 0; ii < comm_size; ii++) {
                    now_pos2_ += now_sizes2[ii];
                }


                if(proDone == 0) {
                    int proStop = producerStop.load();
                    int stops[comm_size];
                    stops[0] = proStop;
                    MPI_Barrier(MPI_COMM_WORLD);
                    if(my_rank) {
                        MPI_Send(&proStop, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                    } else {
                        for(int ii = 1; ii < comm_size; ii++) {
                            MPI_Recv(&(stops[ii]), 1, MPI_INT, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        } 
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                    if(my_rank == 0) {
                        for(int ii = 1; ii < comm_size; ii++) {
                            MPI_Send(stops, comm_size, MPI_INT, ii, 0, MPI_COMM_WORLD);
                        }
                    } else {
                        MPI_Recv(stops, comm_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                    int all_stop = 1;
                    int some_stop = 0;
                    for(int ii = 0; ii < comm_size; ii++) {
                        if(stops[ii] == 0) all_stop = 0;
                        if(stops[ii] == 1) some_stop = 1;
                    }
                    if(all_stop) {
                        int64_t chunk_sizes[comm_size];
                        chunk_sizes[0] = now_chunks;
                        MPI_Barrier(MPI_COMM_WORLD);
                        if(my_rank) {
                            MPI_Send(&now_chunks, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
                        } else {
                            for(int ii = 1; ii < comm_size; ii++) {
                                MPI_Recv(&(chunk_sizes[ii]), 1, MPI_LONG, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            } 
                        }
                        MPI_Barrier(MPI_COMM_WORLD);
                        if(my_rank == 0) {
                            for(int ii = 1; ii < comm_size; ii++) {
                                MPI_Send(chunk_sizes, comm_size, MPI_LONG, ii, 0, MPI_COMM_WORLD);
                            }
                        } else {
                            MPI_Recv(chunk_sizes, comm_size, MPI_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        }
                        for(int ii = 0; ii < comm_size; ii++) {
                            mx_chunks = max(mx_chunks, chunk_sizes[ii]);
                        }
                        consumerCommDone = 1;
                        proDone = 1;
                    } else if(some_stop){
                        int goonn = 1;
                        while(goonn) {
                            usleep(1000);
                            int proStop = producerStop.load();
                            int stops[comm_size];
                            stops[0] = proStop;
                            MPI_Barrier(MPI_COMM_WORLD);
                            if(my_rank) {
                                MPI_Send(&proStop, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                            } else {
                                for(int ii = 1; ii < comm_size; ii++) {
                                    MPI_Recv(&(stops[ii]), 1, MPI_INT, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                } 
                            }
                            MPI_Barrier(MPI_COMM_WORLD);
                            if(my_rank == 0) {
                                for(int ii = 1; ii < comm_size; ii++) {
                                    MPI_Send(stops, comm_size, MPI_INT, ii, 0, MPI_COMM_WORLD);
                                }
                            } else {
                                MPI_Recv(stops, comm_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            }
                            int all_stop = 1;
                            int some_stop = 0;
                            for(int ii = 0; ii < comm_size; ii++) {
                                if(stops[ii] == 0) all_stop = 0;
                                if(stops[ii] == 1) some_stop = 1;
                            }
                            //cerr << "===waiting stop" << my_rank << " " << proStop << endl;
                            //printf("waiting stop%d %d\n", my_rank, all_stop);
                            if(all_stop) {
                                int64_t chunk_sizes[comm_size];
                                chunk_sizes[0] = now_chunks;
                                MPI_Barrier(MPI_COMM_WORLD);
                                if(my_rank) {
                                    MPI_Send(&now_chunks, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
                                } else {
                                    for(int ii = 1; ii < comm_size; ii++) {
                                        MPI_Recv(&(chunk_sizes[ii]), 1, MPI_LONG, ii, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                    } 
                                }
                                MPI_Barrier(MPI_COMM_WORLD);
                                if(my_rank == 0) {
                                    for(int ii = 1; ii < comm_size; ii++) {
                                        MPI_Send(chunk_sizes, comm_size, MPI_LONG, ii, 0, MPI_COMM_WORLD);
                                    }
                                } else {
                                    MPI_Recv(chunk_sizes, comm_size, MPI_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                }
                                for(int ii = 0; ii < comm_size; ii++) {
                                    mx_chunks = max(mx_chunks, chunk_sizes[ii]);
                                }
                                consumerCommDone = 1;
                                proDone = 1;
                                goonn = 0;
                            }

                        }
                    }
                }


                //mylock.lock();
                while (queueNumNow >= queueSizeLim) {
#ifdef Verbose
                    //printf("waiting to push a chunk to out queue1\n");
#endif
                    usleep(1000);
                }
                out_queue_[queueP2++] = make_pair(make_pair(out_data1, make_pair(out_len1, now_pos_base1)), make_pair(out_data2, make_pair(out_len2, now_pos_base2)));
                queueNumNow++;
         
//                //mylock.lock();
//                while (queueNumNow1 >= queueSizeLim1 || queueNumNow2 >= queueSizeLim2) {
//#ifdef Verbose
//                    //printf("waiting to push a chunk to out queue1\n");
//#endif
//                    usleep(1000);
//                }
//                out_queue1_[queue1P2++] = make_pair(out_data1, make_pair(out_len1, now_pos_base1));
//                queueNumNow1++;
//                out_queue2_[queue2P2++] = make_pair(out_data2, make_pair(out_len2, now_pos_base2));
//                queueNumNow2++;
                //printf("push chunk %lld === %d %lld, %d %lld\n", id, out_len1, now_pos_base1, out_len2, now_pos_base2);
                //mylock.unlock();
            }
        }
        if(fqdatachunk != NULL) {
            fastqPool->Release(fqdatachunk->left_part);
            fastqPool->Release(fqdatachunk->right_part);
        }
    }
    printf("NGSnew tot cost %lf\n", GetTime() - t0);
    printf("NGSnew producer cost %lf\n", tsum1);
    printf("NGSnew format cost %lf\n", tsum2);
    printf("NGSnew slave cost %lf\n", tsum3);
    printf("NGSnew write cost %lf\n", tsum4);
    printf("NGSnew write1 cost %lf\n", tsum4_1);
    printf("NGSnew write2 cost %lf\n", tsum4_2);
    done_thread_number_++;
    athread_halt();
    //fprintf(stderr, "rank%d 11111111111111\n", my_rank);
}


void PeQc::ConsumerPeInterFastqTask(ThreadInfo *thread_info, rabbit::fq::FastqDataPool *fastqPool,
        rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {
    printf("interleaved mode TODO...\n");
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
    pair<CIPair, CIPair> now0;
    CIPair now1;
    CIPair now2;
    while (true) {
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
        now0 = out_queue_[queueP1++];
        now1 = now0.first;
        now2 = now0.second;
        queueNumNow--;
        if (out_is_zip_) {
            if (cmd_info_->use_pigz_) {
                printf("use pigz TODO...\n");
            } else {
                int written1 = gzwrite(zip_out_stream1, now1.first, now1.second.first);
                if (written1 != now1.second.first) {
                    printf("gzwrite error\n");
                    exit(0);
                }
                delete[] now1.first;

                int written2 = gzwrite(zip_out_stream2, now2.first, now2.second.first);
                if (written2 != now2.second.first) {
                    printf("gzwrite error\n");
                    exit(0);
                }
                delete[] now2.first;
            }
        } else {

            tot_size1 += now1.second.first;
            if(now1.second.first) {
#ifdef use_out_mem
                //memcpy(OutMemData1 + now1.second.second, now1.first, now1.second.first); 
                memcpy(OutMemData1, now1.first, now1.second.first); 
#else

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

#endif
                delete[] now1.first;
                //fprintf(stderr, "writer %d write del done\n", my_rank);
            }


            tot_size2 += now2.second.first;
            if(now2.second.first) {

#ifdef use_out_mem
                //memcpy(OutMemData2 + now2.second.second, now2.first, now2.second.first); 
                memcpy(OutMemData2, now2.first, now2.second.first); 
#else

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

#endif
                delete[] now2.first;
                //fprintf(stderr, "writer %d write del done\n", my_rank);
            }
            //fprintf(stderr, "writer %d %d done\n", my_rank, cnt++);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    //fprintf(stderr, "writer rank%d write donedone\n", my_rank);
    if (out_is_zip_) {
        if (cmd_info_->use_pigz_) {

        } else {
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
        }

    } else {
#ifdef use_mpi_file
        //fprintf(stderr, "writer%d close file 1\n", my_rank);
        MPI_File_close(&fh1);
        MPI_File_close(&fh2);
        //fprintf(stderr, "writer%d close file 2\n", my_rank);
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
#ifdef Verbose
    printf("writer%d cost %lf, tot size %lld %lld\n", my_rank, GetTime() - t0, tot_size1, tot_size2);
    fprintf(stderr, "writer %d in func done\n", my_rank);
    //printf("write %d cost %.5f --- %lld\n", my_rank, GetTime() - t0, tot_size);
    //
#endif
}



/**
 * @brief a function to write pe data from out_data1 queue to file1
 */
void PeQc::WriteSeFastqTask1() {
#ifdef Verbose
    double t0 = GetTime();
#endif
    int cnt = 0;
    long long tot_size = 0;
    bool overWhile = 0;
    CIPair now;
    while (true) {
        while (queueNumNow1 == 0) {
            if (done_thread_number_ == cmd_info_->thread_number_) {
                overWhile = 1;
                printf("rank %d write done\n", my_rank);
                break;
            }

            //fprintf(stderr, "rank %d write wait\n", my_rank);
            usleep(1000);
        }
        if (overWhile) break;
        now = out_queue1_[queue1P1++];
        queueNumNow1--;
        //fprintf(stderr, "rank %d write get %p\n", my_rank, now.first);
        if (out_is_zip_) {
            if (cmd_info_->use_pigz_) {
                printf("use pigz TODO...\n");
            } else {
                int written = gzwrite(zip_out_stream1, now.first, now.second.first);
                if (written != now.second.first) {
                    printf("gzwrite error\n");
                    exit(0);
                }
                delete[] now.first;
            }
        } else {
            tot_size += now.second.first;
            if(now.second.first) {

#ifdef use_mpi_file
                //fprintf(stderr, "rank %d write seek %lld\n", my_rank, now.second.first);
                MPI_File_seek(fh1, now.second.second, MPI_SEEK_SET);
                //fprintf(stderr, "rank %d write ww %lld\n", my_rank, now.second.first);
                MPI_File_write(fh1, now.first, now.second.first, MPI_CHAR, &status1);
                //fprintf(stderr, "rank %d write ww done\n", my_rank);
#else
                fseek(out_stream1_, now.second.second, SEEK_SET);
                //fprintf(stderr, "rank %d write ww %lld\n", my_rank, now.second.first);
                fwrite(now.first, sizeof(char), now.second.first, out_stream1_);
                //fprintf(stderr, "rank %d write ww done\n", my_rank);
#endif
                delete[] now.first;
                //fprintf(stderr, "rank %d write del done\n", my_rank);
            }
            //fprintf(stderr, "write%d %d done\n", my_rank, cnt++);
        }
    }
    //fprintf(stderr, "111111 rank%d donedone1\n", my_rank);
    MPI_Barrier(MPI_COMM_WORLD);
    //fprintf(stderr, "111111 rank%d donedone2\n", my_rank);
    if (out_is_zip_) {
        if (cmd_info_->use_pigz_) {

        } else {
            if (zip_out_stream1) {
                gzflush(zip_out_stream1, Z_FINISH);
                gzclose(zip_out_stream1);
                zip_out_stream1 = NULL;
            }
        }

    } else {
#ifdef use_mpi_file
        //fprintf(stderr, "111111 rank%d donedone3\n", my_rank);
        MPI_File_close(&fh1);
        //fprintf(stderr, "111111 rank%d donedone4\n", my_rank);
#else
        fclose(out_stream1_);
        if(my_rank == 0) {
            truncate(cmd_info_->out_file_name1_.c_str(), sizeof(char) * now_pos1_);
        }
#endif
    }
#ifdef Verbose
    cerr << "write" << my_rank << " cost " << GetTime() - t0 << ", tot size " << tot_size << endl;
    //printf("write %d cost %.5f --- %lld\n", my_rank, GetTime() - t0, tot_size);
    //
#endif
}

/**
 * @brief a function to write pe data from out_data2 queue to file2
 */
void PeQc::WriteSeFastqTask2() {
#ifdef Verbose
    double t0 = GetTime();
#endif
    int cnt = 0;
    long long tot_size = 0;
    bool overWhile = 0;
    CIPair now;
    while (true) {
        while (queueNumNow2 == 0) {
            if (done_thread_number_ == cmd_info_->thread_number_) {
                overWhile = 1;
                printf("rank %d write done\n", my_rank);
                break;
            }

            //fprintf(stderr, "rank %d write wait\n", my_rank);
            usleep(1000);
        }
        if (overWhile) break;
        now = out_queue2_[queue2P1++];
        queueNumNow2--;
        //fprintf(stderr, "rank %d write get %p\n", my_rank, now.first);
        if (out_is_zip_) {
            if (cmd_info_->use_pigz_) {
                printf("use pigz TODO...\n");
            } else {
                int written = gzwrite(zip_out_stream2, now.first, now.second.first);
                if (written != now.second.first) {
                    printf("gzwrite error\n");
                    exit(0);
                }
                delete[] now.first;
            }
        } else {
            tot_size += now.second.first;
            if(now.second.first) {

#ifdef use_mpi_file
                //fprintf(stderr, "rank %d write seek %lld\n", my_rank, now.second.first);
                MPI_File_seek(fh2, now.second.second, MPI_SEEK_SET);
                //fprintf(stderr, "rank %d write ww %lld\n", my_rank, now.second.first);
                MPI_File_write(fh2, now.first, now.second.first, MPI_CHAR, &status2);
                //fprintf(stderr, "rank %d write ww done\n", my_rank);
#else
                fseek(out_stream2_, now.second.second, SEEK_SET);
                //fprintf(stderr, "rank %d write ww %lld\n", my_rank, now.second.first);
                fwrite(now.first, sizeof(char), now.second.first, out_stream2_);
                //fprintf(stderr, "rank %d write ww done\n", my_rank);
#endif
                delete[] now.first;
                //fprintf(stderr, "rank %d write del done\n", my_rank);
            }
            //fprintf(stderr, "write%d %d done\n", my_rank, cnt++);
        }
    }
    //fprintf(stderr, "222222 rank%d donedone1\n", my_rank);
    MPI_Barrier(MPI_COMM_WORLD);
    //fprintf(stderr, "222222 rank%d donedone2\n", my_rank);
    if (out_is_zip_) {
        if (cmd_info_->use_pigz_) {

        } else {
            if (zip_out_stream2) {
                gzflush(zip_out_stream2, Z_FINISH);
                gzclose(zip_out_stream2);
                zip_out_stream2 = NULL;
            }
        }

    } else {
#ifdef use_mpi_file
        //fprintf(stderr, "222222 rank%d donedone3\n", my_rank);
        MPI_File_close(&fh2);
        //fprintf(stderr, "222222 rank%d donedone4\n", my_rank);
#else
        fclose(out_stream2_);
        if(my_rank == 0) {
            truncate(cmd_info_->out_file_name2_.c_str(), sizeof(char) * now_pos2_);
        }
#endif
    }
#ifdef Verbose
    cerr << "write" << my_rank << " cost " << GetTime() - t0 << ", tot size " << tot_size << endl;
    //printf("write %d cost %.5f --- %lld\n", my_rank, GetTime() - t0, tot_size);
    //
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


void PeQc::ProcessPeFastq() {
    if(my_rank == 0) block_size = 6 * (1 << 20);
    else block_size = 4 * (1 << 20);
    auto *fastqPool = new rabbit::fq::FastqDataPool(Q_lim * 64 * 2, block_size);
    //rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> queue1(Q_lim, 1);
    auto **p_thread_info = new ThreadInfo *[slave_num];
    for (int t = 0; t < slave_num; t++) {
        p_thread_info[t] = new ThreadInfo(cmd_info_, true);
    }
    rabbit::uint32 tmpSize = 1 << 20;
    if (cmd_info_->seq_len_ <= 200) tmpSize = 1 << 14;
    
    fqFileReader = new rabbit::fq::FastqFileReader(cmd_info_->in_file_name1_, *fastqPool, cmd_info_->in_file_name2_, in_is_zip_, tmpSize);
#ifdef use_in_mem
    fqFileReader->MemDataReader();
    fqFileReader->MemDataReader2();
    double ttt = GetTime();
#endif


    //tagg

    //thread *write_thread12;
    //thread *write_thread1;
    //thread *write_thread2;
    //if (cmd_info_->write_data_) {
    //    write_thread12 = new thread(bind(&PeQc::WriteSeFastqTask12, this));
    //    //write_thread1 = new thread(bind(&PeQc::WriteSeFastqTask1, this));
    //    //if (cmd_info_->interleaved_out_ == 0)
    //    //    write_thread2 = new thread(bind(&PeQc::WriteSeFastqTask2, this));
    //}

    //thread producer(bind(&PeQc::ProducerPeFastqTask, this, cmd_info_->in_file_name1_, cmd_info_->in_file_name2_, fastqPool, ref(queue1)));
    thread producer(bind(&PeQc::ProducerPeFastqTask64, this, cmd_info_->in_file_name1_, cmd_info_->in_file_name2_, fastqPool));
    //thread consumer(bind(&PeQc::ConsumerPeFastqTask, this, p_thread_info, fastqPool, ref(queue1)));
    thread consumer(bind(&PeQc::ConsumerPeFastqTask64, this, p_thread_info, fastqPool));

    if (cmd_info_->write_data_) {
        WriteSeFastqTask12();
        printf("write%d done\n", my_rank);
    }

    //bool writeThreadJoined = false;
    //bool consumerThreadJoined = false;
    //bool producerThreadJoined = false;
    //if(cmd_info_->write_data_ == 0) writeThreadJoined = 1;

    //while (!writeThreadJoined || !consumerThreadJoined || !producerThreadJoined) {
    //    if (cmd_info_->write_data_ && !writeThreadJoined && write_thread12->joinable()) {
    //        write_thread12->join();
    //        printf("write%d done\n", my_rank);
    //        writerDone1 = 1;
    //        writeThreadJoined = true;
    //    }

    //    if (!consumerThreadJoined && consumer.joinable()) {
    //        consumer.join();
    //        printf("consumer%d done\n", my_rank);
    //        consumerThreadJoined = true;
    //    }

    //    if (!producerThreadJoined && producer.joinable()) {
    //        producer.join();
    //        printf("producer%d done\n", my_rank);
    //        producerThreadJoined = true;
    //    }
    //    printf("waiting %d %d %d\n", producerThreadJoined, consumerThreadJoined, writeThreadJoined);
    //    usleep(1000);
    //}

    //if (cmd_info_->write_data_) {
    //    //write_thread1->join();
    //    write_thread12->join();
    //    printf("write%d done\n", my_rank);
    //    writerDone1 = 1;
    //    //if (cmd_info_->interleaved_out_ == 0) {
    //    //    write_thread2->join();
    //    //    printf("write%d done2\n", my_rank);
    //    //    writerDone2 = 1;
    //    //}
    //}

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
    auto pre_state_tmp1 = State::MergeStates(pre_vec_state1);
    auto pre_state_tmp2 = State::MergeStates(pre_vec_state2);
    auto aft_state_tmp1 = State::MergeStates(aft_vec_state1);
    auto aft_state_tmp2 = State::MergeStates(aft_vec_state2);
    printf("merge1 done\n");

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
        bool res = checkStates(pre_state_tmp1, tmp_state);
        if(!res) cerr << "GGGGGGGG" << endl;
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
        bool res = checkStates(pre_state_tmp2, tmp_state);
        if(!res) cerr << "GGGGGGGG" << endl;
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
    printf("aft merge done\n");

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
        return;
    }



    auto pre_state1 = State::MergeStates(pre_state_mpis1);
    auto pre_state2 = State::MergeStates(pre_state_mpis2);
    auto aft_state1 = State::MergeStates(aft_state_mpis1);
    auto aft_state2 = State::MergeStates(aft_state_mpis2);


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
#ifdef use_in_mem
    fprintf(stderr, "TOT TIME %lf\n", GetTime() - ttt);
    fqFileReader->ReleaseMemData();
#endif
    delete fqFileReader;
#ifdef use_out_mem
    delete[] OutMemData1;
    delete[] OutMemData2;
#endif

}
