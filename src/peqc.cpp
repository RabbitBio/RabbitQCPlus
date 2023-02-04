// O
// Created by ylf9811 on 2021/7/6.
//
//

#include "peqc.h"
using namespace std;

struct dupInfo{
    uint32_t key;
    uint64_t kmer32;
    uint8_t gc;
};

int block_size;

struct qc_data {
    ThreadInfo **thread_info_;
    CmdInfo *cmd_info_;
    vector <dupInfo> *dups;
    vector <neoReference> *data1_;
    vector <neoReference> *data2_;
    vector <neoReference> *pass_data1_;
    vector <neoReference> *pass_data2_;
    int *cnt;
    int bit_len;
};

extern "C" {
#include <athread.h>
#include <pthread.h>
    void slave_ngspefunc();
}


//extern void slave_ngspefunc(qc_data *para);


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
    out_queue1_ = NULL;
    out_queue2_ = NULL;
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
    int64_t right_pos = GetNextFastq(tmp_chunk, 0, res_size);
    if(my_rank == 0) right_pos = 0;
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
    
    //MPI_Barrier(MPI_COMM_WORLD);
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


    if (cmd_info1->write_data_) {
        out_queue1_ = new CIPair[1 << 20];
        queue1P1 = 0;
        queue1P2 = 0;
        queueNumNow1 = 0;
        queueSizeLim1 = 64;
        if (cmd_info_->interleaved_out_ == 0) {
            out_queue2_ = new CIPair[1 << 20];
            queue2P1 = 0;
            queue2P2 = 0;
            queueNumNow2 = 0;
            queueSizeLim2 = 64;
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
    if (cmd_info_->write_data_) {
        delete out_queue1_;
        if (cmd_info_->interleaved_out_ == 0)
            delete out_queue2_;
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


void PeQc::ProducerPeFastqTask(string file, string file2, rabbit::fq::FastqDataPool *fastqPool,
        rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> &dq) {
#ifdef Verbose
    double t0 = GetTime();
#endif
    //cerr << "part_sizes ";
    //for(int i = 0; i < comm_size; i++) cerr << part_sizes[i] << " ";
    //cerr << endl;
    
    rabbit::fq::FastqFileReader *fqFileReader;
    rabbit::uint32 tmpSize = 1 << 20;
    if (cmd_info_->seq_len_ <= 200) tmpSize = 1 << 14;
    fqFileReader = new rabbit::fq::FastqFileReader(file, *fastqPool, file2, in_is_zip_, tmpSize);
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
     
            fqdatachunk = fqFileReader->readNextPairChunk(offset, part_sizes[my_rank]);
            if (fqdatachunk == NULL) {

                printf("null\n");
                if(cmd_info_->write_data_) {
                    //printf("producer%d stop\n", my_rank);
                    //cerr << "producer" << my_rank << "stop" << endl;
                    now_chunks = n_chunks;
                    producerStop = 1;
                    while(consumerCommDone == 0) {
                        usleep(10000);
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
                printf("producer done === %lld\n", n_chunks);
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
    delete fqFileReader;
    producerDone = 1;
    printf("producer cost %lf\n", GetTime() - t0);
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
                //char *out_data = new char[out_len1 + out_len2];
                //int pos = 0;
                //int len = min(pass_data1.size(), pass_data2.size());
                //for (int i = 0; i < len; i++) {
                //    auto item1 = pass_data1[i];
                //    auto item2 = pass_data2[i];
                //    if(item1.lname == 0) continue;
                //    if(item2.lname == 0) continue;
                //    Read2Chars(item1, out_data, pos);
                //    Read2Chars(item2, out_data, pos);
                //}
                //mylock.lock();
                //while (queueNumNow1 >= queueSizeLim1) {
#ifdef Verbose
                //    //printf("waiting to push a chunk to out queue1\n");
#endif
                //    usleep(100);
                //}
                //out_queue1_[queue1P2++] = {out_data, pos};
                //queueNumNow1++;
                //mylock.unlock();
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
                //MPI_Barrier(MPI_COMM_WORLD);
                if(my_rank == 0) {
                    for(int ii = 1; ii < comm_size; ii++) {
                        MPI_Send(now_sizes1, comm_size, MPI_INT, ii, 0, MPI_COMM_WORLD);
                    }
                } else {
                    MPI_Recv(now_sizes1, comm_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                //MPI_Barrier(MPI_COMM_WORLD);
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
                //MPI_Barrier(MPI_COMM_WORLD);
                if(my_rank == 0) {
                    for(int ii = 1; ii < comm_size; ii++) {
                        MPI_Send(now_sizes2, comm_size, MPI_INT, ii, 0, MPI_COMM_WORLD);
                    }
                } else {
                    MPI_Recv(now_sizes2, comm_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                //MPI_Barrier(MPI_COMM_WORLD);
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
                    //printf("stop%d %d %d\n", my_rank, all_stop, some_stop);
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
                            usleep(100000);
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



                mylock.lock();
                while (queueNumNow1 >= queueSizeLim1 || queueNumNow2 >= queueSizeLim2) {
#ifdef Verbose
                    //printf("waiting to push a chunk to out queue1\n");
#endif
                    usleep(10000);
                }
                out_queue1_[queue1P2++] = make_pair(out_data1, make_pair(out_len1, now_pos_base1));
                queueNumNow1++;
                out_queue2_[queue2P2++] = make_pair(out_data2, make_pair(out_len2, now_pos_base2));
                queueNumNow2++;
                //printf("push chunk %lld === %d %lld, %d %lld\n", id, out_len1, now_pos_base1, out_len2, now_pos_base2);
                mylock.unlock();
            }
        }
        if(fqdatachunk != NULL) {
            fastqPool->Release(fqdatachunk->left_part);
            fastqPool->Release(fqdatachunk->right_part);
        }
    }
    athread_halt();
    printf("NGSnew tot cost %lf\n", GetTime() - t0);
    printf("NGSnew producer cost %lf\n", tsum1);
    printf("NGSnew format cost %lf\n", tsum2);
    printf("NGSnew slave cost %lf\n", tsum3);
    printf("NGSnew write cost %lf\n", tsum4);
    printf("NGSnew write1 cost %lf\n", tsum4_1);
    printf("NGSnew write2 cost %lf\n", tsum4_2);
    done_thread_number_++;
}


void PeQc::ConsumerPeInterFastqTask(ThreadInfo *thread_info, rabbit::fq::FastqDataPool *fastqPool,
        rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {
    printf("interleaved mode TODO...\n");
    /*
       rabbit::int64 id = 0;
       rabbit::fq::FastqDataChunk *fqdatachunk;
       while (dq.Pop(id, fqdatachunk)) {
       vector<neoReference> data;
       vector<neoReference> pass_data1, pass_data2;
       rabbit::fq::chunkFormat(fqdatachunk, data, true);
       ASSERT(data.size() % 2 == 0);
       int out_len1 = 0, out_len2 = 0;
       for (int i = 0; i + 2 <= data.size(); i += 2) {
       auto item1 = data[i];
       auto item2 = data[i + 1];
       thread_info->pre_state1_->StateInfo(item1);
       thread_info->pre_state2_->StateInfo(item2);
       if (cmd_info_->state_duplicate_) {
       duplicate_->statPair(item1, item2);
       }
       if (cmd_info_->add_umi_) {
       umier_->ProcessPe(item1, item2);
       }

    //do pe sequence trim
    bool trim_res1 = filter_->TrimSeq(item1, cmd_info_->trim_front1_, cmd_info_->trim_tail1_);
    bool trim_res2 = filter_->TrimSeq(item2, cmd_info_->trim_front2_, cmd_info_->trim_tail2_);

    if (trim_res1 && trim_res2 && cmd_info_->trim_polyg_) {
    PolyX::trimPolyG(item1, item2, cmd_info_->trim_poly_len_);
    }

    if (trim_res1 && trim_res2 && cmd_info_->trim_polyx_) {
    PolyX::trimPolyX(item1, item2, cmd_info_->trim_poly_len_);
    }

    //do pe overlap analyze
    OverlapRes overlap_res;
    if (trim_res1 && trim_res2 && cmd_info_->analyze_overlap_) {
    overlap_res = Adapter::AnalyzeOverlap(item1, item2, cmd_info_->overlap_diff_limit_,
    cmd_info_->overlap_require_);
    int now_size = cmd_info_->max_insert_size_;
    if (overlap_res.overlapped) {
    if (overlap_res.offset > 0)
    now_size = item1.lseq + item2.lseq - overlap_res.overlap_len;
    else
    now_size = overlap_res.overlap_len;
    }
    now_size = min(now_size, cmd_info_->max_insert_size_);
    thread_info->insert_size_dist_[now_size]++;
    }
    if (trim_res1 && trim_res2 && cmd_info_->correct_data_) {
    Adapter::CorrectData(item1, item2, overlap_res, cmd_info_->isPhred64_);
    }
    if (trim_res1 && trim_res2 && cmd_info_->trim_adapter_) {
    int trimmed;
    if (cmd_info_->print_what_trimmed_) {
    trimmed = Adapter::TrimAdapter(item1, item2, overlap_res.offset, overlap_res.overlap_len,
    thread_info->aft_state1_->adapter_map_,
    thread_info->aft_state2_->adapter_map_, cmd_info_->adapter_len_lim_);
    } else {
    trimmed = Adapter::TrimAdapter(item1, item2, overlap_res.offset, overlap_res.overlap_len);
    }
    if (trimmed) {
    thread_info->aft_state1_->AddTrimAdapter();
    thread_info->aft_state1_->AddTrimAdapter();
    thread_info->aft_state1_->AddTrimAdapterBase(trimmed);
    }
    if (!trimmed) {
    int res1, res2;
    if (cmd_info_->detect_adapter1_) {
    int res1;
    if (cmd_info_->print_what_trimmed_) {
    res1 = Adapter::TrimAdapter(item1, cmd_info_->adapter_seq1_,
    thread_info->aft_state1_->adapter_map_,
        cmd_info_->adapter_len_lim_, false);
} else {
    res1 = Adapter::TrimAdapter(item1, cmd_info_->adapter_seq1_, false);
}
}
if (cmd_info_->detect_adapter2_) {
    int res2;
    if (cmd_info_->print_what_trimmed_) {
        res2 = Adapter::TrimAdapter(item2, cmd_info_->adapter_seq2_,
                thread_info->aft_state2_->adapter_map_,
                cmd_info_->adapter_len_lim_, true);
    } else {
        res2 = Adapter::TrimAdapter(item2, cmd_info_->adapter_seq2_, true);
    }
}
if (res1) {
    thread_info->aft_state1_->AddTrimAdapter();
    thread_info->aft_state1_->AddTrimAdapterBase(res1);
}
if (res2) {
    thread_info->aft_state1_->AddTrimAdapter();
    thread_info->aft_state1_->AddTrimAdapterBase(res2);
}
}
}


//do filer in refs
int filter_res1 = filter_->ReadFiltering(item1, trim_res1, cmd_info_->isPhred64_);
int filter_res2 = filter_->ReadFiltering(item2, trim_res2, cmd_info_->isPhred64_);


int filter_res = max(filter_res1, filter_res2);


if (filter_res == 0) {
    thread_info->aft_state1_->StateInfo(item1);
    thread_info->aft_state2_->StateInfo(item2);
    thread_info->aft_state1_->AddPassReads();
    thread_info->aft_state1_->AddPassReads();
    if (cmd_info_->write_data_) {
        pass_data1.push_back(item1);
        pass_data2.push_back(item2);
        out_len1 += item1.lname + item1.lseq + item1.lstrand + item1.lqual + 4;
        out_len2 += item2.lname + item2.lseq + item2.lstrand + item2.lqual + 4;
    }
} else if (filter_res == 2) {
    thread_info->aft_state1_->AddFailShort();
    thread_info->aft_state1_->AddFailShort();
} else if (filter_res == 3) {
    thread_info->aft_state1_->AddFailLong();
    thread_info->aft_state1_->AddFailLong();
} else if (filter_res == 4) {
    thread_info->aft_state1_->AddFailLowq();
    thread_info->aft_state1_->AddFailLowq();
} else if (filter_res == 1) {
    thread_info->aft_state1_->AddFailN();
    thread_info->aft_state1_->AddFailN();
}
}
if (cmd_info_->write_data_) {
    if (cmd_info_->interleaved_out_) {
        char *out_data = new char[out_len1 + out_len2];
        int pos = 0;
        int len = min(pass_data1.size(), pass_data2.size());
        for (int i = 0; i < len; i++) {
            auto item1 = pass_data1[i];
            auto item2 = pass_data2[i];
            Read2Chars(item1, out_data, pos);
            Read2Chars(item2, out_data, pos);
        }
        mylock.lock();
        while (queueNumNow1 >= queueSizeLim1) {
#ifdef Verbose
            //printf("waiting to push a chunk to out queue1\n");
#endif
            usleep(100);
        }
        out_queue1_[queue1P2++] = {out_data, pos};
        queueNumNow1++;
        mylock.unlock();
    } else {
        if (pass_data1.size() > 0 && pass_data2.size() > 0) {
            char *out_data1 = new char[out_len1];
            int pos = 0;
            for (auto item: pass_data1) {
                Read2Chars(item, out_data1, pos);
            }
            ASSERT(pos == out_len1);
            char *out_data2 = new char[out_len2];
            pos = 0;
            for (auto item: pass_data2) {
                Read2Chars(item, out_data2, pos);
            }
            ASSERT(pos == out_len2);

            mylock.lock();
            while (queueNumNow1 >= queueSizeLim1 || queueNumNow2 >= queueSizeLim2) {
#ifdef Verbose
                //printf("waiting to push a chunk to out queue1\n");
#endif
                usleep(100);
            }
            out_queue1_[queue1P2++] = {out_data1, out_len1};
            queueNumNow1++;
            out_queue2_[queue2P2++] = {out_data2, out_len2};
            queueNumNow2++;
            mylock.unlock();
        }
    }
}

fastqPool->Release(fqdatachunk);
}
done_thread_number_++;
*/
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
    //CIPair now2;
    while (true) {
        while (queueNumNow1 == 0) {
            if (done_thread_number_ == cmd_info_->thread_number_) {
                overWhile = 1;
                break;
            }
            usleep(10000);
        }
        if (overWhile) break;
        now = out_queue1_[queue1P1++];
        //now2 = out_queue2_[queue2P1++];
        queueNumNow1--;
        //queueNumNow2--;
        if (out_is_zip_) {
            if (cmd_info_->use_pigz_) {
                printf("use pigz TODO...\n");
                //                while (pigzQueueNumNow1 > pigzQueueSizeLim1) {
                //
                //#ifdef Verbose
                //                    //printf("waiting to push a chunk to pigz queue1\n");
                //#endif
                //                    usleep(100);
                //                }
                //                pigzQueue1->enqueue(now);
                //                pigzQueueNumNow1++;
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
                fseek(out_stream1_, now.second.second, SEEK_SET);
                fwrite(now.first, sizeof(char), now.second.first, out_stream1_);
                delete[] now.first;
            }

            //tot_size += now2.second.first;
            //if(now2.second.first) {
            //    fseek(out_stream2_, now2.second.second, SEEK_SET);
            //    fwrite(now2.first, sizeof(char), now2.second.first, out_stream2_);
            //    delete[] now2.first;
            //}
            //printf("write%d %d %d %lld done0\n", my_rank, cnt++, now.second.first, now.second.second);

        }
    }
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
        fclose(out_stream1_);
        //fclose(out_stream2_);
        if(my_rank == 0) {
            truncate(cmd_info_->out_file_name1_.c_str(), sizeof(char) * now_pos1_);
            //truncate(cmd_info_->out_file_name2_.c_str(), sizeof(char) * now_pos2_);
        }
    }
#ifdef Verbose
    printf("write1 cost %.5f\n", GetTime() - t0);

    //cerr << "write" << my_rank << " cost " << GetTime() - t0 << ", tot size " << tot_size << endl;
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
                break;
            }
            usleep(10000);
        }
        if (overWhile) break;
        now = out_queue2_[queue2P1++];
        queueNumNow2--;
        if (out_is_zip_) {
            if (cmd_info_->use_pigz_) {
                printf("use pigz TODO...\n");
                //                while (pigzQueueNumNow2 > pigzQueueSizeLim2) {
                //
                //#ifdef Verbose
                //                    //printf("waiting to push a chunk to pigz queue2\n");
                //#endif
                //                    usleep(100);
                //                }
                //                pigzQueue2->enqueue(now);
                //                pigzQueueNumNow2++;
            } else {
                int written = gzwrite(zip_out_stream2, now.first, now.second.first);
                if (written != now.second.first) {
                    printf("GG");
                    exit(0);
                }

                delete[] now.first;
            }
        } else {
            tot_size += now.second.first;
            if(now.second.first) {
                fseek(out_stream2_, now.second.second, SEEK_SET);
                fwrite(now.first, sizeof(char), now.second.first, out_stream2_);
                delete[] now.first;
            }
            //printf("write%d %d %d %lld done1\n", my_rank, cnt++, now.second.first, now.second.second);
        }
    }

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
        fclose(out_stream2_);
        if(my_rank == 0) {
            truncate(cmd_info_->out_file_name2_.c_str(), sizeof(char) * now_pos2_);
        }
    }
#ifdef Verbose
    printf("write1 cost %.5f\n", GetTime() - t0);
    //cerr << "write" << my_rank << " cost " << GetTime() - t0 << ", tot size " << tot_size << endl;
#endif
}

void PeQc::PugzTask1() {
#ifdef Verbose
    printf("pugz1 start\n");
    double t0 = GetTime();
#endif
    main_pugz(cmd_info_->in_file_name1_, cmd_info_->pugz_threads_, pugzQueue1, &producerDone);
    pugzDone1 = 1;
#ifdef Verbose
    printf("pugz1 done, cost %.6f\n", GetTime() - t0);
#endif
}

void PeQc::PugzTask2() {
#ifdef Verbose
    printf("pugz2 start\n");
    double t0 = GetTime();
#endif
    main_pugz(cmd_info_->in_file_name2_, cmd_info_->pugz_threads_, pugzQueue2, &producerDone);
    pugzDone2 = 1;
#ifdef Verbose
    printf("pugz2 done, cost %.6f\n", GetTime() - t0);
#endif
}


void PeQc::PigzTask1() {
    int cnt = 10;

    char **infos = new char *[10];
    infos[0] = "./pigz";
    infos[1] = "-p";
    int th_num = cmd_info_->pigz_threads_;
    //    printf("th num is %d\n", th_num);
    string th_num_s = to_string(th_num);
    //    printf("th num s is %s\n", th_num_s.c_str());
    //    printf("th num s len is %d\n", th_num_s.length());

    infos[2] = new char[th_num_s.length() + 1];
    memcpy(infos[2], th_num_s.c_str(), th_num_s.length());
    infos[2][th_num_s.length()] = '\0';
    infos[3] = "-k";

    string tmp_level = to_string(cmd_info_->compression_level_);
    tmp_level = "-" + tmp_level;
    infos[4] = new char[tmp_level.length() + 1];
    memcpy(infos[4], tmp_level.c_str(), tmp_level.length());
    infos[4][tmp_level.length()] = '\0';


    infos[5] = "-f";
    infos[6] = "-b";
    infos[7] = "4096";
    string out_name1 = cmd_info_->out_file_name1_;
    string out_file = out_name1.substr(0, out_name1.find(".gz"));
    //    printf("th out_file is %s\n", out_file.c_str());
    //    printf("th out_file len is %d\n", out_file.length());
    infos[8] = new char[out_file.length() + 1];
    memcpy(infos[8], out_file.c_str(), out_file.length());
    infos[8][out_file.length()] = '\0';
    infos[9] = "-v";
    main_pigz(cnt, infos, pigzQueue1, &writerDone1, pigzLast1, &pigzQueueNumNow1);
#ifdef Verbose
    printf("pigz1 done\n");
#endif
}

void PeQc::PigzTask2() {

    int cnt = 10;

    char **infos = new char *[10];
    infos[0] = "./pigz";
    infos[1] = "-p";
    int th_num = cmd_info_->pigz_threads_;
    //    printf("th num is %d\n", th_num);
    string th_num_s = to_string(th_num);
    //    printf("th num s is %s\n", th_num_s.c_str());
    //    printf("th num s len is %d\n", th_num_s.length());

    infos[2] = new char[th_num_s.length() + 1];
    memcpy(infos[2], th_num_s.c_str(), th_num_s.length());
    infos[2][th_num_s.length()] = '\0';
    infos[3] = "-k";

    string tmp_level = to_string(cmd_info_->compression_level_);
    tmp_level = "-" + tmp_level;
    infos[4] = new char[tmp_level.length() + 1];
    memcpy(infos[4], tmp_level.c_str(), tmp_level.length());
    infos[4][tmp_level.length()] = '\0';

    infos[5] = "-f";
    infos[6] = "-b";
    infos[7] = "4096";
    string out_name2 = cmd_info_->out_file_name2_;
    string out_file = out_name2.substr(0, out_name2.find(".gz"));
    //    printf("th out_file is %s\n", out_file.c_str());
    //    printf("th out_file len is %d\n", out_file.length());
    infos[8] = new char[out_file.length() + 1];
    memcpy(infos[8], out_file.c_str(), out_file.length());
    infos[8][out_file.length()] = '\0';
    infos[9] = "-v";
    main_pigz(cnt, infos, pigzQueue2, &writerDone2, pigzLast2, &pigzQueueNumNow2);
#ifdef Verbose
    printf("pigz2 done\n");
#endif
}

void PeQc::NGSTask(std::string file, std::string file2, rabbit::fq::FastqDataPool *fastq_data_pool, ThreadInfo **thread_infos){
    rabbit::fq::FastqFileReader *fqFileReader;
    rabbit::uint32 tmpSize = 1 << 20;
    if (cmd_info_->seq_len_ <= 200) tmpSize = 1 << 14;
    fqFileReader = new rabbit::fq::FastqFileReader(file, *fastq_data_pool, file2, in_is_zip_, tmpSize);
    int64_t n_chunks = 0;
    qc_data para;
    para.cmd_info_ = cmd_info_;
    para.thread_info_ = thread_infos;
    para.bit_len = 0;
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
    char enter_char[1];
    enter_char[0] = '\n';
    int nums = 0;
    while (true) {
        double tt0 = GetTime();
        rabbit::fq::FastqDataPairChunk *fqdatachunk;
        fqdatachunk = fqFileReader->readNextPairChunk();
        if (fqdatachunk == NULL) break;
        n_chunks++;
        tsum1 += GetTime() - tt0;
        tt0 = GetTime();
        vector <neoReference> data1;
        vector <neoReference> data2;
        rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk *) (fqdatachunk->left_part), data1, true);
        rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk *) (fqdatachunk->right_part), data2, true);
        if(data1.size() != data2.size()) printf("GG size pe\n");
        //printf("num %d\n", data1.size());
        nums += data1.size();
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
        tsum2 += GetTime() - tt0;
        tt0 = GetTime();

        para.data1_ = &data1;
        para.data2_ = &data2;
        para.pass_data1_ = &pass_data1;
        para.pass_data2_ = &pass_data2;
        para.dups = &dups;
        //athread_spawn_tupled(slave_ngspefunc, &para);
        __real_athread_spawn((void *)slave_ngspefunc, &para, 1);
        athread_join();
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
        if(cmd_info_->write_data_) {

            int tot_len1 = 0;
            double tt00 = GetTime();
            for(auto item : pass_data1){
                if(item.lname == 0) continue;
                tot_len1 += item.lname + item.lseq + item.lstrand + item.lqual + 4;
            }
            char *tmp_out1 = new char[tot_len1];
            char *now_out1 = tmp_out1;
            for(auto item : pass_data1) {
                if(item.lname == 0) continue;
                memcpy(now_out1, (char *)item.base + item.pname, item.lname);
                now_out1 += item.lname;
                memcpy(now_out1, enter_char, 1);
                now_out1++;
                memcpy(now_out1, (char *)item.base + item.pseq, item.lseq);
                now_out1 += item.lseq;
                memcpy(now_out1, enter_char, 1);
                now_out1++;
                memcpy(now_out1, (char *)item.base + item.pstrand, item.lstrand);
                now_out1 += item.lstrand;
                memcpy(now_out1, enter_char, 1);
                now_out1++;
                memcpy(now_out1, (char *)item.base + item.pqual, item.lqual);
                now_out1 += item.lqual;
                memcpy(now_out1, enter_char, 1);
                now_out1++;
            }
            tsum4_1 += GetTime() - tt00;
            tt00 = GetTime();
            //out_stream1_.write(tmp_out1, tot_len1);
            fwrite(tmp_out1, sizeof(char), tot_len1, out_stream1_);
            delete[] tmp_out1;



            int tot_len2 = 0;
            tt00 = GetTime();
            for(auto item : pass_data2){
                if(item.lname == 0) continue;
                tot_len2 += item.lname + item.lseq + item.lstrand + item.lqual + 4;
            }
            char *tmp_out2 = new char[tot_len2];
            char *now_out2 = tmp_out2;
            for(auto item : pass_data2) {
                if(item.lname == 0) continue;
                memcpy(now_out2, (char *)item.base + item.pname, item.lname);
                now_out2 += item.lname;
                memcpy(now_out2, enter_char, 1);
                now_out2++;
                memcpy(now_out2, (char *)item.base + item.pseq, item.lseq);
                now_out2 += item.lseq;
                memcpy(now_out2, enter_char, 1);
                now_out2++;
                memcpy(now_out2, (char *)item.base + item.pstrand, item.lstrand);
                now_out2 += item.lstrand;
                memcpy(now_out2, enter_char, 1);
                now_out2++;
                memcpy(now_out2, (char *)item.base + item.pqual, item.lqual);
                now_out2 += item.lqual;
                memcpy(now_out2, enter_char, 1);
                now_out2++;
            }
            tsum4_1 += GetTime() - tt00;
            tt00 = GetTime();
            //out_stream2_.write(tmp_out2, tot_len2);
            fwrite(tmp_out2, sizeof(char), tot_len2, out_stream2_);
            delete[] tmp_out2;

            tsum4_2 += GetTime() - tt00;
        }
        tsum4 += GetTime() - tt0;
        fastq_data_pool->Release(fqdatachunk->left_part);
        fastq_data_pool->Release(fqdatachunk->right_part);
    }
    athread_halt();
    printf("NGS tot cost %lf\n", GetTime() - t0);
    printf("NGS producer cost %lf\n", tsum1);
    printf("NGS format cost %lf\n", tsum2);
    printf("NGS slave cost %lf\n", tsum3);
    printf("NGS write cost %lf\n", tsum4);
    printf("NGS write1 cost %lf\n", tsum4_1);
    printf("NGS write2 cost %lf\n", tsum4_2);
    printf("nums %d\n", nums);
    delete fqFileReader;
}


/**
 * @brief do QC for pair-end data
 */
void PeQc::ProcessPeFastqOneThread() {
    auto *fastqPool = new rabbit::fq::FastqDataPool(2, 1 << 26);
    auto **p_thread_info = new ThreadInfo *[slave_num];
    for (int t = 0; t < slave_num; t++) {
        p_thread_info[t] = new ThreadInfo(cmd_info_, true);
    }
    NGSTask(cmd_info_->in_file_name1_, cmd_info_->in_file_name2_, fastqPool, p_thread_info);
#ifdef Verbose
    printf("all thrad done\n");
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
    auto pre_state1 = State::MergeStates(pre_vec_state1);
    auto pre_state2 = State::MergeStates(pre_vec_state2);
    auto aft_state1 = State::MergeStates(aft_vec_state1);
    auto aft_state2 = State::MergeStates(aft_vec_state2);
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
}
void PeQc::ProcessPeFastq() {
    if(my_rank == 0) block_size = 6 * (1 << 20);
    else block_size = 4 * (1 << 20);
    auto *fastqPool = new rabbit::fq::FastqDataPool(256, block_size);
    //auto *fastqPool = new rabbit::fq::FastqDataPool(256, 1 << 22);
    rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> queue1(256, 1);
    auto **p_thread_info = new ThreadInfo *[slave_num];
    for (int t = 0; t < slave_num; t++) {
        p_thread_info[t] = new ThreadInfo(cmd_info_, true);
    }
    thread *write_thread1;
    thread *write_thread2;
    if (cmd_info_->write_data_) {
        write_thread1 = new thread(bind(&PeQc::WriteSeFastqTask1, this));
        if (cmd_info_->interleaved_out_ == 0)
            write_thread2 = new thread(bind(&PeQc::WriteSeFastqTask2, this));
    }
    //tagg
    thread producer(bind(&PeQc::ProducerPeFastqTask, this, cmd_info_->in_file_name1_, cmd_info_->in_file_name2_, fastqPool, ref(queue1)));
    thread consumer(bind(&PeQc::ConsumerPeFastqTask, this, p_thread_info, fastqPool, ref(queue1)));
    printf("===%d\n", my_rank);
    producer.join();
    printf("producer%d done\n", my_rank);
    consumer.join();
    printf("consumer%d done\n", my_rank);


    if (cmd_info_->write_data_) {
        write_thread1->join();
        writerDone1 = 1;
        if (cmd_info_->interleaved_out_ == 0) {
            write_thread2->join();
            writerDone2 = 1;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printf("==%d all pro write done\n", my_rank);


#ifdef Verbose
    printf("all thrad done\n");
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
        //State* tmp_state = new State(pre_state_is.c_str(), pre_state_is.length(), cmd_info_, cmd_info_->seq_len_, cmd_info_->qul_range_, false);
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
        //State* tmp_state = new State(pre_state_is.c_str(), pre_state_is.length(), cmd_info_, cmd_info_->seq_len_, cmd_info_->qul_range_, false);
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
}
