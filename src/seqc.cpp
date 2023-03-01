//
// Created by ylf9811 on 2021/7/6.
//
#include "seqc.h"
using namespace std;

/**
 * @brief Construct function
 * @param cmd_info1 : cmd information
 */
SeQc::SeQc(CmdInfo *cmd_info1) {
    cmd_info_ = cmd_info1;
    filter_ = new Filter(cmd_info1);
    done_thread_number_ = 0;
    int out_block_nums = int(1.0 * cmd_info1->in_file_size1_ / cmd_info1->out_block_size_);
    out_queue_ = NULL;

    nowChunkId = 1;
    in_is_zip_ = cmd_info1->in_file_name1_.find(".gz") != string::npos;
    out_is_zip_ = cmd_info1->out_file_name1_.find(".gz") != string::npos;
    if (cmd_info1->write_data_) {
        out_queue_ = new CIPair[1 << 20];
        queueP1 = 0;
        queueP2 = 0;
        queueNumNow = 0;
        queueSizeLim = 1 << 5;
        if (out_is_zip_) {
            if (cmd_info1->use_pigz_) {
                pigzQueueNumNow = 0;
                pigzQueueSizeLim = 1 << 5;
                string out_name1 = cmd_info1->out_file_name1_;
                out_name1 = out_name1.substr(0, out_name1.find(".gz"));
                out_stream_.open(out_name1);
                out_stream_.close();
#ifdef Verbose
                printf("now use pigz to compress output data\n");
#endif
            } else {
#ifdef Verbose
                printf("open gzip stream %s\n", cmd_info1->out_file_name1_.c_str());
#endif
                zip_out_stream = gzopen(cmd_info1->out_file_name1_.c_str(), "w");
                gzsetparams(zip_out_stream, cmd_info1->compression_level_, Z_DEFAULT_STRATEGY);
                gzbuffer(zip_out_stream, 1024 * 1024);
            }
        } else {
#ifdef Verbose
            printf("open stream %s\n", cmd_info1->out_file_name1_.c_str());
#endif
            out_stream_.open(cmd_info1->out_file_name1_);
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
        pugzQueue = new moodycamel::ReaderWriterQueue<pair<char *, int>>(1 << 10);
    }
    if (cmd_info1->use_pigz_) {
        pigzQueue = new moodycamel::ReaderWriterQueue<pair<char *, int>>;
        pigzLast.first = new char[1 << 24];
        pigzLast.second = 0;
    }
    pugzDone = 0;
    producerDone = 0;
    writerDone = 0;
}

SeQc::~SeQc() {
    delete filter_;
    if (cmd_info_->write_data_) {
        delete out_queue_;
    }

    if (cmd_info_->state_duplicate_) {
        delete duplicate_;
    }
    if (cmd_info_->add_umi_) {
        delete umier_;
    }
}

void SeQc::careProcess() {
    
    printf("pairmode %s\n", cmd_info_->pairmode_.c_str());
    printf("coverage %d\n", cmd_info_->coverage_);
    vector<char*>paras;
    paras.push_back("./RabbitQCPlus");
    paras.push_back("-i");
    paras.push_back((char*)(cmd_info_->in_file_name1_.data()));
    paras.push_back("-d");
    paras.push_back("./");
    paras.push_back("-o");
    paras.push_back("tmp.fq");
    paras.push_back("-c");
    string str_coverage = to_string(cmd_info_->coverage_);
    printf("str_coverage %s\n", str_coverage.c_str());
    paras.push_back((char*)(str_coverage.data()));
    paras.push_back("-t");
    string str_thread = to_string(cmd_info_->thread_number_);
    printf("str_thread %s\n", str_thread.c_str());
    paras.push_back((char*)(str_thread.data()));
    paras.push_back("--pairmode");
    if(cmd_info_->pairmode_ == "SE") paras.push_back("SE");
    else paras.push_back("PE");
    for(int i = 0; i < paras.size(); i++) {
        printf("%s ", paras[i]);
    }
    printf("\n");

    printf("start care part...\n");

    main_correction(paras.size(), &(paras[0]));

    printf("care end\n");

}

/**
 * @brief get fastq data chunk from fastq_data_pool and put it into data queue
 * @param file : fastq file name, which is also the input file
 * @param fastq_data_pool : fastq data pool
 * @param dq : a data queue
 */

void SeQc::ProducerSeFastqTask(string file, rabbit::fq::FastqDataPool *fastq_data_pool,
                               rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {

    rabbit::fq::FastqFileReader *fqFileReader;
    rabbit::uint32 tmpSize = 1 << 20;
    if (cmd_info_->seq_len_ <= 200) tmpSize = 1 << 14;
    fqFileReader = new rabbit::fq::FastqFileReader(file, *fastq_data_pool, "", in_is_zip_, tmpSize);
    int64_t n_chunks = 0;

    if (cmd_info_->use_pugz_) {
        pair<char *, int> last_info;
        last_info.first = new char[1 << 20];
        last_info.second = 0;
        while (true) {
            rabbit::fq::FastqDataChunk *fqdatachunk;
            fqdatachunk = fqFileReader->readNextChunk(pugzQueue, &pugzDone, last_info);
            if (fqdatachunk == NULL) break;
            n_chunks++;
            dq.Push(n_chunks, fqdatachunk);
        }
        delete[] last_info.first;
    } else {
        while (true) {
            rabbit::fq::FastqDataChunk *fqdatachunk;
            fqdatachunk = fqFileReader->readNextChunk();
            if (fqdatachunk == NULL) break;
            n_chunks++;
            dq.Push(n_chunks, fqdatachunk);
        }
    }


    dq.SetCompleted();
    delete fqFileReader;
    producerDone = 1;
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

/**
 * @brief get fastq data chunks from the data queue and do QC for them
 * @param thread_info : thread information
 * @param fastq_data_pool :a fastq data pool, it will be used to release data chunk
 * @param dq : data queue
 */
void SeQc::ConsumerSeFastqTask(ThreadInfo *thread_info, rabbit::fq::FastqDataPool *fastq_data_pool,
                               rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {
    rabbit::int64 id = 0;
    rabbit::fq::FastqDataChunk *fqdatachunk;

    if (cmd_info_->is_TGS_) {
        while (dq.Pop(id, fqdatachunk)) {
            vector<neoReference> data;
            rabbit::fq::chunkFormat(fqdatachunk, data, true);
            for (auto item: data) {
                thread_info->TGS_state_->tgsStatRead(item, cmd_info_->isPhred64_);
            }
            fastq_data_pool->Release(fqdatachunk);
        }
    } else {
        while (dq.Pop(id, fqdatachunk)) {
            vector<neoReference> data;
            vector<neoReference> pass_data;
            rabbit::fq::chunkFormat(fqdatachunk, data, true);
            int out_len = 0;
            for (auto item: data) {
                thread_info->pre_state1_->StateInfo(item);
                if (cmd_info_->state_duplicate_) {
                    duplicate_->statRead(item);
                }

                if (cmd_info_->add_umi_) {
                    umier_->ProcessSe(item);
                }
                bool trim_res = filter_->TrimSeq(item, cmd_info_->trim_front1_, cmd_info_->trim_tail1_);

                if (trim_res && cmd_info_->trim_polyg_) {
                    PolyX::trimPolyG(item, cmd_info_->trim_poly_len_);
                }

                if (trim_res && cmd_info_->trim_polyx_) {
                    PolyX::trimPolyX(item, cmd_info_->trim_poly_len_);
                }

                if (trim_res && cmd_info_->trim_adapter_) {
                    int res = 0;
                    bool is_trimmed = false;
                    if(cmd_info_->detect_adapter1_) {
                        if (cmd_info_->print_what_trimmed_) {
                            res = Adapter::TrimAdapter(item, cmd_info_->adapter_seq1_,
                                    thread_info->aft_state1_->adapter_map_, cmd_info_->adapter_len_lim_, false);
                        } else {
                            res = Adapter::TrimAdapter(item, cmd_info_->adapter_seq1_, false);
                        }
                        if (res) {
                            is_trimmed = true;
                            thread_info->aft_state1_->AddTrimAdapter();
                            thread_info->aft_state1_->AddTrimAdapterBase(res);
                        }
                    }
                    if(cmd_info_->adapter_from_fasta_.size() > 0) {
                        if (cmd_info_->print_what_trimmed_) {
                            res = Adapter::TrimAdapters(item, cmd_info_->adapter_from_fasta_,
                                    thread_info->aft_state1_->adapter_map_, cmd_info_->adapter_len_lim_, false);
                        } else {
                            res = Adapter::TrimAdapters(item, cmd_info_->adapter_from_fasta_, false);
                        }
                        if (res) {
                            if(!is_trimmed) thread_info->aft_state1_->AddTrimAdapter();
                            thread_info->aft_state1_->AddTrimAdapterBase(res);
                        }

                    }
                }
                
                int filter_res = filter_->ReadFiltering(item, trim_res, cmd_info_->isPhred64_);
                if (filter_res == 0) {
                    thread_info->aft_state1_->StateInfo(item);
                    thread_info->aft_state1_->AddPassReads();
                    if (cmd_info_->write_data_) {
                        pass_data.push_back(item);
                        out_len += item.lname + item.lseq + item.lstrand + item.lqual + 4;
                    }
                } else if (filter_res == 1) {
                    thread_info->aft_state1_->AddFailN();
                } else if (filter_res == 2) {
                    thread_info->aft_state1_->AddFailShort();
                } else if (filter_res == 3) {
                    thread_info->aft_state1_->AddFailLong();
                } else if (filter_res == 4) {
                    thread_info->aft_state1_->AddFailLowq();
                }
            }

            if (cmd_info_->write_data_) {
                if (pass_data.size() > 0) {
                    char *out_data = new char[out_len];
                    int pos = 0;
                    for (auto item: pass_data) {
                        Read2Chars(item, out_data, pos);
                    }
                    ASSERT(pos == out_len);
                    if(cmd_info_->notKeepOrder_ == 0) {
                        while(nowChunkId != id) {
                            usleep(100);
                        }
                    }

                    mylock.lock();
                    while (queueNumNow >= queueSizeLim) {
#ifdef Verbose
                        //printf("waiting to push a chunk to out queue %d\n",out_len);
#endif
                        usleep(100);
                    }
                    
                    out_queue_[queueP2++] = {out_data, pos};
                    queueNumNow++;
                    nowChunkId++;
                    mylock.unlock();
                }
            }

            fastq_data_pool->Release(fqdatachunk);
        }
    }
    done_thread_number_++;
}

/**
 * @brief a function to write data from out_data queue to file
 */
void SeQc::WriteSeFastqTask() {
#ifdef Verbose
    double t0 = GetTime();
#endif
    int cnt = 0;
    bool overWhile = 0;
    pair<char *, int> now;
    while (true) {
        while (queueNumNow == 0) {
            if (done_thread_number_ == cmd_info_->thread_number_) {
                overWhile = 1;
                break;
            }
            usleep(100);
        }
        if (overWhile) break;
        now = out_queue_[queueP1++];
        queueNumNow--;
        if (out_is_zip_) {
            if (cmd_info_->use_pigz_) {
                while (pigzQueueNumNow > pigzQueueSizeLim) {
#ifdef Verbose
                    //printf("waiting to push a chunk to pigz queue\n");
#endif
                    usleep(100);
                }
                pigzQueue->enqueue(now);
                pigzQueueNumNow++;
            } else {
                int written = gzwrite(zip_out_stream, now.first, now.second);
                if (written != now.second) {
                    printf("gzwrite error\n");
                    exit(0);
                }
                delete[] now.first;
            }
        } else {
            out_stream_.write(now.first, now.second);
            delete[] now.first;
        }
    }

    if (out_is_zip_) {
        if (cmd_info_->use_pigz_) {
            //            out_stream_.close();

        } else {
            if (zip_out_stream) {
                gzflush(zip_out_stream, Z_FINISH);
                gzclose(zip_out_stream);
                zip_out_stream = NULL;
            }
        }
    } else {
        out_stream_.close();
    }
#ifdef Verbose
    printf("write cost %.5f\n", GetTime() - t0);
#endif
}


/**
 * @brief do pugz
 */

void SeQc::PugzTask() {
#ifdef Verbose
    printf("pugz start\n");
    auto t0 = GetTime();
#endif
    main_pugz(cmd_info_->in_file_name1_, cmd_info_->pugz_threads_, pugzQueue, &producerDone);
#ifdef Verbose
    printf("pugz cost %.5f\n", GetTime() - t0);
#endif
    pugzDone = 1;
}


void SeQc::PigzTask() {
    /*
       argc 9
       argv ./pigz
       argv -p
       argv 16
       argv -k
       argv -4
       argv -f
       argv -b
       argv -4096
       argv p.fq
       */
    int cnt = 9;

    char **infos = new char *[9];
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
    main_pigz(cnt, infos, pigzQueue, &writerDone, pigzLast, &pigzQueueNumNow);
#ifdef Verbose
    printf("pigz done\n");
#endif
}


/**
 * @brief do QC for single-end data
 */

void SeQc::ProcessSeFastq() {
#ifdef Verbose
    double t0 = GetTime();
#endif
    thread *carer;
    if(cmd_info_->do_correction_with_care_) {
        carer = new thread(bind(&SeQc::careProcess, this));
        carer->join();
        delete carer;
        cmd_info_->in_file_name1_ = "./tmp.fq";
    }


    thread *pugzer;

    if (cmd_info_->use_pugz_) {
        pugzer = new thread(bind(&::SeQc::PugzTask, this));
    }


    auto *fastqPool = new rabbit::fq::FastqDataPool(32, 1 << 22);
    rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> queue1(32, 1);

    auto **p_thread_info = new ThreadInfo *[cmd_info_->thread_number_];
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        p_thread_info[t] = new ThreadInfo(cmd_info_, false);
    }
    thread *write_thread;
    if (cmd_info_->write_data_) {
        write_thread = new thread(bind(&SeQc::WriteSeFastqTask, this));
    }
    thread *pigzer;
    if (cmd_info_->use_pigz_) {
        pigzer = new thread(bind(&SeQc::PigzTask, this));
    }

    thread producer(
            bind(&SeQc::ProducerSeFastqTask, this, cmd_info_->in_file_name1_, fastqPool, ref(queue1)));
    auto **threads = new thread *[cmd_info_->thread_number_];
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        threads[t] = new thread(
                bind(&SeQc::ConsumerSeFastqTask, this, p_thread_info[t], fastqPool, ref(queue1)));
    }
    if (cmd_info_->use_pugz_) {
        pugzer->join();
    }

    producer.join();
#ifdef Verbose
    printf("producer cost %.4f\n", GetTime() - t0);
#endif
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        threads[t]->join();
    }
#ifdef Verbose
    printf("consumer cost %.4f\n", GetTime() - t0);
#endif
    if (cmd_info_->write_data_) {
        write_thread->join();
        writerDone = 1;
    }
    if (cmd_info_->use_pigz_) {
        pigzer->join();
    }
#ifdef Verbose
    printf("all thrad done\n");
    printf("now merge thread info\n");
#endif
    vector<State *> pre_vec_state;
    vector<State *> aft_vec_state;

    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        pre_vec_state.push_back(p_thread_info[t]->pre_state1_);
        aft_vec_state.push_back(p_thread_info[t]->aft_state1_);
    }
    auto pre_state = State::MergeStates(pre_vec_state);
    auto aft_state = State::MergeStates(aft_vec_state);
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
    string srr_name = cmd_info_->in_file_name1_;
    srr_name = PaseFileName(srr_name);
    Repoter::ReportHtmlSe(srr_name + "_RabbitQCPlus.html", pre_state, aft_state, cmd_info_->in_file_name1_,
                          dupRate * 100.0);


    delete pre_state;
    delete aft_state;

    delete fastqPool;


    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        delete threads[t];
        delete p_thread_info[t];
    }

    delete[] threads;
    delete[] p_thread_info;
    if (cmd_info_->write_data_) {
        delete write_thread;
    }
}

void SeQc::ProcessSeTGS() {

    thread *pugzer;
    if (cmd_info_->use_pugz_) {
        pugzer = new thread(bind(&::SeQc::PugzTask, this));
    }


    auto *fastqPool = new rabbit::fq::FastqDataPool(32, 1 << 22);
    rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> queue1(32, 1);

    auto **p_thread_info = new ThreadInfo *[cmd_info_->thread_number_];
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        p_thread_info[t] = new ThreadInfo(cmd_info_, false);
    }
    thread producer(
            bind(&SeQc::ProducerSeFastqTask, this, cmd_info_->in_file_name1_, fastqPool, ref(queue1)));
    auto **threads = new thread *[cmd_info_->thread_number_];
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        threads[t] = new thread(
                bind(&SeQc::ConsumerSeFastqTask, this, p_thread_info[t], fastqPool, ref(queue1)));
    }

    if (cmd_info_->use_pugz_) {
        pugzer->join();
    }

    producer.join();
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        threads[t]->join();
    }
#ifdef Verbose
    printf("all thrad done\n");
    printf("now merge thread info\n");
#endif
    vector<TGSStats *> vec_state;

    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        vec_state.push_back(p_thread_info[t]->TGS_state_);
    }
    auto mer_state = TGSStats::merge(vec_state);
#ifdef Verbose
    printf("merge done\n");
#endif
    printf("\nprint TGS state info :\n");

    //    report3(mer_state);
    mer_state->CalReadsLens();

    mer_state->print();
    string srr_name = cmd_info_->in_file_name1_;
    srr_name = PaseFileName(srr_name);
    string command = cmd_info_->command_;
    Repoter::ReportHtmlTGS(srr_name + "_RabbitQCPlus.html", command, mer_state, cmd_info_->in_file_name1_);

    delete fastqPool;
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        delete threads[t];
        delete p_thread_info[t];
    }

    delete[] threads;
    delete[] p_thread_info;
}
