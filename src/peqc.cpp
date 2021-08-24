//
// Created by ylf9811 on 2021/7/6.
//

#include "peqc.h"

/**
 * @brief Construct function
 * @param cmd_info1 : cmd information
 */
PeQc::PeQc(CmdInfo *cmd_info1) {
    cmd_info_ = cmd_info1;
    filter_ = new Filter(cmd_info1);
    done_thread_number_ = 0;
    int out_block_nums = int(1.0 * cmd_info1->in_file_size1_ / cmd_info1->out_block_size_);
    printf("out_block_nums %d\n", out_block_nums);
    out_queue1_ = NULL;
    out_queue2_ = NULL;
    in_is_zip_ = cmd_info1->in_file_name1_.find(".gz") != std::string::npos;
    out_is_zip_ = cmd_info1->out_file_name1_.find(".gz") != std::string::npos;

    if (cmd_info1->write_data_) {
        out_queue1_ = new moodycamel::ConcurrentQueue<std::pair<char *, int>>
                (out_block_nums + 1000);
        if (cmd_info_->interleaved_out_ == 0)
            out_queue2_ = new moodycamel::ConcurrentQueue<std::pair<char *, int>>
                    (out_block_nums + 1000);

        if (out_is_zip_) {
            printf("open gzip stream1 %s\n", cmd_info1->out_file_name1_.c_str());
            printf("open gzip stream2 %s\n", cmd_info1->out_file_name2_.c_str());
            zip_out_stream1 = gzopen(cmd_info1->out_file_name1_.c_str(), "w");
            gzsetparams(zip_out_stream1, cmd_info1->compression_level_, Z_DEFAULT_STRATEGY);
            gzbuffer(zip_out_stream1, 1024 * 1024);
            zip_out_stream2 = gzopen(cmd_info1->out_file_name2_.c_str(), "w");
            gzsetparams(zip_out_stream2, cmd_info1->compression_level_, Z_DEFAULT_STRATEGY);
            gzbuffer(zip_out_stream2, 1024 * 1024);
        } else {
            printf("open stream1 %s\n", cmd_info1->out_file_name1_.c_str());
            if (cmd_info_->interleaved_out_ == 0)
                printf("open stream2 %s\n", cmd_info1->out_file_name2_.c_str());
            out_stream1_ = std::fstream(cmd_info1->out_file_name1_, std::ios::out | std::ios::binary);
            if (cmd_info_->interleaved_out_ == 0)
                out_stream2_ = std::fstream(cmd_info1->out_file_name2_, std::ios::out | std::ios::binary);
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
}

PeQc::~PeQc() {
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

//
//void PeQc::PrintRead(neoReference &ref) {
//    std::cout << std::string((char *) ref.base + ref.pname, ref.lname) << std::endl;
//    std::cout << std::string((char *) ref.base + ref.pseq, ref.lseq) << std::endl;
//    std::cout << std::string((char *) ref.base + ref.pstrand, ref.lstrand) << std::endl;
//    std::cout << std::string((char *) ref.base + ref.pqual, ref.lqual) << std::endl;
//}

std::string PeQc::Read2String(neoReference &ref) {
    return std::string((char *) ref.base + ref.pname, ref.lname) + "\n" +
           std::string((char *) ref.base + ref.pseq, ref.lseq) + "\n" +
           std::string((char *) ref.base + ref.pstrand, ref.lstrand) + "\n" +
           std::string((char *) ref.base + ref.pqual, ref.lqual) + "\n";
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

void PeQc::ProducerPeInterFastqTask(std::string file, rabbit::fq::FastqDataPool *fastq_data_pool,
                                    rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {
    double t0 = GetTime();

    rabbit::fq::FastqFileReader *fqFileReader;
    fqFileReader = new rabbit::fq::FastqFileReader(file, fastq_data_pool, "", in_is_zip_);
    int64_t n_chunks = 0;
    while (true) {
        rabbit::fq::FastqDataChunk *fqdatachunk;
        fqdatachunk = fqFileReader->readNextChunk();
        if (fqdatachunk == NULL) break;
        n_chunks++;
        //std::cout << "readed chunk: " << n_chunks << std::endl;
        dq.Push(n_chunks, fqdatachunk);
    }

    dq.SetCompleted();
    delete fqFileReader;
    std::cout << "file " << file << " has " << n_chunks << " chunks" << std::endl;
    printf("producer cost %.3f\n", GetTime() - t0);

}


void PeQc::ProducerPeFastqTask(std::string file, std::string file2, rabbit::fq::FastqDataPool *fastqPool,
                               rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> &dq) {
    double t0 = GetTime();

    rabbit::fq::FastqFileReader *fqFileReader;
    fqFileReader = new rabbit::fq::FastqFileReader(file, fastqPool, file2, in_is_zip_);
    int n_chunks = 0;

    rabbit::fq::FastqDataPairChunk *fqdatachunk;
    while (true) {
        fqdatachunk = fqFileReader->readNextPairChunk();
        if (fqdatachunk == NULL) break;
        n_chunks++;
        dq.Push(n_chunks, fqdatachunk);
    }

    dq.SetCompleted();
    delete fqFileReader;
    std::cout << "file " << file << " has " << n_chunks << " chunks" << std::endl;
    printf("producer cost %.3f\n", GetTime() - t0);
}

/**
 * @brief get fastq data chunks from the data queue and do QC for them
 * @param thread_info : thread information
 * @param fastq_data_pool :a fastq data pool, it will be used to release data chunk
 * @param dq : data queue
 */
void PeQc::ConsumerPeFastqTask(ThreadInfo *thread_info, rabbit::fq::FastqDataPool *fastqPool,
                               rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> &dq) {
    rabbit::int64 id = 0;
    rabbit::fq::FastqDataPairChunk *fqdatachunk;
    while (dq.Pop(id, fqdatachunk)) {
        std::vector<neoReference> data1, data2;
        std::vector<neoReference> pass_data1, pass_data2;
        rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk *) (fqdatachunk->left_part), data1, true);
        rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk *) (fqdatachunk->right_part), data2, true);
        ASSERT(data1.size() == data2.size());
        int out_len1 = 0, out_len2 = 0;
        for (int i = 0; i < data1.size(); i++) {
            auto item1 = data1[i];
            auto item2 = data2[i];
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
                Adapter::CorrectData(item1, item2, overlap_res);
            }
            if (trim_res1 && trim_res2 && cmd_info_->trim_adapter_) {
                bool trimmed = Adapter::TrimAdapter(item1, item2, overlap_res.offset, overlap_res.overlap_len);
                if (!trimmed) {
                    if (cmd_info_->detect_adapter1_)
                        Adapter::TrimAdapter(item1, cmd_info_->adapter_seq1_, false);
                    if (cmd_info_->detect_adapter2_)
                        Adapter::TrimAdapter(item2, cmd_info_->adapter_seq2_, true);
                }
            }


            //do filer in refs
            int filter_res1 = filter_->ReadFiltering(item1, trim_res1);
            int filter_res2 = filter_->ReadFiltering(item2, trim_res2);

            if (filter_res1 == 0 && filter_res2 == 0) {
                thread_info->aft_state1_->StateInfo(item1);
                thread_info->aft_state2_->StateInfo(item2);
                if (cmd_info_->write_data_) {
                    pass_data1.push_back(item1);
                    pass_data2.push_back(item2);
                    out_len1 += item1.lname + item1.lseq + item1.lstrand + item1.lqual + 4;
                    out_len2 += item2.lname + item2.lseq + item2.lstrand + item2.lqual + 4;
                }
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
                out_queue1_->enqueue({out_data, pos});
            } else {
                if (pass_data1.size() > 0) {
                    char *out_data1 = new char[out_len1];
                    int pos = 0;
                    for (auto item:pass_data1) {
                        Read2Chars(item, out_data1, pos);
                    }
                    ASSERT(pos == out_len1);
                    out_queue1_->enqueue({out_data1, out_len1});
                }
                if (pass_data2.size() > 0) {
                    char *out_data2 = new char[out_len2];
                    int pos = 0;
                    for (auto item:pass_data2) {
                        Read2Chars(item, out_data2, pos);
                    }
                    ASSERT(pos == out_len2);
                    out_queue2_->enqueue({out_data2, out_len2});
                }
            }

        }
        fastqPool->Release(fqdatachunk->left_part);
        fastqPool->Release(fqdatachunk->right_part);
    }
    done_thread_number_++;
}

void PeQc::ConsumerPeInterFastqTask(ThreadInfo *thread_info, rabbit::fq::FastqDataPool *fastqPool,
                                    rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {
    rabbit::int64 id = 0;
    rabbit::fq::FastqDataChunk *fqdatachunk;
    rabbit::fq::FastqDataChunk *laterFqdatachunk;
    int now = 0;
    neoReference item1, item2;
    bool laterRelease = 0;

    while (dq.Pop(id, fqdatachunk)) {
        std::vector<neoReference> data;
        std::vector<neoReference> pass_data1, pass_data2;
        rabbit::fq::chunkFormat(fqdatachunk, data, true);

        int out_len1 = 0, out_len2 = 0;
        for (int i = 0; i < data.size(); i++) {
            if (now == 0) {
                item1 = data[i];
                now ^= 1;
                continue;
            }
            item2 = data[i];
            now ^= 1;
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
                Adapter::CorrectData(item1, item2, overlap_res);
            }
            if (trim_res1 && trim_res2 && cmd_info_->trim_adapter_) {
                bool trimmed = Adapter::TrimAdapter(item1, item2, overlap_res.offset, overlap_res.overlap_len);
                if (!trimmed) {
                    if (cmd_info_->detect_adapter1_)
                        Adapter::TrimAdapter(item1, cmd_info_->adapter_seq1_, false);
                    if (cmd_info_->detect_adapter2_)
                        Adapter::TrimAdapter(item2, cmd_info_->adapter_seq2_, true);
                }
            }


            //do filer in refs
            int filter_res1 = filter_->ReadFiltering(item1, trim_res1);
            int filter_res2 = filter_->ReadFiltering(item2, trim_res2);

            if (filter_res1 == 0 && filter_res2 == 0) {
                thread_info->aft_state1_->StateInfo(item1);
                thread_info->aft_state2_->StateInfo(item2);
                if (cmd_info_->write_data_) {
                    pass_data1.push_back(item1);
                    pass_data2.push_back(item2);
                    out_len1 += item1.lname + item1.lseq + item1.lstrand + item1.lqual + 4;
                    out_len2 += item2.lname + item2.lseq + item2.lstrand + item2.lqual + 4;
                }
            }
        }
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
            out_queue1_->enqueue({out_data, pos});
        } else {
            if (pass_data1.size() > 0) {
                char *out_data1 = new char[out_len1];
                int pos = 0;
                for (auto item:pass_data1) {
                    Read2Chars(item, out_data1, pos);
                }
                ASSERT(pos == out_len1);
                out_queue1_->enqueue({out_data1, out_len1});
            }
            if (pass_data2.size() > 0) {
                char *out_data2 = new char[out_len2];
                int pos = 0;
                for (auto item:pass_data2) {
                    Read2Chars(item, out_data2, pos);
                }
                ASSERT(pos == out_len2);
                out_queue2_->enqueue({out_data2, out_len2});
            }
        }
        if (laterRelease) {
            fastqPool->Release(laterFqdatachunk);
        }
        if (now == 0) {
            fastqPool->Release(fqdatachunk);
            laterRelease = 0;
        } else {
            laterRelease = 1;
            laterFqdatachunk = fqdatachunk;
        }

    }
    done_thread_number_++;
}


/**
 * @brief a function to write pe data from out_data1 queue to file1
 */
void PeQc::WriteSeFastqTask1() {
    double t0 = GetTime();
    int cnt = 0;
    while (true) {
        if (out_queue1_->size_approx() == 0 && done_thread_number_ == cmd_info_->thread_number_) {
            break;
        }
        if (out_queue1_->size_approx() == 0) {
            usleep(100);
        }
        std::pair<char *, int> now;
        while (out_queue1_->size_approx()) {
            out_queue1_->try_dequeue(now);
//            printf("write thread working, write %d %d\n", now.second, cnt++);
            if (out_is_zip_) {
                int written = gzwrite(zip_out_stream1, now.first, now.second);
                if (written != now.second) {
                    printf("GG");
                    exit(0);
                }

            } else {
                out_stream1_.write(now.first, now.second);
            }
            delete[] now.first;
        }
    }
    if (out_is_zip_) {
        if (zip_out_stream1) {
            gzflush(zip_out_stream1, Z_FINISH);
            gzclose(zip_out_stream1);
            zip_out_stream1 = NULL;
        }
    } else {
        out_stream1_.close();
    }
    printf("write cost %.5f\n", GetTime() - t0);
}

/**
 * @brief a function to write pe data from out_data2 queue to file2
 */
void PeQc::WriteSeFastqTask2() {
    double t0 = GetTime();
    int cnt = 0;
    while (true) {
        if (out_queue2_->size_approx() == 0 && done_thread_number_ == cmd_info_->thread_number_) {
            break;
        }
        if (out_queue2_->size_approx() == 0) {
            usleep(100);
        }
        std::pair<char *, int> now;
        while (out_queue2_->size_approx()) {
            out_queue2_->try_dequeue(now);
//            printf("write thread working, write %d %d\n", now.second, cnt++);
            if (out_is_zip_) {
                int written = gzwrite(zip_out_stream2, now.first, now.second);
                if (written != now.second) {
                    printf("GG");
                    exit(0);
                }

            } else {
                out_stream2_.write(now.first, now.second);
            }
            delete[] now.first;
        }
    }
    if (out_is_zip_) {
        if (zip_out_stream2) {
            gzflush(zip_out_stream2, Z_FINISH);
            gzclose(zip_out_stream2);
            zip_out_stream2 = NULL;
        }
    } else {
        out_stream2_.close();
    }
    printf("write cost %.5f\n", GetTime() - t0);

}

/**
 * @brief do QC for pair-end data
 */
void PeQc::ProcessPeFastq() {
    if (cmd_info_->interleaved_in_) {
        auto *fastqPool = new rabbit::fq::FastqDataPool(256, 1 << 22);
        //TODO replace this queue
        rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> queue1(128, 1);
        auto **p_thread_info = new ThreadInfo *[cmd_info_->thread_number_];
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            p_thread_info[t] = new ThreadInfo(cmd_info_, true);
        }
        std::thread *write_thread1;
        std::thread *write_thread2;
        if (cmd_info_->write_data_) {
            write_thread1 = new std::thread(std::bind(&PeQc::WriteSeFastqTask1, this));
            if (cmd_info_->interleaved_out_ == 0)
                write_thread2 = new std::thread(std::bind(&PeQc::WriteSeFastqTask2, this));
        }

        //TODO bind ?
        std::thread producer(
                std::bind(&PeQc::ProducerPeInterFastqTask, this, cmd_info_->in_file_name1_, fastqPool,
                          std::ref(queue1)));
        auto **threads = new std::thread *[cmd_info_->thread_number_];
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            threads[t] = new std::thread(
                    std::bind(&PeQc::ConsumerPeInterFastqTask, this, p_thread_info[t], fastqPool, std::ref(queue1)));
        }
        producer.join();
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            threads[t]->join();
        }
        if (cmd_info_->write_data_) {
            write_thread1->join();
            if (cmd_info_->interleaved_out_ == 0)
                write_thread2->join();
        }

        printf("all thrad done\n");
        printf("now merge thread info\n");

        std::vector<State *> pre_vec_state1;
        std::vector<State *> pre_vec_state2;
        std::vector<State *> aft_vec_state1;
        std::vector<State *> aft_vec_state2;

        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            pre_vec_state1.push_back(p_thread_info[t]->pre_state1_);
            pre_vec_state2.push_back(p_thread_info[t]->pre_state2_);
            aft_vec_state1.push_back(p_thread_info[t]->aft_state1_);
            aft_vec_state2.push_back(p_thread_info[t]->aft_state2_);
        }
        auto pre_state1 = State::MergeStates(pre_vec_state1);
        auto pre_state2 = State::MergeStates(pre_vec_state2);
        auto aft_state1 = State::MergeStates(aft_vec_state1);
        auto aft_state2 = State::MergeStates(aft_vec_state2);


        printf("merge done\n");
        printf("print pre state1 info :\n\n");
        State::PrintStates(pre_state1);
        printf("print pre state2 info :\n\n");
        State::PrintStates(pre_state2);
        printf("print aft state1 info :\n\n");
        State::PrintStates(aft_state1);
        printf("print aft state2 info :\n\n");
        State::PrintStates(aft_state2);

        if (cmd_info_->do_overrepresentation_) {

            auto pre_hash_graph1 = pre_state1->GetHashGraph();
            int pre_hash_num1 = pre_state1->GetHashNum();
            auto pre_hash_graph2 = pre_state2->GetHashGraph();
            int pre_hash_num2 = pre_state2->GetHashNum();

            auto aft_hash_graph1 = aft_state1->GetHashGraph();
            int aft_hash_num1 = aft_state1->GetHashNum();
            auto aft_hash_graph2 = aft_state2->GetHashGraph();
            int aft_hash_num2 = aft_state2->GetHashNum();
            ofstream ofs;
            ofs.open("ORP2.log", ifstream::out);
            for (int i = 0; i < pre_hash_num1; i++) {
                ofs << pre_hash_graph1[i].seq << " " << pre_hash_graph1[i].cnt << "\n";
            }
            for (int i = 0; i < pre_hash_num2; i++) {
                ofs << pre_hash_graph2[i].seq << " " << pre_hash_graph2[i].cnt << "\n";
            }
            ofs.close();
            ofs.open("ORP3.log", ifstream::out);
            for (int i = 0; i < aft_hash_num1; i++) {
                ofs << aft_hash_graph1[i].seq << " " << aft_hash_graph1[i].cnt << "\n";
            }
            for (int i = 0; i < aft_hash_num2; i++) {
                ofs << aft_hash_graph2[i].seq << " " << aft_hash_graph2[i].cnt << "\n";
            }
            ofs.close();
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

        int64_t *merge_insert_size;
        if (!cmd_info_->no_insert_size_) {
            merge_insert_size = new int64_t[cmd_info_->max_insert_size_ + 1];
            memset(merge_insert_size, 0, sizeof(int64_t) * (cmd_info_->max_insert_size_ + 1));

            for (int t = 0; t < cmd_info_->thread_number_; t++) {
                for (int i = 0; i <= cmd_info_->max_insert_size_; i++) {
                    merge_insert_size[i] += p_thread_info[t]->insert_size_dist_[i];
                }
            }
            int mx_id = 0;
            for (int i = 0; i < cmd_info_->max_insert_size_; i++) {
                if (merge_insert_size[i] > merge_insert_size[mx_id])mx_id = i;
            }
            printf("Insert size peak (evaluated by paired-end reads): %d\n", mx_id);
        }

        Repoter::ReportHtmlPe(pre_state1, pre_state2, aft_state1, aft_state2, cmd_info_->in_file_name1_,
                              cmd_info_->in_file_name2_, dupRate * 100.0, merge_insert_size);

        printf("report done\n");

        delete pre_state1;
        delete pre_state2;
        delete aft_state1;
        delete aft_state2;
        if (!cmd_info_->no_insert_size_)
            delete[] merge_insert_size;


        delete fastqPool;
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            delete p_thread_info[t];
            delete threads[t];
        }
        delete[] threads;
        delete[] p_thread_info;
        if (cmd_info_->write_data_) {
            delete write_thread1;
            if (cmd_info_->interleaved_out_ == 0)
                delete write_thread2;
        }
    } else {
        auto *fastqPool = new rabbit::fq::FastqDataPool(256, 1 << 22);
        //TODO replace this queue
        rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> queue1(128, 1);
        auto **p_thread_info = new ThreadInfo *[cmd_info_->thread_number_];
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            p_thread_info[t] = new ThreadInfo(cmd_info_, true);
        }
        std::thread *write_thread1;
        std::thread *write_thread2;
        if (cmd_info_->write_data_) {
            write_thread1 = new std::thread(std::bind(&PeQc::WriteSeFastqTask1, this));
            if (cmd_info_->interleaved_out_ == 0)
                write_thread2 = new std::thread(std::bind(&PeQc::WriteSeFastqTask2, this));
        }

        //TODO bind ?
        std::thread producer(
                std::bind(&PeQc::ProducerPeFastqTask, this, cmd_info_->in_file_name1_, cmd_info_->in_file_name2_,
                          fastqPool, std::ref(queue1)));
        auto **threads = new std::thread *[cmd_info_->thread_number_];
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            threads[t] = new std::thread(
                    std::bind(&PeQc::ConsumerPeFastqTask, this, p_thread_info[t], fastqPool, std::ref(queue1)));
        }
        producer.join();
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            threads[t]->join();
        }
        if (cmd_info_->write_data_) {
            write_thread1->join();
            if (cmd_info_->interleaved_out_ == 0)
                write_thread2->join();
        }

        printf("all thrad done\n");
        printf("now merge thread info\n");

        std::vector<State *> pre_vec_state1;
        std::vector<State *> pre_vec_state2;
        std::vector<State *> aft_vec_state1;
        std::vector<State *> aft_vec_state2;

        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            pre_vec_state1.push_back(p_thread_info[t]->pre_state1_);
            pre_vec_state2.push_back(p_thread_info[t]->pre_state2_);
            aft_vec_state1.push_back(p_thread_info[t]->aft_state1_);
            aft_vec_state2.push_back(p_thread_info[t]->aft_state2_);
        }
        auto pre_state1 = State::MergeStates(pre_vec_state1);
        auto pre_state2 = State::MergeStates(pre_vec_state2);
        auto aft_state1 = State::MergeStates(aft_vec_state1);
        auto aft_state2 = State::MergeStates(aft_vec_state2);


        printf("merge done\n");
        printf("print pre state1 info :\n\n");
        State::PrintStates(pre_state1);
        printf("print pre state2 info :\n\n");
        State::PrintStates(pre_state2);
        printf("print aft state1 info :\n\n");
        State::PrintStates(aft_state1);
        printf("print aft state2 info :\n\n");
        State::PrintStates(aft_state2);

        if (cmd_info_->do_overrepresentation_) {

            auto pre_hash_graph1 = pre_state1->GetHashGraph();
            int pre_hash_num1 = pre_state1->GetHashNum();
            auto pre_hash_graph2 = pre_state2->GetHashGraph();
            int pre_hash_num2 = pre_state2->GetHashNum();

            auto aft_hash_graph1 = aft_state1->GetHashGraph();
            int aft_hash_num1 = aft_state1->GetHashNum();
            auto aft_hash_graph2 = aft_state2->GetHashGraph();
            int aft_hash_num2 = aft_state2->GetHashNum();
            ofstream ofs;
            ofs.open("ORP2.log", ifstream::out);
            for (int i = 0; i < pre_hash_num1; i++) {
                ofs << pre_hash_graph1[i].seq << " " << pre_hash_graph1[i].cnt << "\n";
            }
            for (int i = 0; i < pre_hash_num2; i++) {
                ofs << pre_hash_graph2[i].seq << " " << pre_hash_graph2[i].cnt << "\n";
            }
            ofs.close();
            ofs.open("ORP3.log", ifstream::out);
            for (int i = 0; i < aft_hash_num1; i++) {
                ofs << aft_hash_graph1[i].seq << " " << aft_hash_graph1[i].cnt << "\n";
            }
            for (int i = 0; i < aft_hash_num2; i++) {
                ofs << aft_hash_graph2[i].seq << " " << aft_hash_graph2[i].cnt << "\n";
            }
            ofs.close();
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

        int64_t *merge_insert_size;
        if (!cmd_info_->no_insert_size_) {
            merge_insert_size = new int64_t[cmd_info_->max_insert_size_ + 1];
            memset(merge_insert_size, 0, sizeof(int64_t) * (cmd_info_->max_insert_size_ + 1));

            for (int t = 0; t < cmd_info_->thread_number_; t++) {
                for (int i = 0; i <= cmd_info_->max_insert_size_; i++) {
                    merge_insert_size[i] += p_thread_info[t]->insert_size_dist_[i];
                }
            }
            int mx_id = 0;
            for (int i = 0; i < cmd_info_->max_insert_size_; i++) {
                if (merge_insert_size[i] > merge_insert_size[mx_id])mx_id = i;
            }
            printf("Insert size peak (evaluated by paired-end reads): %d\n", mx_id);
        }

        Repoter::ReportHtmlPe(pre_state1, pre_state2, aft_state1, aft_state2, cmd_info_->in_file_name1_,
                              cmd_info_->in_file_name2_, dupRate * 100.0, merge_insert_size);

        printf("report done\n");

        delete pre_state1;
        delete pre_state2;
        delete aft_state1;
        delete aft_state2;
        if (!cmd_info_->no_insert_size_)
            delete[] merge_insert_size;


        delete fastqPool;
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            delete p_thread_info[t];
            delete threads[t];
        }
        delete[] threads;
        delete[] p_thread_info;
        if (cmd_info_->write_data_) {
            delete write_thread1;
            if (cmd_info_->interleaved_out_ == 0)
                delete write_thread2;
        }
    }


}
