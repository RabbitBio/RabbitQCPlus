// Created by ylf9811 on 2021/7/6.
//
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
    //printf("out_block_nums %d\n", out_block_nums);
    out_queue1_ = NULL;
    out_queue2_ = NULL;
    in_is_zip_ = cmd_info1->in_file_name1_.find(".gz") != std::string::npos;
    out_is_zip_ = cmd_info1->out_file_name1_.find(".gz") != std::string::npos;
    if (cmd_info1->write_data_) {
        out_queue1_ = new CIPair[1 << 25];
        queue1P1 = 0;
        queue1P2 = 0;
        queueNumNow1 = 0;
        queueSizeLim1 = 1 << 6;
        if (cmd_info_->interleaved_out_ == 0) {
            out_queue2_ = new CIPair[1 << 25];
            queue2P1 = 0;
            queue2P2 = 0;
            queueNumNow2 = 0;
            queueSizeLim2 = 1 << 6;
        }
        if (out_is_zip_) {
            if (cmd_info1->use_pigz_) {
                pigzQueueNumNow1 = 0;
                pigzQueueSizeLim1 = 1 << 6;
                string out_name1 = cmd_info1->out_file_name1_;
                out_name1 = out_name1.substr(0, out_name1.find(".gz"));
                //out_stream1_ = std::fstream(out_name1, std::ios::out | std::ios::binary);
                out_stream1_.open(out_name1);
                out_stream1_.close();

                pigzQueueNumNow2 = 0;
                pigzQueueSizeLim2 = 1 << 6;
                string out_name2 = cmd_info1->out_file_name2_;
                out_name2 = out_name2.substr(0, out_name2.find(".gz"));
                //out_stream2_ = std::fstream(out_name2, std::ios::out | std::ios::binary);
                out_stream2_.open(out_name2);
                out_stream2_.close();
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
            //out_stream1_ = std::fstream(cmd_info1->out_file_name1_, std::ios::out | std::ios::binary);
            out_stream1_.open(cmd_info1->out_file_name1_);
            if (cmd_info_->interleaved_out_ == 0){
                //out_stream2_ = std::fstream(cmd_info1->out_file_name2_, std::ios::out | std::ios::binary);
                out_stream2_.open(cmd_info1->out_file_name2_);
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
        //std::cout << "readed chunk: " << n_chunks << std::endl;
        dq.Push(n_chunks, fqdatachunk);
    }

    dq.SetCompleted();
    delete fqFileReader;
#ifdef Verbose
    std::cout << "file " << file << " has " << n_chunks << " chunks" << std::endl;
    printf("producer cost %.3f\n", GetTime() - t0);
#endif
}


void PeQc::ProducerPeFastqTask(std::string file, std::string file2, rabbit::fq::FastqDataPool *fastqPool,
                               rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> &dq) {
#ifdef Verbose
    double t0 = GetTime();
#endif
    rabbit::fq::FastqFileReader *fqFileReader;
    rabbit::uint32 tmpSize = 1 << 20;
    if (cmd_info_->seq_len_ <= 200) tmpSize = 1 << 14;
    fqFileReader = new rabbit::fq::FastqFileReader(file, *fastqPool, file2, in_is_zip_, tmpSize);
    int n_chunks = 0;


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
        while (true) {
            rabbit::fq::FastqDataPairChunk *fqdatachunk;
            fqdatachunk = fqFileReader->readNextPairChunkParallel();
            if (fqdatachunk == NULL) break;
            n_chunks++;
            dq.Push(n_chunks, fqdatachunk);
        }
    }


    dq.SetCompleted();
    delete fqFileReader;
#ifdef Verbose
    std::cout << "file " << file << " has " << n_chunks << " chunks" << std::endl;
    printf("producer cost %.3f\n", GetTime() - t0);
#endif
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
        int b_size = min(data1.size(), data2.size());
        for (int i = 0; i < b_size; i++) {
            auto item1 = data1[i];
            auto item2 = data2[i];
            auto name1 = std::string((char *) item1.base + item1.pname, item1.lname);
            auto name2 = std::string((char *) item2.base + item2.pname, item2.lname);
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
                while (queueNumNow1 > queueSizeLim1) {
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
                    while (queueNumNow1 > queueSizeLim1) {
#ifdef Verbose
                        //printf("waiting to push a chunk to out queue1\n");
#endif
                        usleep(100);
                    }
                    while (queueNumNow2 > queueSizeLim2) {
#ifdef Verbose
                        //printf("waiting to push a chunk to out queue2\n");
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
        fastqPool->Release(fqdatachunk->left_part);
        fastqPool->Release(fqdatachunk->right_part);
    }
    done_thread_number_++;
}


void PeQc::ConsumerPeInterFastqTask(ThreadInfo *thread_info, rabbit::fq::FastqDataPool *fastqPool,
                                    rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {
    rabbit::int64 id = 0;
    rabbit::fq::FastqDataChunk *fqdatachunk;
    while (dq.Pop(id, fqdatachunk)) {
        std::vector<neoReference> data;
        std::vector<neoReference> pass_data1, pass_data2;
        rabbit::fq::chunkFormat(fqdatachunk, data, true);
        ASSERT(0);
        ASSERT(data.size() % 2 == 0);
        int out_len1 = 0, out_len2 = 0;
        for (int i = 0; i + 2 <= data.size(); i += 2) {
            auto item1 = data[i];
            auto item2 = data[i + 1];
            auto name1 = std::string((char *) item1.base + item1.pname, item1.lname);
            auto name2 = std::string((char *) item2.base + item2.pname, item2.lname);
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
                while (queueNumNow1 > queueSizeLim1) {
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
                    while (queueNumNow1 > queueSizeLim1) {
#ifdef Verbose
                        //printf("waiting to push a chunk to out queue1\n");
#endif
                        usleep(100);
                    }
                    while (queueNumNow2 > queueSizeLim2) {
#ifdef Verbose
                        //printf("waiting to push a chunk to out queue2\n");
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
}


/**
 * @brief a function to write pe data from out_data1 queue to file1
 */
void PeQc::WriteSeFastqTask1() {
#ifdef Verbose
    double t0 = GetTime();
#endif
    int cnt = 0;
    bool overWhile = 0;
    std::pair<char *, int> now;
    while (true) {
        while (queueNumNow1 == 0) {
            if (done_thread_number_ == cmd_info_->thread_number_) {
                overWhile = 1;
                break;
            }
            usleep(100);
        }
        if (overWhile) break;
        now = out_queue1_[queue1P1++];
        queueNumNow1--;
        if (out_is_zip_) {
            if (cmd_info_->use_pigz_) {
                while (pigzQueueNumNow1 > pigzQueueSizeLim1) {

#ifdef Verbose
                    //printf("waiting to push a chunk to pigz queue1\n");
#endif
                    usleep(100);
                }
                pigzQueue1->enqueue(now);
                pigzQueueNumNow1++;
            } else {
                int written = gzwrite(zip_out_stream1, now.first, now.second);
                if (written != now.second) {
                    printf("gzwrite error\n");
                    exit(0);
                }
                delete[] now.first;
            }
        } else {
            out_stream1_.write(now.first, now.second);
            delete[] now.first;
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
        out_stream1_.close();
    }
#ifdef Verbose
    printf("write1 cost %.5f\n", GetTime() - t0);
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
    bool overWhile = 0;
    std::pair<char *, int> now;
    while (true) {
        while (queueNumNow2 == 0) {
            if (done_thread_number_ == cmd_info_->thread_number_) {
                overWhile = 1;
                break;
            }
            usleep(100);
        }
        if (overWhile) break;
        now = out_queue2_[queue2P1++];
        queueNumNow2--;
        if (out_is_zip_) {
            if (cmd_info_->use_pigz_) {
                while (pigzQueueNumNow2 > pigzQueueSizeLim2) {

#ifdef Verbose
                    //printf("waiting to push a chunk to pigz queue2\n");
#endif
                    usleep(1000000);
                }
                pigzQueue2->enqueue(now);
                pigzQueueNumNow2++;
            } else {
                int written = gzwrite(zip_out_stream2, now.first, now.second);
                if (written != now.second) {
                    printf("GG");
                    exit(0);
                }

                delete[] now.first;
            }
        } else {

            out_stream2_.write(now.first, now.second);
            delete[] now.first;
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
        out_stream2_.close();
    }
#ifdef Verbose
    printf("write2 cost %.5f\n", GetTime() - t0);
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
    infos[4] = "-4";
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
    infos[4] = "-4";
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

/**
 * @brief do QC for pair-end data
 */
void PeQc::ProcessPeFastq() {
    if (cmd_info_->interleaved_in_) {
        auto *fastqPool = new rabbit::fq::FastqDataPool(128, 1 << 22);
        rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> queue1(64, 1);
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
#ifdef Verbose
        printf("all thrad done\n");
        printf("now merge thread info\n");
#endif
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

#ifdef Verbose
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
                if (merge_insert_size[i] > merge_insert_size[mx_id]) mx_id = i;
            }
            //printf("Insert size peak (evaluated by paired-end reads): %d\n", mx_id);
            printf("Insert size peak (based on PE overlap analyze): %d\n", mx_id);
        }
        std::string srr_name1 = cmd_info_->in_file_name1_;
        srr_name1 = PaseFileName(srr_name1);
        std::string srr_name2 = cmd_info_->in_file_name2_;
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

        thread *pugzer1;
        thread *pugzer2;

        if (cmd_info_->use_pugz_) {
            pugzer1 = new thread(bind(&::PeQc::PugzTask1, this));
            pugzer2 = new thread(bind(&::PeQc::PugzTask2, this));
        }


        auto *fastqPool = new rabbit::fq::FastqDataPool(128, 1 << 22);
        rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> queue1(64, 1);
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
        thread *pigzer1;
        thread *pigzer2;
        if (cmd_info_->use_pigz_) {
            pigzer1 = new thread(bind(&PeQc::PigzTask1, this));
            if (cmd_info_->interleaved_out_ == 0)
                pigzer2 = new thread(bind(&PeQc::PigzTask2, this));
        }
        std::thread producer(
                std::bind(&PeQc::ProducerPeFastqTask, this, cmd_info_->in_file_name1_, cmd_info_->in_file_name2_,
                          fastqPool, std::ref(queue1)));
        auto **threads = new std::thread *[cmd_info_->thread_number_];
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            threads[t] = new std::thread(
                    std::bind(&PeQc::ConsumerPeFastqTask, this, p_thread_info[t], fastqPool, std::ref(queue1)));
        }

        if (cmd_info_->use_pugz_) {
            pugzer1->join();
            pugzer2->join();
        }

        producer.join();
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            threads[t]->join();
        }
        if (cmd_info_->write_data_) {
            write_thread1->join();
            writerDone1 = 1;
            if (cmd_info_->interleaved_out_ == 0) {

                write_thread2->join();
                writerDone2 = 1;
            }
        }

        if (cmd_info_->use_pigz_) {
            pigzer1->join();
            if (cmd_info_->interleaved_out_ == 0)
                pigzer2->join();
        }
#ifdef Verbose
        printf("all thrad done\n");
        printf("now merge thread info\n");
#endif
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
                if (merge_insert_size[i] > merge_insert_size[mx_id]) mx_id = i;
            }
            //printf("Insert size peak (evaluated by paired-end reads): %d\n", mx_id);
            printf("Insert size peak (based on PE overlap analyze): %d\n", mx_id);
        }
        std::string srr_name1 = cmd_info_->in_file_name1_;
        srr_name1 = PaseFileName(srr_name1);
        std::string srr_name2 = cmd_info_->in_file_name2_;
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
