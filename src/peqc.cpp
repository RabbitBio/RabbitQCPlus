// Created by ylf9811 on 2021/7/6.
//
//

#include "peqc.h"
using namespace std;

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
    nowChunkId = 1;
    in_is_zip_ = cmd_info1->in_file_name1_.find(".gz") != string::npos;
    out_is_zip_ = cmd_info1->out_file_name1_.find(".gz") != string::npos;
    if (cmd_info1->write_data_) {
        out_queue1_ = new CIPair[1 << 20];
        queue1P1 = 0;
        queue1P2 = 0;
        queueNumNow1 = 0;
        queueSizeLim1 = 1 << 5;
        if (cmd_info_->interleaved_out_ == 0) {
            out_queue2_ = new CIPair[1 << 20];
            queue2P1 = 0;
            queue2P2 = 0;
            queueNumNow2 = 0;
            queueSizeLim2 = 1 << 5;
        }
        if (out_is_zip_) {
            if (cmd_info1->use_pigz_) {
                pigzQueueNumNow1 = 0;
                pigzQueueSizeLim1 = 1 << 5;
                string out_name1 = cmd_info1->out_file_name1_;
                out_name1 = out_name1.substr(0, out_name1.find(".gz"));
                out_stream1_.open(out_name1);
                out_stream1_.close();

                pigzQueueNumNow2 = 0;
                pigzQueueSizeLim2 = 1 << 5;
                string out_name2 = cmd_info1->out_file_name2_;
                out_name2 = out_name2.substr(0, out_name2.find(".gz"));
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
            out_stream1_.open(cmd_info1->out_file_name1_);
            if (cmd_info_->interleaved_out_ == 0){
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
    if(cmd_info1->do_correction_with_care_) {
        careQueue1 = new moodycamel::ReaderWriterQueue<pair<char *, int>>(1 << 20);
        careQueue2 = new moodycamel::ReaderWriterQueue<pair<char *, int>>(1 << 20);
    }
    changeNum = 0;
    careStartWrite = 0;
    careDone1 = 0;
    careDone2 = 0;
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

void PeQc::careProcess() {
    
    printf("pairmode %s\n", cmd_info_->pairmode_.c_str());
    printf("coverage %d\n", cmd_info_->coverage_);
    vector<char*>paras;
    paras.push_back("./RabbitQCPlus");
    paras.push_back("-i");
    paras.push_back((char*)(cmd_info_->in_file_name1_.data()));
    paras.push_back("-i");
    paras.push_back((char*)(cmd_info_->in_file_name2_.data()));

    paras.push_back("-d");
    paras.push_back("./");

    paras.push_back("-o");
    paras.push_back("tmp1.fq");
    
    paras.push_back("-o");
    paras.push_back("tmp2.fq");

    paras.push_back("-c");
    string str_coverage = to_string(cmd_info_->coverage_);
    printf("str_coverage %s\n", str_coverage.c_str());
    paras.push_back((char*)(str_coverage.data()));

    paras.push_back("-t");
    string str_thread = to_string(cmd_info_->correct_threadnum_);
    printf("str_thread %s\n", str_thread.c_str());
    paras.push_back((char*)(str_thread.data()));

    paras.push_back("--pairmode");
    if(cmd_info_->pairmode_ == "SE") paras.push_back("SE");
    else paras.push_back("PE");


    string str_hashmaps = to_string(cmd_info_->hashmaps_);
    if(cmd_info_->hashmaps_ != 48) {
        paras.push_back("--hashmaps");
        paras.push_back((char*)(str_hashmaps.data()));
    }
    string str_kmerlength = to_string(cmd_info_->kmerlength_);
    if(cmd_info_->kmerlength_ != 0) {
        paras.push_back("--kmerlength");
        paras.push_back((char*)(str_kmerlength.data()));
    }
    if(cmd_info_->enforceHashmapCount_ != false) {
        paras.push_back("--enforceHashmapCount");
    }
    if(cmd_info_->useQualityScores_ != false) {
        paras.push_back("--useQualityScores");
    }
    string str_qualityScoreBits = to_string(cmd_info_->qualityScoreBits_);
    if(cmd_info_->qualityScoreBits_ != 8) {
        paras.push_back("--qualityScoreBits");
        paras.push_back((char*)(str_qualityScoreBits.data()));
    }
    if(cmd_info_->excludeAmbiguous_ != false) {
        paras.push_back("--excludeAmbiguous");
    }
    string str_maxmismatchratio = to_string(cmd_info_->maxmismatchratio_);
    if(cmd_info_->maxmismatchratio_ != 0.200000) {
        paras.push_back("--maxmismatchratio");
        paras.push_back((char*)(str_maxmismatchratio.data()));
    }
    string str_minalignmentoverlap = to_string(cmd_info_->minalignmentoverlap_);
    if(cmd_info_->minalignmentoverlap_ != 30) {
        paras.push_back("--minalignmentoverlap");
        paras.push_back((char*)(str_minalignmentoverlap.data()));
    }
    string str_minalignmentoverlapratio = to_string(cmd_info_->minalignmentoverlapratio_);
    if(cmd_info_->minalignmentoverlapratio_ != 0.300000) {
        paras.push_back("--minalignmentoverlapratio");
        paras.push_back((char*)(str_minalignmentoverlapratio.data()));
    }
    string str_errorfactortuning = to_string(cmd_info_->errorfactortuning_);
    if(cmd_info_->errorfactortuning_ != 0.060000) {
        paras.push_back("--errorfactortuning");
        paras.push_back((char*)(str_errorfactortuning.data()));
    }
    string str_coveragefactortuning = to_string(cmd_info_->coveragefactortuning_);
    if(cmd_info_->coveragefactortuning_ != 0.600000) {
        paras.push_back("--coveragefactortuning");
        paras.push_back((char*)(str_coveragefactortuning.data()));
    }
    if(cmd_info_->showProgress_ != false) {
        paras.push_back("--showProgress");
    }
    string str_tempdir = cmd_info_->tempdir_;
    if(cmd_info_->tempdir_ != "") {
        paras.push_back("--tempdir");
        paras.push_back((char*)(str_tempdir.data()));
    }
    string str_save_preprocessedreads_to = cmd_info_->save_preprocessedreads_to_;
    if(cmd_info_->save_preprocessedreads_to_ != "") {
        paras.push_back("--save_preprocessedreads_to");
        paras.push_back((char*)(str_save_preprocessedreads_to.data()));
    }
    string str_load_preprocessedreads_from = cmd_info_->load_preprocessedreads_from_;
    if(cmd_info_->load_preprocessedreads_from_ != "") {
        paras.push_back("--load_preprocessedreads_from");
        paras.push_back((char*)(str_load_preprocessedreads_from.data()));
    }
    string str_save_hashtables_to = cmd_info_->save_hashtables_to_;
    if(cmd_info_->save_hashtables_to_ != "") {
        paras.push_back("--save_hashtables_to");
        paras.push_back((char*)(str_save_hashtables_to.data()));
    }
    string str_load_hashtables_from = cmd_info_->load_hashtables_from_;
    if(cmd_info_->load_hashtables_from_ != "") {
        paras.push_back("--load_hashtables_from");
        paras.push_back((char*)(str_load_hashtables_from.data()));
    }
    string str_memHashtables = cmd_info_->memHashtables_;
    if(cmd_info_->memHashtables_ != "") {
        paras.push_back("--memHashtables");
        paras.push_back((char*)(str_memHashtables.data()));
    }
    string str_memTotal = cmd_info_->memTotal_;
    if(cmd_info_->memTotal_ != "") {
        paras.push_back("--memTotal");
        paras.push_back((char*)(str_memTotal.data()));
    }
    string str_hashloadfactor = to_string(cmd_info_->hashloadfactor_);
    if(cmd_info_->hashloadfactor_ != 0.800000) {
        paras.push_back("--hashloadfactor");
        paras.push_back((char*)(str_hashloadfactor.data()));
    }
    string str_fixedNumberOfReads = to_string(cmd_info_->fixedNumberOfReads_);
    if(cmd_info_->fixedNumberOfReads_ != 0) {
        paras.push_back("--fixedNumberOfReads");
        paras.push_back((char*)(str_fixedNumberOfReads.data()));
    }
    if(cmd_info_-> singlehash_!= false) {
        paras.push_back("--singlehash");
    }
    if(cmd_info_->correctionQualityLabels_ != false) {
        paras.push_back("--correctionQualityLabels");
    }
    if(cmd_info_->candidateCorrection_ != false) {
        paras.push_back("--candidateCorrection");
    }
    string str_candidateCorrectionNewColumns = to_string(cmd_info_->candidateCorrectionNewColumns_);
    if(cmd_info_->candidateCorrectionNewColumns_ != 15) {
        paras.push_back("--candidateCorrectionNewColumns");
        paras.push_back((char*)(str_candidateCorrectionNewColumns.data()));
    }
    string str_correctionType = to_string(cmd_info_->correctionType_);
    if(cmd_info_->correctionType_ != 0) {
        paras.push_back("--correctionType");
        paras.push_back((char*)(str_correctionType.data()));
    }
    string str_correctionTypeCands = to_string(cmd_info_->correctionTypeCands_);
    if(cmd_info_->correctionTypeCands_ != 0) {
        paras.push_back("--correctionTypeCands");
        paras.push_back((char*)(str_correctionTypeCands.data()));
    }






    for(int i = 0; i < paras.size(); i++) {
        printf("%s ", paras[i]);
    }
    printf("\n");

    printf("start care part...\n");

    printf("now output to queue, %p %p %p\n", careQueue1, careQueue2, &producerDone);
    main_correction(paras.size(), &(paras[0]), careQueue1, careQueue2, &producerDone, &careStartWrite, &changeNum);

    printf("care end\n");
    printf("care1 queue size %d\n", careQueue1->size_approx());
    printf("care2 queue size %d\n", careQueue2->size_approx());
    printf("care change size %d\n", changeNum);

    careDone1 = 1;
    careDone2 = 1;

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
    } else if(cmd_info_->do_correction_with_care_){
        pair<char *, int> last1;
        pair<char *, int> last2;
        last1.first = new char[1 << 20];
        last1.second = 0;
        last2.first = new char[1 << 20];
        last2.second = 0;
        while (true) {
            rabbit::fq::FastqDataPairChunk *fqdatachunk;
            fqdatachunk = fqFileReader->readNextPairChunkParallel(careQueue1, careQueue2, &careDone1, &careDone2, last1,
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
    cout << "file " << file << " has " << n_chunks << " chunks" << endl;
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
        vector<neoReference> data1, data2;
        vector<neoReference> pass_data1, pass_data2;
        rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk *) (fqdatachunk->left_part), data1, true);
        rabbit::fq::chunkFormat((rabbit::fq::FastqDataChunk *) (fqdatachunk->right_part), data2, true);
        ASSERT(data1.size() == data2.size());
        int out_len1 = 0, out_len2 = 0;
        int b_size = min(data1.size(), data2.size());
        for (int i = 0; i < b_size; i++) {
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
                Adapter::CorrectData(item1, item2, overlap_res, cmd_info_->isPhred64_);
            }
            if (trim_res1 && trim_res2 && cmd_info_->trim_adapter_) {
                int trimmed= false;
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
                int res1 = 0, res2 = 0;
                bool is_trimmed1 = trimmed;
                bool is_trimmed2 = trimmed;

                if (!trimmed) {
                    if (cmd_info_->detect_adapter1_) {
                        int res1;
                        if (cmd_info_->print_what_trimmed_) {
                            res1 = Adapter::TrimAdapter(item1, cmd_info_->adapter_seq1_,
                                    thread_info->aft_state1_->adapter_map_,
                                    cmd_info_->adapter_len_lim_, false);
                        } else {
                            res1 = Adapter::TrimAdapter(item1, cmd_info_->adapter_seq1_, false);
                        }
                        if (res1) {
                            is_trimmed1 = true;
                            thread_info->aft_state1_->AddTrimAdapter();
                            thread_info->aft_state1_->AddTrimAdapterBase(res1);
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
                        if (res2) {
                            is_trimmed2 = true;
                            thread_info->aft_state1_->AddTrimAdapter();
                            thread_info->aft_state1_->AddTrimAdapterBase(res2);
                        }

                    }


                }

                if(cmd_info_->adapter_from_fasta_.size() > 0) {
                    if (cmd_info_->print_what_trimmed_) {
                        res1 = Adapter::TrimAdapters(item1, cmd_info_->adapter_from_fasta_,
                                thread_info->aft_state1_->adapter_map_, cmd_info_->adapter_len_lim_, false);
                    } else {
                        res1 = Adapter::TrimAdapters(item1, cmd_info_->adapter_from_fasta_, false);
                    }
                    if (res1) {
                        if(!is_trimmed1) thread_info->aft_state1_->AddTrimAdapter();
                        thread_info->aft_state1_->AddTrimAdapterBase(res1);
                    }

                    if (cmd_info_->print_what_trimmed_) {
                        res2 = Adapter::TrimAdapters(item2, cmd_info_->adapter_from_fasta_,
                                thread_info->aft_state2_->adapter_map_, cmd_info_->adapter_len_lim_, false);
                    } else {
                        res2 = Adapter::TrimAdapters(item2, cmd_info_->adapter_from_fasta_, false);
                    }
                    if (res2) {
                        if(!is_trimmed2) thread_info->aft_state1_->AddTrimAdapter();
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

                if(cmd_info_->notKeepOrder_ == 0) {
                    while(nowChunkId != id) {
                        usleep(100);
                    }
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
                nowChunkId++;
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

                    if(cmd_info_->notKeepOrder_ == 0) {
                        while(nowChunkId != id) {
                            usleep(100);
                        }
                    }
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
                    nowChunkId++;
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
    pair<char *, int> now;
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
    pair<char *, int> now;
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
                    usleep(100);
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

/**
 * @brief do QC for pair-end data
 */
void PeQc::ProcessPeFastq() {
    if (cmd_info_->interleaved_in_) {
        auto *fastqPool = new rabbit::fq::FastqDataPool(64, 1 << 22);
        rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> queue1(64, 1);
        auto **p_thread_info = new ThreadInfo *[cmd_info_->thread_number_];
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            p_thread_info[t] = new ThreadInfo(cmd_info_, true);
        }
        thread *write_thread1;
        thread *write_thread2;
        if (cmd_info_->write_data_) {
            write_thread1 = new thread(bind(&PeQc::WriteSeFastqTask1, this));
            if (cmd_info_->interleaved_out_ == 0)
                write_thread2 = new thread(bind(&PeQc::WriteSeFastqTask2, this));
        }

        thread producer(
                bind(&PeQc::ProducerPeInterFastqTask, this, cmd_info_->in_file_name1_, fastqPool,
                    ref(queue1)));
        auto **threads = new thread *[cmd_info_->thread_number_];
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            threads[t] = new thread(
                    bind(&PeQc::ConsumerPeInterFastqTask, this, p_thread_info[t], fastqPool, ref(queue1)));
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
        vector<State *> pre_vec_state1;
        vector<State *> pre_vec_state2;
        vector<State *> aft_vec_state1;
        vector<State *> aft_vec_state2;

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
        thread *carer;
        if(cmd_info_->do_correction_with_care_) {
            carer = new thread(bind(&PeQc::careProcess, this));
            //cmd_info_->in_file_name1_ = "./tmp.fq";
            while(careStartWrite == 0) {
                usleep(10000);
            }
            printf("QC start...\n");
            if(changeNum == 0) {
                carer->join();
                delete carer;
                cmd_info_->do_correction_with_care_ = 0;
            } else {
                cmd_info_->use_pugz_ = 0;
            }
        }
        if(cmd_info_->use_pigz_) {
            if(cmd_info_->do_correction_with_care_) {
                carer->join();
                delete carer;
            }
        }


        thread *pugzer1;
        thread *pugzer2;

        if (cmd_info_->use_pugz_) {
            pugzer1 = new thread(bind(&::PeQc::PugzTask1, this));
            pugzer2 = new thread(bind(&::PeQc::PugzTask2, this));
        }


        auto *fastqPool = new rabbit::fq::FastqDataPool(64, 1 << 22);
        rabbit::core::TDataQueue<rabbit::fq::FastqDataPairChunk> queue1(64, 1);
        auto **p_thread_info = new ThreadInfo *[cmd_info_->thread_number_];
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            p_thread_info[t] = new ThreadInfo(cmd_info_, true);
        }
        thread *write_thread1;
        thread *write_thread2;
        if (cmd_info_->write_data_) {
            write_thread1 = new thread(bind(&PeQc::WriteSeFastqTask1, this));
            if (cmd_info_->interleaved_out_ == 0)
                write_thread2 = new thread(bind(&PeQc::WriteSeFastqTask2, this));
        }
        thread *pigzer1;
        thread *pigzer2;
        if (cmd_info_->use_pigz_) {
            pigzer1 = new thread(bind(&PeQc::PigzTask1, this));
            if (cmd_info_->interleaved_out_ == 0)
                pigzer2 = new thread(bind(&PeQc::PigzTask2, this));
        }
        thread producer(
                bind(&PeQc::ProducerPeFastqTask, this, cmd_info_->in_file_name1_, cmd_info_->in_file_name2_,
                    fastqPool, ref(queue1)));
        auto **threads = new thread *[cmd_info_->thread_number_];
        for (int t = 0; t < cmd_info_->thread_number_; t++) {
            threads[t] = new thread(
                    bind(&PeQc::ConsumerPeFastqTask, this, p_thread_info[t], fastqPool, ref(queue1)));
        }

        if (cmd_info_->use_pugz_) {
            pugzer1->join();
            pugzer2->join();
        }
        if(cmd_info_->use_pigz_ == 0) {
            if(cmd_info_->do_correction_with_care_) {
                carer->join();
                delete carer;
            }
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
        vector<State *> pre_vec_state1;
        vector<State *> pre_vec_state2;
        vector<State *> aft_vec_state1;
        vector<State *> aft_vec_state2;

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
