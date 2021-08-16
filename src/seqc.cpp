//
// Created by ylf9811 on 2021/7/6.
//
#include "seqc.h"

/**
 * @brief Construct function
 * @param cmd_info1 : cmd information
 */
SeQc::SeQc(CmdInfo *cmd_info1) {
    cmd_info_ = cmd_info1;
    filter_ = new Filter(cmd_info1);
    done_thread_number_ = 0;
    int out_block_nums = int(1.0 * cmd_info1->in_file_size1_ / cmd_info1->out_block_size_);
    printf("out_block_nums %d\n", out_block_nums);
    out_queue_ = NULL;
    if (cmd_info1->write_data_) {
        out_queue_ = new moodycamel::ConcurrentQueue<std::pair<char *, int>>
                (out_block_nums + 1000);
        printf("open stream %s\n", cmd_info1->out_file_name1_.c_str());
        out_stream_ = std::fstream(cmd_info1->out_file_name1_, std::ios::out | std::ios::binary);
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


/**
 * @brief get fastq data chunk from fastq_data_pool and put it into data queue
 * @param file : fastq file name, which is also the input file
 * @param fastq_data_pool : fastq data pool
 * @param dq : a data queue
 */

void SeQc::ProducerSeFastqTask(std::string file, rabbit::fq::FastqDataPool *fastq_data_pool,
                               rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> &dq) {
    rabbit::fq::FastqFileReader *fqFileReader;
    fqFileReader = new rabbit::fq::FastqFileReader(file, fastq_data_pool);
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
}
//
//void SeQc::PrintRead(neoReference &ref) {
//    std::cout << std::string((char *) ref.base + ref.pname, ref.lname) << std::endl;
//    std::cout << std::string((char *) ref.base + ref.pseq, ref.lseq) << std::endl;
//    std::cout << std::string((char *) ref.base + ref.pstrand, ref.lstrand) << std::endl;
//    std::cout << std::string((char *) ref.base + ref.pqual, ref.lqual) << std::endl;
//}

std::string SeQc::Read2String(neoReference &ref) {
    return std::string((char *) ref.base + ref.pname, ref.lname) + "\n" +
           std::string((char *) ref.base + ref.pseq, ref.lseq) + "\n" +
           std::string((char *) ref.base + ref.pstrand, ref.lstrand) + "\n" +
           std::string((char *) ref.base + ref.pqual, ref.lqual) + "\n";
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
            std::vector<neoReference> data;
            rabbit::fq::chunkFormat(fqdatachunk, data, true);
            for (auto item:data) {
                thread_info->TGS_state_->tgsStatRead(item);
            }
            fastq_data_pool->Release(fqdatachunk);
        }
    } else {
        while (dq.Pop(id, fqdatachunk)) {
            std::vector<neoReference> data;
            std::vector<neoReference> pass_data;
            rabbit::fq::chunkFormat(fqdatachunk, data, true);
            int out_len = 0;
            for (auto item:data) {
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


                if (trim_res && cmd_info_->trim_adapter_ && cmd_info_->detect_adapter1_) {
                    Adapter::TrimAdapter(item, cmd_info_->adapter_seq1_, false);
                }
                int filter_res = filter_->ReadFiltering(item, trim_res);
                if (filter_res == 0) {
                    thread_info->aft_state1_->StateInfo(item);
                    if (cmd_info_->write_data_) {
                        pass_data.push_back(item);
                        out_len += item.lname + item.lseq + item.lstrand + item.lqual + 4;
                    }
                }
            }

            if (cmd_info_->write_data_) {
                if (pass_data.size() > 0) {
                    char *out_data = new char[out_len];
                    int pos = 0;
                    for (auto item:pass_data) {
                        //TODO delete name
                        Read2Chars(item, out_data, pos);
                    }
                    ASSERT(pos == out_len);
                    out_queue_->enqueue({out_data, out_len});
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
    int cnt = 0;
    while (true) {
        if (out_queue_->size_approx() == 0 && done_thread_number_ == cmd_info_->thread_number_) {
            break;
        }
        if (out_queue_->size_approx() == 0) {
            usleep(100);
        }
        std::pair<char *, int> now;
        while (out_queue_->size_approx()) {
            out_queue_->try_dequeue(now);
//            printf("write thread working, write %d %d\n", now.second, cnt++);
            out_stream_.write(now.first, now.second);
            delete[] now.first;
        }
    }
    out_stream_.close();
}

/**
 * @brief do QC for single-end data
 */

void SeQc::ProcessSeFastq() {
    auto *fastqPool = new rabbit::fq::FastqDataPool(256, 1 << 22);
    //TODO replace this queue
    rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> queue1(128, 1);

    auto **p_thread_info = new ThreadInfo *[cmd_info_->thread_number_];
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        p_thread_info[t] = new ThreadInfo(cmd_info_);
    }
    std::thread *write_thread;
    if (cmd_info_->write_data_) {
        write_thread = new std::thread(std::bind(&SeQc::WriteSeFastqTask, this));
    }
    //TODO bind ?
    std::thread producer(
            std::bind(&SeQc::ProducerSeFastqTask, this, cmd_info_->in_file_name1_, fastqPool, std::ref(queue1)));
    auto **threads = new std::thread *[cmd_info_->thread_number_];
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        threads[t] = new std::thread(
                std::bind(&SeQc::ConsumerSeFastqTask, this, p_thread_info[t], fastqPool, std::ref(queue1)));
    }
    producer.join();
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        threads[t]->join();
    }
    if (cmd_info_->write_data_) {
        write_thread->join();
    }

    printf("all thrad done\n");
    printf("now merge thread info\n");

    std::vector<State *> pre_vec_state;
    std::vector<State *> aft_vec_state;

    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        pre_vec_state.push_back(p_thread_info[t]->pre_state1_);
        aft_vec_state.push_back(p_thread_info[t]->aft_state1_);
    }
    auto pre_state = State::MergeStates(pre_vec_state);
    auto aft_state = State::MergeStates(aft_vec_state);
//
//    for (auto item:pre_vec_state) {
//        printf("now pre state hash graph size %d\n", item->GetHashNum());
//    }
//    for (auto item:aft_vec_state) {
//        printf("now aft state hash graph size %d\n", item->GetHashNum());
//    }

    printf("merge done\n");
    printf("print pre state info :\n");
    State::PrintStates(pre_state);
    printf("print aft state info :\n");
    State::PrintStates(aft_state);


//    auto OverRepSeq1 = pre_state->GetHotSeqsInfo();
    auto hash_graph1 = pre_state->GetHashGraph();
    int hash_num1 = pre_state->GetHashNum();
//    cout << "=============OverRepSeq1=============" << endl;
//    for (auto it:OverRepSeq1) {
//        cout << it.first << " " << it.second << endl;
//    }
//    cout << "=====================================" << endl;

//    auto OverRepSeq2 = aft_state->GetHotSeqsInfo();
    auto hash_graph2 = aft_state->GetHashGraph();
    int hash_num2 = aft_state->GetHashNum();
//    cout << "=============OverRepSeq2=============" << endl;
//    for (auto it:OverRepSeq2) {
//        cout << it.first << " " << it.second << endl;
//    }
//    cout << "=====================================" << endl;


//    ofstream ofs;
//    ofs.open("ORP2.log", ifstream::out);
//    for (auto it:OverRepSeq1) {
//        ofs << it.first << " " << it.second << "\n";
//    }
//    ofs.close();
//    ofs.open("ORP3.log", ifstream::out);
//    for (auto it:OverRepSeq2) {
//        ofs << it.first << " " << it.second << "\n";
//    }
//    ofs.close();
//
    ofstream ofs;
    ofs.open("ORP2.log", ifstream::out);
    for (int i = 0; i < hash_num1; i++) {
        ofs << hash_graph1[i].seq << " " << hash_graph1[i].cnt << "\n";
    }
    ofs.close();
    ofs.open("ORP3.log", ifstream::out);
    for (int i = 0; i < hash_num2; i++) {
        ofs << hash_graph2[i].seq << " " << hash_graph2[i].cnt << "\n";
    }
    ofs.close();


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

    Repoter::ReportHtmlSe(pre_state, aft_state, cmd_info_->in_file_name1_, dupRate * 100.0);


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


void printCSS(ofstream &ofs) {
    ofs << "<style type=\"text/css\">" << endl;
    ofs << "td {border:1px solid #dddddd;padding:5px;font-size:12px;}" << endl;
    ofs << "table {border:1px solid #999999;padding:2x;border-collapse:collapse; width:800px}" << endl;
    ofs << ".col1 {width:240px; font-weight:bold;}" << endl;
    ofs << ".adapter_col {width:500px; font-size:10px;}" << endl;
    ofs << "img {padding:30px;}" << endl;
    ofs << "#menu {font-family:Consolas, 'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs
            << "#menu a {color:#0366d6; font-size:18px;font-weight:600;line-height:28px;text-decoration:none;font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol'}"
            << endl;
    ofs << "a:visited {color: #999999}" << endl;
    ofs << ".alignleft {text-align:left;}" << endl;
    ofs << ".alignright {text-align:right;}" << endl;
    ofs << ".figure {width:800px;height:600px;}" << endl;
    ofs << ".header {color:#ffffff;padding:1px;height:20px;background:#000000;}" << endl;
    ofs
            << ".section_title {color:#ffffff;font-size:20px;padding:5px;text-align:left;background:#663355; margin-top:10px;}"
            << endl;
    ofs << ".subsection_title {font-size:16px;padding:5px;margin-top:10px;text-align:left;color:#663355}" << endl;
    ofs
            << "#container {text-align:center;padding:3px 3px 3px 10px;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}"
            << endl;
    ofs << ".menu_item {text-align:left;padding-top:5px;font-size:18px;}" << endl;
    ofs << ".highlight {text-align:left;padding-top:30px;padding-bottom:30px;font-size:20px;line-height:35px;}" << endl;
    ofs << "#helper {text-align:left;border:1px dotted #fafafa;color:#777777;font-size:12px;}" << endl;
    ofs
            << "#footer {text-align:left;padding:15px;color:#ffffff;font-size:10px;background:#663355;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}"
            << endl;
    ofs << ".kmer_table {text-align:center;font-size:8px;padding:2px;}" << endl;
    ofs << ".kmer_table td{text-align:center;font-size:8px;padding:0px;color:#ffffff}" << endl;
    ofs << ".sub_section_tips {color:#999999;font-size:10px;padding-left:5px;padding-bottom:3px;}" << endl;
    ofs << "</style>" << endl;
}

void printJS(ofstream &ofs) {
    ofs << "<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>" << endl;
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << "    function showOrHide(divname) {" << endl;
    ofs << "        div = document.getElementById(divname);" << endl;
    ofs << "        if(div.style.display == 'none')" << endl;
    ofs << "            div.style.display = 'block';" << endl;
    ofs << "        else" << endl;
    ofs << "            div.style.display = 'none';" << endl;
    ofs << "    }" << endl;
    ofs << "</script>" << endl;
}


void printHeader(ofstream &ofs) {
    ofs << "<html><head><meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />";
    //ofs << "<title>RabbitQC report at " + getCurrentSystemTime() + " </title>";
    ofs << "<title>RabbitQC report at 111 </title>";
    printJS(ofs);
    printCSS(ofs);
    ofs << "</head>";
    ofs << "<body><div id='container'>";
}

const string getCurrentSystemTime() {
    auto tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    struct tm *ptm = localtime(&tt);
    char date[60] = {0};
    sprintf(date, "%d-%02d-%02d      %02d:%02d:%02d",
            (int) ptm->tm_year + 1900, (int) ptm->tm_mon + 1, (int) ptm->tm_mday,
            (int) ptm->tm_hour, (int) ptm->tm_min, (int) ptm->tm_sec);
    return std::string(date);
}

void printFooter(ofstream &ofs) {
    ofs << "\n</div>" << endl;
    ofs << "<div id='footer'> ";
    ofs << "<p>" << "" << "</p>";
    ofs << "RabbitQC " << "" << ", at " << getCurrentSystemTime() << " </div>";
    ofs << "</body></html>";
}

void report3(TGSStats *preStats1) {
    ofstream ofs;
    ofs.open("TGS.html", ifstream::out);

    printHeader(ofs);

    //printSummary(ofs, result, preStats1, postStats1, preStats2, postStats2);

    ofs << "<div class='section_div'>\n";
    ofs
            << "<div class='section_title' onclick=showOrHide('QC_information')><a name='summary'>QC information</a></div>\n";
    ofs << "<div id='QC_information'>\n";

    if (preStats1) {
        preStats1->reportHtml(ofs, "QC information", "read1");
    }

    //if(preStats2) {
    //    preStats2 -> reportHtml(ofs, "Before filtering", "read2");
    //}

    ofs << "</div>\n";
    ofs << "</div>\n";

    //ofs << "<div class='section_div'>\n";
    //ofs << "<div class='section_title' onclick=showOrHide('after_filtering')><a name='summary'>After filtering</a></div>\n";
    //ofs << "<div id='after_filtering'>\n";

    //if(postStats1) {
    //    postStats1 -> reportHtml(ofs, "After filtering", "read1");
    //}

    //if(postStats2) {
    //    postStats2 -> reportHtml(ofs, "After filtering", "read2");
    //}

    //ofs << "</div>\n";
    //ofs << "</div>\n";

    printFooter(ofs);
    ofs.close();

}

void SeQc::ProcessSeTGS() {
    auto *fastqPool = new rabbit::fq::FastqDataPool(256, 1 << 22);
    //TODO replace this queue
    rabbit::core::TDataQueue<rabbit::fq::FastqDataChunk> queue1(128, 1);

    auto **p_thread_info = new ThreadInfo *[cmd_info_->thread_number_];
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        p_thread_info[t] = new ThreadInfo(cmd_info_);
    }
    //TODO bind ?
    std::thread producer(
            std::bind(&SeQc::ProducerSeFastqTask, this, cmd_info_->in_file_name1_, fastqPool, std::ref(queue1)));
    auto **threads = new std::thread *[cmd_info_->thread_number_];
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        threads[t] = new std::thread(
                std::bind(&SeQc::ConsumerSeFastqTask, this, p_thread_info[t], fastqPool, std::ref(queue1)));
    }
    producer.join();
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        threads[t]->join();
    }

    printf("all thrad done\n");
    printf("now merge thread info\n");

    std::vector<TGSStats *> vec_state;

    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        vec_state.push_back(p_thread_info[t]->TGS_state_);
    }
    auto mer_state = TGSStats::merge(vec_state);

    printf("merge done\n");
    printf("print TGS state info :\n");
    mer_state->print();

//    report3(mer_state);

    Repoter::ReportHtmlTGS(mer_state, cmd_info_->in_file_name1_);

    delete fastqPool;
    for (int t = 0; t < cmd_info_->thread_number_; t++) {
        delete threads[t];
        delete p_thread_info[t];
    }

    delete[] threads;
    delete[] p_thread_info;
}