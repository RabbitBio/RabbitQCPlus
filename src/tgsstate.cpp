//
// Created by ylf9811 on 2021/8/11.
//

#include "tgsstate.h"


TGSStats::TGSStats(int minLen) {
    mMinlen = minLen;
    mHalfMinlen = mMinlen >> 1;
    ///mLengths;
    int size_range = mMinlen >> 1;
    head_qual_sum = new int64_t[size_range];
    tail_qual_sum = new int64_t[size_range];
    int i;
    for (i = 0; i < 4; i++) {
        head_seq_pos_count[i] = new int64_t[size_range];
        tail_seq_pos_count[i] = new int64_t[size_range];
        memset(head_seq_pos_count[i], 0, sizeof(int64_t) * size_range);
        memset(tail_seq_pos_count[i], 0, sizeof(int64_t) * size_range);
    }
    //init
    memset(head_qual_sum, 0, size_range * sizeof(int64_t));
    memset(tail_qual_sum, 0, size_range * sizeof(int64_t));


}

TGSStats::~TGSStats() {
    delete[] head_qual_sum;
    delete[] tail_qual_sum;
    int i;
    for (i = 0; i < 4; i++) {
        delete[] head_seq_pos_count[i];
        delete[] tail_seq_pos_count[i];
    }
}

void TGSStats::tgsStatRead(neoReference &ref, bool isPhred64) {
    const int rlen = ref.lseq;
    const char *seq = reinterpret_cast<const char *>(ref.base + ref.pseq);
    const char *quality = reinterpret_cast<const char *>(ref.base + ref.pqual);
    int size_range = mMinlen >> 1;
    int i;
    char c1, c2;
    mTotalReadsLen.push_back(rlen);
    int phredSub = 33;
    if (isPhred64)phredSub = 64;
    if (rlen > mMinlen) {
        //[1] stats lengths
        mLengths.push_back(rlen);
        //--head
        for (i = 0; i < size_range; ++i) {
            c1 = seq[i];
            //[2] quality sum
            head_qual_sum[i] += (quality[i] - phredSub);
            //[3] sequence count
            if (c1 == 'A') {
                head_seq_pos_count[0][i]++;
            } else if (c1 == 'T') {
                head_seq_pos_count[1][i]++;
            } else if (c1 == 'C') {
                head_seq_pos_count[2][i]++;
            } else if (c1 == 'G') {
                head_seq_pos_count[3][i]++;
            }

        }
        //--tail
        for (i = rlen - 1; i >= rlen - size_range; --i) {
            c2 = seq[i];
            tail_qual_sum[i - (rlen - size_range)] += (quality[i] - phredSub);

            if (c2 == 'A') {
                tail_seq_pos_count[0][i - (rlen - size_range)]++;
            } else if (c2 == 'T') {
                tail_seq_pos_count[1][i - (rlen - size_range)]++;
            } else if (c2 == 'C') {
                tail_seq_pos_count[2][i - (rlen - size_range)]++;
            } else if (c2 == 'G') {
                tail_seq_pos_count[3][i - (rlen - size_range)]++;
            }
        }
    }
}

void TGSStats::print() {
    //cerr << "nothing here" << std::endl;
    int i;
    std::cout << "head A    T   C   G" << std::endl;
    for (i = 0; i < (mMinlen >> 1); ++i) {
        std::cout << head_seq_pos_count[0][i] << " ";
        std::cout << head_seq_pos_count[1][i] << " ";
        std::cout << head_seq_pos_count[2][i] << " ";
        std::cout << head_seq_pos_count[3][i] << " ";
        std::cout << std::endl;
    }
    std::cout << "tail A    T   C   G" << std::endl;
    //std::cout << std::endl;
    for (i = 0; i < (mMinlen >> 1); ++i) {
        std::cout << tail_seq_pos_count[0][i] << " ";
        std::cout << tail_seq_pos_count[1][i] << " ";
        std::cout << tail_seq_pos_count[2][i] << " ";
        std::cout << tail_seq_pos_count[3][i] << " ";
        std::cout << std::endl;
    }
    std::cout << "head and tail quality" << std::endl;
    for (i = 0; i < (mMinlen >> 1); ++i) {
        std::cout << head_qual_sum[i] << " " << tail_qual_sum[i];
        std::cout << std::endl;
    }

}

TGSStats *TGSStats::merge(std::vector<TGSStats *> &list) {
    int i;
    const int minLen = list[0]->mMinlen;
    TGSStats *s = new TGSStats(minLen);
    for (TGSStats *ds : list) {
        //mLengths
        //s->mLengths.push_back(ds->mLengths);
        //ds->print();
        s->mLengths.insert(s->mLengths.end(), ds->mLengths.begin(), ds->mLengths.end());
        s->mTotalReadsLen.insert(s->mTotalReadsLen.end(), ds->mTotalReadsLen.begin(), ds->mTotalReadsLen.end());
        for (i = 0; i < (minLen >> 1); ++i) {
            //head_seq
            s->head_seq_pos_count[0][i] += ds->head_seq_pos_count[0][i];
            s->head_seq_pos_count[1][i] += ds->head_seq_pos_count[1][i];
            s->head_seq_pos_count[2][i] += ds->head_seq_pos_count[2][i];
            s->head_seq_pos_count[3][i] += ds->head_seq_pos_count[3][i];
            //tail_seq
            s->tail_seq_pos_count[0][i] += ds->tail_seq_pos_count[0][i];
            s->tail_seq_pos_count[1][i] += ds->tail_seq_pos_count[1][i];
            s->tail_seq_pos_count[2][i] += ds->tail_seq_pos_count[2][i];
            s->tail_seq_pos_count[3][i] += ds->tail_seq_pos_count[3][i];
        }
        for (i = 0; i < (minLen >> 1); ++i) {
            //head_qual
            s->head_qual_sum[i] += ds->head_qual_sum[i];
            //tail_qual
            s->tail_qual_sum[i] += ds->tail_qual_sum[i];
        }
    }
    //s->head_qual_mean = s->head_qual_sum/s->mLengths.size()
    return s;

}


//generate html data

bool TGSStats::isLongRead() {
    return mHalfMinlen > 300;
}

int TGSStats::base2num(std::string base) {
    if (base == "A")
        return 0;
    else if (base == "C")
        return 1;
    else if (base == "G")
        return 2;
    else if (base == "T")
        return 3;

    return -1;//fix warning (no ACGT)
}


std::string TGSStats::list2string(double *list, int size) {
    std::stringstream ss;
    for (int i = 0; i < size; i++) {
        ss << list[i];
        if (i < size - 1)
            ss << ",";
    }
    return ss.str();
}

std::string TGSStats::list2string(double *list, int size, int64_t *coords) {
    std::stringstream ss;
    for (int i = 0; i < size; i++) {
        // coords is 1,2,3,...
        int64_t start = 0;
        if (i > 0)
            start = coords[i - 1];
        int64_t end = coords[i];

        double total = 0.0;
        for (int k = start; k < end; k++)
            total += list[k];

        // get average
        if (end == start)
            ss << "0.0";
        else
            ss << total / (end - start);
        //ss << list[coords[i]-1];
        if (i < size - 1)
            ss << ",";
    }
    return ss.str();
}

std::string TGSStats::list2string(int64_t *list, int size) {
    std::stringstream ss;
    for (int i = 0; i < size; i++) {
        ss << list[i];
        if (i < size - 1)
            ss << ",";
    }
    return ss.str();
}

std::string TGSStats::list2stringReversedOrder(int64_t *list, int size) {
    std::stringstream ss;
    for (int i = size - 1; i >= 0; --i) {
        ss << list[i];
        if (i > 0)
            ss << ",";
    }
    return ss.str();
}


void TGSStats::reportHtml(std::ofstream &ofs, std::string filteringType, std::string seqFileName) {
    int64_t seq_count = mLengths.size();
    double *head_base_pos_persent[4];
    double *tail_base_pos_persent[4];
    double *head_quality_mean = new double[mHalfMinlen];
    double *tail_quality_mean = new double[mHalfMinlen];
    for (int i = 0; i < 4; ++i) {
        head_base_pos_persent[i] = new double[mHalfMinlen];
        tail_base_pos_persent[i] = new double[mHalfMinlen];
        memset(head_base_pos_persent[i], 0, sizeof(double) * mHalfMinlen);
        memset(tail_base_pos_persent[i], 0, sizeof(double) * mHalfMinlen);
    }
    memset(head_quality_mean, 0, sizeof(double) * mHalfMinlen);
    memset(tail_quality_mean, 0, sizeof(double) * mHalfMinlen);

    for (int i = 0; i < mHalfMinlen; ++i) {
        int64_t head_total_base_per_pos = head_seq_pos_count[0][i]
                                          + head_seq_pos_count[1][i]
                                          + head_seq_pos_count[2][i]
                                          + head_seq_pos_count[3][i];
        head_base_pos_persent[0][i] = head_seq_pos_count[0][i] * (1.0 / head_total_base_per_pos);
        head_base_pos_persent[1][i] = head_seq_pos_count[1][i] * (1.0 / head_total_base_per_pos);
        head_base_pos_persent[2][i] = head_seq_pos_count[2][i] * (1.0 / head_total_base_per_pos);
        head_base_pos_persent[3][i] = head_seq_pos_count[3][i] * (1.0 / head_total_base_per_pos);
        int64_t tail_total_base_per_pos = tail_seq_pos_count[0][i]
                                          + tail_seq_pos_count[1][i]
                                          + tail_seq_pos_count[2][i]
                                          + tail_seq_pos_count[3][i];
        tail_base_pos_persent[0][i] = tail_seq_pos_count[0][i] * (1.0 / tail_total_base_per_pos);
        tail_base_pos_persent[1][i] = tail_seq_pos_count[1][i] * (1.0 / tail_total_base_per_pos);
        tail_base_pos_persent[2][i] = tail_seq_pos_count[2][i] * (1.0 / tail_total_base_per_pos);
        tail_base_pos_persent[3][i] = tail_seq_pos_count[3][i] * (1.0 / tail_total_base_per_pos);
        head_quality_mean[i] = head_qual_sum[i] * (1.0 / seq_count); // /total_count
        tail_quality_mean[i] = tail_qual_sum[i] * (1.0 / seq_count); // /total_count
    }
    reportHistogram(ofs);
    reportHtmlContents(ofs, seqFileName, 0, "Position in read from start", "Frequency of nucleotide in read",
                       head_base_pos_persent);
    reportHtmlContents(ofs, seqFileName, 1, "Position in read from end", "Frequency of nucleotide in read",
                       tail_base_pos_persent);
    reportHtmlQuality(ofs, seqFileName, 0, "Position in read from start", "Mean quality score of base calls",
                      head_quality_mean);
    reportHtmlQuality(ofs, seqFileName, 1, "Position in read from end", "Mean quality score of base calls",
                      tail_quality_mean);

    //delete

}

//figure ID is made up of seqFileName and figureTitle, so seqFileName+figureTitle must be different for different figure
void TGSStats::reportHtmlQuality(std::ofstream &ofs, std::string seqFileName, bool isTail, std::string xAxisName,
                                 std::string yAxisName, double *statsData) {

    // quality
    std::string figureTitle = (isTail) ? "tail" : "head";
    figureTitle += "_qual_sum";
    std::string subsection = seqFileName + ": " + figureTitle + ": quality";
    std::string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    std::string title = "";

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName
        << "')>" + subsection + "</a></div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n";
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";
    ofs << "</div>\n";

    std::string alphabets[5] = {"A", "T", "C", "G", "mean"};
    std::string colors[5] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", "rgba(0,0,255,1.0)",
                             "rgba(20,20,20,1.0)"};
    ofs << "\n<script type=\"text/javascript\">" << std::endl;
    std::string json_str = "var data=[";

    int64_t *x = new int64_t[mHalfMinlen];
    int total = 0;
    //remove log-style in int64_t reads to be identical with nanoQC by liumy @2017/7/2
    //if(!isLongRead()) {
    for (int i = 0; i < mHalfMinlen; i++) {
        x[total] = i + 1;
        total++;
    }
    //remove log-style in int64_t reads to be identical with nanoQC by liumy @2017/7/2
    //} else {
    //    const int fullSampling = 40;
    //    for(int i=0; i<fullSampling && i<mHalfMinlen; i++){
    //        x[total] = i+1;
    //        total++;
    //    }
    //    // down sampling if it's too int64_t
    //    if(mHalfMinlen>fullSampling) {
    //        double pos = fullSampling;
    //        while(true){
    //            pos *= 1.05;
    //            if(pos >= mHalfMinlen)
    //                break;
    //            x[total] = (int)pos;
    //            total++;
    //        }
    //        // make sure lsat one is contained
    //        if(x[total-1] != mHalfMinlen){
    //            x[total] = mHalfMinlen;
    //            total++;
    //        }
    //    }
    //}
    // four bases
    //for (int b = 0; b<5; b++) {
    std::string base = "TEST"; //alphabets[b];
    json_str += "{";
    if (isTail)
        json_str += "x:[" + list2stringReversedOrder(x, total) + "],";
    else
        json_str += "x:[" + list2string(x, total) + "],";
    //json_str += "y:[" + list2string(mQualityCurves[base], total, x) + "],";
    //json_str += "y:[" + list2string(head_qual_sum, total, x) + "],";
    json_str += "y:[" + list2string(statsData, total, x) + "],";
    json_str += "name: '" + base + "',";
    json_str += "mode:'lines',";
    json_str += "line:{color:'" + colors[0] + "', width:1}\n";
    json_str += "},";
    //}
    json_str += "];\n";
    json_str += "var layout={title:'" + title + "', xaxis:{title:'" + xAxisName + "'";
    if (isTail)
        json_str += ",autorange:'reversed'";
    // use log plot if it's too int64_t
    //remove log-style in int64_t reads to be identical with nanoQC by liumy @2017/7/2
    //if(isLongRead()) {
    //    json_str += ",type:'log'";
    //}
    json_str += "}, yaxis:{title:'" + yAxisName + "'}};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << std::endl;

    delete[] x;
}


//figure ID is made up of seqFileName and figureTitle, so seqFileName+figureTitle must be different for different figure
void TGSStats::reportHtmlContents(std::ofstream &ofs, std::string seqFileName, bool isTail, std::string xAxisName,
                                  std::string yAxisName, double **statsData) {

    // content
    std::string figureTitle = (isTail) ? "tail" : "head";
    figureTitle += "_seq_pos_count";
    std::string subsection = seqFileName + ": " + figureTitle + ": base contents";
    std::string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    std::string title = "";

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName
        << "')>" + subsection + "</a></div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n";
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";
    ofs << "</div>\n";

    std::string alphabets[6] = {"A", "T", "C", "G", "N", "GC"};
    std::string colors[6] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", "rgba(0,0,255,1.0)",
                             "rgba(255, 0, 0, 1.0)", "rgba(20,20,20,1.0)"};
    ofs << "\n<script type=\"text/javascript\">" << std::endl;
    std::string json_str = "var data=[";

    int64_t *x = new int64_t[mHalfMinlen];
    int total = 0;
    //remove log-style in int64_t reads to be identical with nanoQC by liumy @2017/7/2
    //if(!isLongRead()) {
    for (int i = 0; i < mHalfMinlen; i++) {
        x[total] = i + 1;
        total++;
    }
    //remove log-style in int64_t reads to be identical with nanoQC by liumy @2017/7/2
    //} else {
    //    const int fullSampling = 40;
    //    for(int i=0; i<fullSampling && i<mHalfMinlen; i++){
    //        x[total] = i+1;
    //        total++;
    //    }
    //    // down sampling if it's too int64_t
    //    if(mHalfMinlen>fullSampling) {
    //        double pos = fullSampling;
    //        while(true){
    //            pos *= 1.05;
    //            if(pos >= mHalfMinlen)
    //                break;
    //            x[total] = (int)pos;
    //            total++;
    //        }
    //        // make sure lsat one is contained
    //        if(x[total-1] != mHalfMinlen){
    //            x[total] = mHalfMinlen;
    //            total++;
    //        }
    //    }
    //}
    // four bases
    //for (int b = 0; b<6; b++) {
    for (int b = 0; b < 4; b++) {
        std::string base = alphabets[b];
        //int64_t count = 0; //total base of A/C/G/T/N/CG
        //if(base.size()==1) {
        //    //need to modify
        //    //char b = base[0] & 0x07;
        //	//count = mBaseContents[b];
        //	count = head_seq_pos_count[][];

        //} else { // for CG
        //    count = mBaseContents['G' & 0x07] + mBaseContents['C' & 0x07] ;
        //}
        //std::string percentage = to_string((double)count * 100.0 / mBases);
        //if(percentage.length()>5)
        //    percentage = percentage.substr(0,5);
        //std::string name = base + "(" + percentage + "%)";
        std::string name = base; //+ "(" + percentage + "%)";

        json_str += "{";
        if (isTail)
            json_str += "x:[" + list2stringReversedOrder(x, total) + "],";
        else
            json_str += "x:[" + list2string(x, total) + "],";
        //json_str += "y:[" + list2string(mContentCurves[base], total, x) + "],";
        //json_str += "y:[" + list2string(head_seq_pos_count[base2num(base)], total, x) + "],";
        json_str += "y:[" + list2string(statsData[base2num(base)], total, x) + "],";
        json_str += "name: '" + name + "',";
        json_str += "mode:'lines',";
        json_str += "line:{color:'" + colors[b] + "', width:1}\n";
        json_str += "},";
    }
    json_str += "];\n";
    json_str += "var layout={title:'" + title + "', xaxis:{title:'" + xAxisName + "'";
    if (isTail)
        json_str += ",autorange:'reversed'";
    // use log plot if it's too int64_t
    //remove log-style in int64_t reads to be identical with nanoQC by liumy @2017/7/2
    //if(isLongRead()) {
    //    json_str += ",type:'log'";
    //}
    json_str += "}, yaxis:{title:'" + yAxisName + "'}};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << std::endl;

    delete[] x;
}

void TGSStats::reportHistogram(std::ofstream &ofs) {
    ofs << "<div id='Length_figure'>\n";
    ofs << "<div class='figure' id='plot_length' style='height:400px;'></div>\n";
    ofs << "</div>\n";

    ofs << "\n<script type=\"text/javascript\">" << std::endl;
    std::string json_str = "var data=[";

    //
    //int64_t total = Stats1->mLengths.size();
    //vector<int> x = Stats1->mLengths;
    int64_t total = mTotalReadsLen.size();
    std::vector<int> x = mTotalReadsLen;

    std::unordered_map<int, int64_t> count;
    //statistics length;
    for (int xi : x) {
        count[xi]++;
    }
    //int* len = new int[count.size()];
    int64_t len[count.size()];
    //int64_t* percents = new int64_t[count.size()];
    int64_t percents[count.size()];
    //memset(len, 0, sizeof(int)*count.size());
    //memset(percents, 0, sizeof(int64_t)*count.size());
    int i = 0;
    //for(auto xi = count.begin(); xi != count.end(); ++xi) {
    int maxLen = 0;
    for (auto xi : count) {
        len[i] = xi.first;
        percents[i] = xi.second;
        maxLen = xi.first > maxLen ? xi.first : maxLen;
        i++;
    }
    //==============end
    json_str += "{";
    json_str += "x:[" + list2string(len, count.size()) + "],";
    json_str += "y:[" + list2string(percents, count.size()) + "],";
    json_str += "name: 'Read length',";
    json_str += "type:'histogram', histfunc:'sum', autobinx:false, xbins:{size:" + std::to_string(maxLen / 1000) +
                " },"; // histogram setting
    json_str += "line:{color:'rgba(128,0,128,1.0)', width:1}\n";
    json_str += "}";

    json_str += "];\n";

    json_str += "var layout={title:'Reads Length Distribution', xaxis:{title:'Length of reads'}, yaxis:{title:'Number of reads'}};\n";
    json_str += "Plotly.newPlot('plot_length', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << std::endl;

    //delete[] len;
    //delete[] percents;
}

int64_t *const *TGSStats::GetHeadSeqPosCount() const {
    return head_seq_pos_count;
}

int64_t *const *TGSStats::GetTailSeqPosCount() const {
    return tail_seq_pos_count;
}

int64_t *TGSStats::GetHeadQualSum() const {
    return head_qual_sum;
}

int64_t *TGSStats::GetTailQualSum() const {
    return tail_qual_sum;
}


