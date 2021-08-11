//
// Created by ylf9811 on 2021/7/12.
//

#include "repoter.h"

void Repoter::PrintRef(neoReference &ref) {
    printf("%s\n", std::string((char *) ref.base + ref.pname, ref.lname).c_str());
    printf("%s\n", std::string((char *) ref.base + ref.pseq, ref.lseq).c_str());
    printf("%s\n", std::string((char *) ref.base + ref.pstrand, ref.lstrand).c_str());
    printf("%s\n", std::string((char *) ref.base + ref.pqual, ref.lqual).c_str());
}

void Repoter::PrintRef(Reference &ref) {
    printf("%s\n", ref.name.c_str());
    printf("%s\n", ref.seq.c_str());
    printf("%s\n", ref.strand.c_str());
    printf("%s\n", ref.quality.c_str());
}

Reference Repoter::GetRevRef(neoReference &ref) {
    Reference res;
    res.name = std::string((char *) ref.base + ref.pname, ref.lname);
    res.seq = std::string((char *) ref.base + ref.pseq, ref.lseq);
    res.strand = std::string((char *) ref.base + ref.pstrand, ref.lstrand);
    res.quality = std::string((char *) ref.base + ref.pqual, ref.lqual);
    reverse(res.quality.begin(), res.quality.end());
    std::string tmp = std::string(ref.lseq, 'A');
    for (int i = 0; i < ref.lseq; i++) {
        tmp[i] = rabbit::reMap[res.seq[ref.lseq - i - 1]];
    }
    res.seq = tmp;
    return res;
}


std::string HTMLHeader() {
    return std::string("<!DOCTYPE html>\n"
                       "<html lang=\"en\">\n"
                       "<head>\n"
                       "    <meta charset=\"UTF-8\">\n"
                       "    <title>Title</title>\n"
                       "    <script src=\"https://cdn.jsdelivr.net/npm/echarts@5.0.2/dist/echarts.min.js\"></script>\n"
                       "</head>\n"
                       "<body>\n");
}

std::string HTMLCss() {
    return "<style type=\"text/css\">\n"
           " @media screen {\n"
           "  div.summary {\n"
           "    width: 18em;\n"
           "    position:fixed;\n"
           "    top: 3em;\n"
           "    margin:1em 0 0 1em;\n"
           "  }\n"
           "  \n"
           "  div.main {\n"
           "    display:block;\n"
           "    position:absolute;\n"
           "    overflow:auto;\n"
           "    height:auto;\n"
           "    width:auto;\n"
           "    top:4.5em;\n"
           "    bottom:2.3em;\n"
           "    left:18em;\n"
           "    right:0;\n"
           "    border-left: 1px solid #CCC;\n"
           "    padding:0 0 0 1em;\n"
           "    background-color: white;\n"
           "    z-index:1;\n"
           "  }\n"
           "  \n"
           "  div.header {\n"
           "    background-color: #EEE;\n"
           "    border:0;\n"
           "    margin:0;\n"
           "    padding: 0.5em;\n"
           "    font-size: 200%;\n"
           "    font-weight: bold;\n"
           "    position:fixed;\n"
           "    width:100%;\n"
           "    top:0;\n"
           "    left:0;\n"
           "    z-index:2;\n"
           "  }\n"
           "\n"
           "  div.footer {\n"
           "    background-color: #EEE;\n"
           "    border:0;\n"
           "    margin:0;\n"
           "\tpadding:0.5em;\n"
           "    height: 1.3em;\n"
           "\toverflow:hidden;\n"
           "    font-size: 100%;\n"
           "    font-weight: bold;\n"
           "    position:fixed;\n"
           "    bottom:0;\n"
           "    width:100%;\n"
           "    z-index:2;\n"
           "  }\n"
           "  \n"
           "  img.indented {\n"
           "    margin-left: 3em;\n"
           "  }\n"
           " }\n"
           " \n"
           " @media print {\n"
           "\timg {\n"
           "\t\tmax-width:100% !important;\n"
           "\t\tpage-break-inside: avoid;\n"
           "\t}\n"
           "\th2, h3 {\n"
           "\t\tpage-break-after: avoid;\n"
           "\t}\n"
           "\tdiv.header {\n"
           "      background-color: #FFF;\n"
           "    }\n"
           "\t\n"
           " }\n"
           " \n"
           " body {    \n"
           "  font-family: sans-serif;   \n"
           "  color: #000;   \n"
           "  background-color: #FFF;\n"
           "  border: 0;\n"
           "  margin: 0;\n"
           "  padding: 0;\n"
           "  }\n"
           "  \n"
           "  div.header {\n"
           "  border:0;\n"
           "  margin:0;\n"
           "  padding: 0.5em;\n"
           "  font-size: 200%;\n"
           "  font-weight: bold;\n"
           "  width:100%;\n"
           "  }    \n"
           "  \n"
           "  #header_title {\n"
           "  display:inline-block;\n"
           "  float:left;\n"
           "  clear:left;\n"
           "  }\n"
           "  #header_filename {\n"
           "  display:inline-block;\n"
           "  float:right;\n"
           "  clear:right;\n"
           "  font-size: 50%;\n"
           "  margin-right:2em;\n"
           "  text-align: right;\n"
           "  }\n"
           "\n"
           "  div.header h3 {\n"
           "  font-size: 50%;\n"
           "  margin-bottom: 0;\n"
           "  }\n"
           "  \n"
           "  div.summary ul {\n"
           "  padding-left:0;\n"
           "  list-style-type:none;\n"
           "  }\n"
           "  \n"
           "  div.summary ul li img {\n"
           "  margin-bottom:-0.5em;\n"
           "  margin-top:0.5em;\n"
           "  }\n"
           "\t  \n"
           "  div.main {\n"
           "  background-color: white;\n"
           "  }\n"
           "      \n"
           "  div.module {\n"
           "  padding-bottom:1.5em;\n"
           "  padding-top:1.5em;\n"
           "  }\n"
           "\t  \n"
           "  div.footer {\n"
           "  background-color: #EEE;\n"
           "  border:0;\n"
           "  margin:0;\n"
           "  padding: 0.5em;\n"
           "  font-size: 100%;\n"
           "  font-weight: bold;\n"
           "  width:100%;\n"
           "  }\n"
           "\n"
           "\n"
           "  a {\n"
           "  color: #000080;\n"
           "  }\n"
           "\n"
           "  a:hover {\n"
           "  color: #800000;\n"
           "  }\n"
           "      \n"
           "  h2 {\n"
           "  color: #800000;\n"
           "  padding-bottom: 0;\n"
           "  margin-bottom: 0;\n"
           "  clear:left;\n"
           "  }\n"
           "\n"
           "  table { \n"
           "  margin-left: 3em;\n"
           "  text-align: center;\n"
           "  }\n"
           "  \n"
           "  th { \n"
           "  text-align: center;\n"
           "  background-color: #000080;\n"
           "  color: #FFF;\n"
           "  padding: 0.4em;\n"
           "  }      \n"
           "  \n"
           "  td { \n"
           "  font-family: monospace; \n"
           "  text-align: left;\n"
           "  background-color: #EEE;\n"
           "  color: #000;\n"
           "  padding: 0.4em;\n"
           "  }\n"
           "\n"
           "  img {\n"
           "  padding-top: 0;\n"
           "  margin-top: 0;\n"
           "  border-top: 0;\n"
           "  }\n"
           "\n"
           "  \n"
           "  p {\n"
           "  padding-top: 0;\n"
           "  margin-top: 0;\n"
           "  }\n"
           "</style>\n";
}

std::string insertDiv(std::string id) {
    return "<div id=\"" + id + "\" style=\"width: 800px;height:600px;\"></div>\n";
}

std::string insertTooltip() {
    return "tooltip: {\n"
           "        trigger: 'axis',\n"
           "        axisPointer: {\n"
           "            type: 'cross',\n"
           "            label: {\n"
           "                backgroundColor: '#6a7985'\n"
           "            }\n"
           "        }\n"
           "    },\n";
}

std::string insertChart(std::string id) {
    return "\nvar " + id + "Chart = echarts.init(document.getElementById(\'" + id + "\'));\n";
}

std::string insertChartOption(std::string id) {
    return id + "Chart.setOption(" + id + "Option);\n";
}

std::string insertOptionBegin(std::string id) {
    return "\nvar " + id + "Option = {\n";
}

std::string insertOptionEnd() {
    return "}\n";
}

std::string insertTitle(std::string text) {
    return "title: {\n"
           "            text: \'" + text + "\',\n"
                                           "            left:\'center\',\n"
                                           "        },\n";
}

std::string insertxAxis(int len, int interval = 1) {
    std::string out("xAxis: {\n"
                    "            type: 'category',\n"
                    "            data: [");
    for (int i = 0; i < len; i += interval) {
        out.append(std::to_string(i) + ',');
    }
    out.append("]\n},\n");
    return out;
}

std::string insertxAxis(std::string name, int len, int interval = 1) {
    std::string out("xAxis: {\n"
                    "            type: 'category',\n"
                    "            name: \'" + name + "\',\n"
                                                    "            nameLocation:'center',\n"
                                                    "nameTextStyle: {\n"
                                                    "                lineHeight: 50,\n"
                                                    "                fontSize: 13,\n"
                                                    "                fontFamily: \"monospace\",\n"
                                                    "                fontWeight: \"bold\"\n"
                                                    "            },\n"
                                                    "            data: [");
    for (int i = 0; i < len; i += interval) {
        out.append(std::to_string(i) + ',');
    }
    out.append("]\n},\n");
    return out;
}

std::string insertyAxis(std::string type, std::string _min, std::string _max) {
    return "yAxis: {\n"
           "            type: \'" + type + "\',\n"
                                           "            min:" + _min + ",\n"
                                                                       "            max:" + _max + "\n"
                                                                                                   "        },\n";
}

std::string insertyAxis(std::string type) {
    return "yAxis: {\n"
           "            type: \'" + type + "\',\n"
                                           "            axisLabel: {\n"
                                           "                formatter: function (value) {\n"
                                           "                    if (value < 100000) return value;\n"
                                           "                    var len=0;\n"
                                           "                    var last=value;\n"
                                           "                    while (value>9){\n"
                                           "                        value=value/10;\n"
                                           "                        len++;\n"
                                           "                    }\n"
                                           "                    return last/Math.pow(10,len)+'E+'+len;\n"
                                           "              },\n"
                                           "           },\n"
                                           "        },\n";
}

std::string insertSeriesBegin() {
    return "series: [\n";
}

std::string insertSeriesEnd() {
    return "],\n";
}

std::string insertSeriesData(std::string type, int *data, int len, int interval = 1) {
    std::string out("{\n"
                    "            type: \'" + type + "\',\n"
                                                    "            data: [\n");
    for (int i = 0; i < len; i += interval) {
        out.append(std::to_string(data[i]) + ',');
    }
    out.append("\n]},\n");
    return out;
}

std::string insertSeriesData(std::string type, std::string name, int *data, int len, int interval = 1) {
    std::string out("{\n"
                    "            type: \'" + type + "\',\n"
                                                    "            name: \'" + name + "\',\n"
                                                                                    "            data: [\n             ");
    for (int i = 0; i < len; i += interval) {
        out.append(std::to_string(data[i]) + ',');
    }
    out.append("\n]},\n");
    return out;
}

std::string insertSeriesData(std::string type, double *data, int len, int interval = 1) {
    std::string out("{\n"
                    "            type: \'" + type + "\',\n"
                                                    "            data: [\n            ");
    for (int i = 0; i < len; i += interval) {
        out.append(std::to_string(data[i]) + ',');
    }
    out.append("\n]},\n");
    return out;
}

std::string insertSeriesData(std::string type, std::string name, double *data, int len, int interval = 1) {
    std::string out("{\n"
                    "            type: \'" + type + "\',\n"
                                                    "            name: \'" + name + "\',\n"
                                                                                    "            data: [\n");
    for (int i = 0; i < len; i += interval) {
        out.append(std::to_string(data[i]) + ',');
    }
    out.append("\n]},\n");
    return out;
}

std::string insertSeriesMultiDataBegin(std::string type, std::string name) {
    return "{\n"
           "            type: \'" + type + "\',\n"
                                           "            name:\'" + name + "\'"
                                                                          "            data: [\n";
}

std::string insertSeriesMultiDataBegin(std::string type) {
    return "{\n"
           "            type: \'" + type + "\',\n"
                                           "            data: [\n";
}

std::string insertSeriesMultiDataEnd() {
    return "],\nitemStyle: {\n"
           "        color: \"#d7ab82\"\n"
           "    },\n},\n";
}

std::string insertSeriesOneData(int *data, int len) {
    std::string out("[");
    for (int i = 0; i < len; i++) {
        out.append(std::to_string(data[i]) + ',');
    }
    out.append("],\n");
    return out;
}

std::string insertDataZoom() {
    return "dataZoom: [\n"
           "        {\n"
           "            id: 'dataZoomX',\n"
           "            type: 'slider',\n"
           "            xAxisIndex: [0],\n"
           "            filterMode: 'filter', // 设定为 'filter' 从而 X 的窗口变化会影响 Y 的范围。\n"
           "            start: 0,\n"
           "            end: 100\n"
           "        },\n"
           "        {\n"
           "            id: 'dataZoomY',\n"
           "            type: 'slider',\n"
           "            yAxisIndex: [0],\n"
           "            filterMode: 'empty',\n"
           "            start: 0,\n"
           "            end: 100\n"
           "        }\n"
           "    ],\n";
}

std::string insertLegend(std::string data) {
    return "legend: {\n"
           "        data: [" + data + "],\n"
                                      "        left:\'85%\',\n"
                                      "    },\n";
}

std::string insertTableBegin() {
    return "<table>\n";
}

std::string insertTableEnd() {
    return "</table>\n";
}

std::string insertTableTitle(std::string str1, std::string str2) {
    return "<thead>\n"
           "    <tr>\n"
           "        <th>" + str1 + "</th>\n"
                                   "        <th>" + str2 + "</th>\n"
                                                           "    </tr>\n"
                                                           "</thead>\n";
}

std::string insertTableTbobyBegin() {
    return "<tboby>\n";
}

std::string insertTableTbobyEnd() {
    return "</tboby>\n";
}

std::string insertTableTr(std::string str1, std::string str2) {
    return "    <tr>\n"
           "        <td>" + str1 + "</td>\n"
                                   "        <td>" + str2 + "</td>\n"
                                                           "    </tr>\n";
}

void Repoter::ReportHtmlSe(State *state1, State *state2, std::string file_name, double dup) {


    printf("tot bases %lld\n", state1->GetTotBases());
    printf("gc bases %lld\n", state1->GetGcBases());
    printf("tot bases %lld\n", state2->GetTotBases());
    printf("gc bases %lld\n", state2->GetGcBases());

    printf("report html se data\n");

    std::string outhtml;

    int mx_len = std::max(state1->GetRealSeqLen(), state2->GetRealSeqLen());


    int *tmp_int = new int[mx_len];
    double *tmp_double = new double[mx_len];


    outhtml.append(HTMLHeader());
    outhtml.append(HTMLCss());
    //Basic Status
    outhtml.append(insertTableBegin());
    outhtml.append(insertTableTitle("General", ""));
    outhtml.append(insertTableTbobyBegin());
    outhtml.append(insertTableTr("RabbitQC version:", "0.0.2"));
    outhtml.append(insertTableTr("FileName", file_name));
    outhtml.append(
            insertTableTr("duplication rate:", std::to_string(dup) + "(may be overestimated since this is SE data)"));
    outhtml.append(insertTableTbobyEnd());
    outhtml.append(insertTableEnd());

    //Before filtering information
    outhtml.append(insertTableBegin());
    outhtml.append(insertTableTitle("Before filtering", ""));
    outhtml.append(insertTableTbobyBegin());
    outhtml.append(insertTableTr("total reads:", std::to_string(state1->GetLines())));
    outhtml.append(insertTableTr("total bases:", std::to_string(state1->GetTotBases())));
    outhtml.append(insertTableTr("Q20 bases:", std::to_string(state1->GetQ20Bases())));
    outhtml.append(insertTableTr("Q30 bases:", std::to_string(state1->GetQ30Bases())));
    outhtml.append(insertTableTr("%GC", std::to_string(1.0 * state1->GetGcBases() / state1->GetTotBases())));
    outhtml.append(insertTableTbobyEnd());
    outhtml.append(insertTableEnd());

    //After filtering information
    outhtml.append(insertTableBegin());
    outhtml.append(insertTableTitle("After filtering", ""));
    outhtml.append(insertTableTbobyBegin());

    outhtml.append(insertTableTr("total reads:", std::to_string(state2->GetLines())));
    outhtml.append(insertTableTr("total bases:", std::to_string(state2->GetTotBases())));
    outhtml.append(insertTableTr("Q20 bases:", std::to_string(state2->GetQ20Bases())));
    outhtml.append(insertTableTr("Q30 bases:", std::to_string(state2->GetQ30Bases())));
    outhtml.append(insertTableTr("%GC", std::to_string(1.0 * state2->GetGcBases() / state2->GetTotBases())));
    outhtml.append(insertTableTbobyEnd());
    outhtml.append(insertTableEnd());


    std::string PositionQuality1("PositionQuality1");
    std::string PositionQuality2("PositionQuality2");
    std::string PositionContent1("PositionContent1");
    std::string PositionContent2("PositionContent2");
    std::string MeanQuality1("MeanQuality1");
    std::string MeanQuality2("MeanQuality2");
    std::string GCContent1("GCContent1");
    std::string GCContent2("GCContent2");


    outhtml.append(insertDiv(PositionQuality1));
    outhtml.append(insertDiv(PositionQuality2));
    outhtml.append(insertDiv(PositionContent1));
    outhtml.append(insertDiv(PositionContent2));


    outhtml.append(insertDiv(MeanQuality1));
    outhtml.append(insertDiv(MeanQuality2));

    outhtml.append(insertDiv(GCContent1));
    outhtml.append(insertDiv(GCContent2));


    outhtml.append("</body>\n");

    //js
    outhtml.append("<script type=\"text/javascript\">\n");

    // Quality Scores cross all bases

    outhtml.append(insertChart(PositionQuality1));
    //option
    outhtml.append(insertOptionBegin(PositionQuality1));
    outhtml.append(insertTitle("Quality Scores cross all bases (Before filtering)"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("Position in read(bp)", mx_len + 1));
    outhtml.append(insertyAxis("value", std::to_string(0), std::to_string(41)));
    outhtml.append(insertSeriesBegin());

    int64_t *pos_qul_ = state1->GetPosQul();
    int64_t *pos_cnt_ = state1->GetPosCnt();
    for (int i = 0; i < mx_len; i++) {
        int64_t sum_qul = 0;
        int64_t sum_cnt = 0;
        for (int j = 0; j < 8; j++) {
            sum_qul += pos_qul_[i * 8 + j];
            sum_cnt += pos_cnt_[i * 8 + j];
        }
        tmp_double[i] = 1.0 * sum_qul / sum_cnt;
    }
    outhtml.append(insertSeriesData("line", tmp_double, mx_len));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(PositionQuality1));


    outhtml.append(insertChart(PositionQuality2));
    //option
    outhtml.append(insertOptionBegin(PositionQuality2));
    outhtml.append(insertTitle("Quality Scores cross all bases (After filtering)"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("Position in read(bp)", mx_len + 1));
    outhtml.append(insertyAxis("value", std::to_string(0), std::to_string(41)));
    outhtml.append(insertSeriesBegin());

    pos_qul_ = state2->GetPosQul();
    pos_cnt_ = state2->GetPosCnt();
    for (int i = 0; i < mx_len; i++) {
        int64_t sum_qul = 0;
        int64_t sum_cnt = 0;
        for (int j = 0; j < 8; j++) {
            sum_qul += pos_qul_[i * 8 + j];
            sum_cnt += pos_cnt_[i * 8 + j];
        }
        tmp_double[i] = 1.0 * sum_qul / sum_cnt;
    }
    outhtml.append(insertSeriesData("line", tmp_double, mx_len));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(PositionQuality2));






//
//    //lengthlist
//    outhtml.append(insertChart(LengthList));
//    //option
//    outhtml.append(insertOptionBegin(LengthList));
//    outhtml.append(insertTitle("Seqence Length List"));
//    outhtml.append(insertTooltip());
//    outhtml.append(insertDataZoom());
//    outhtml.append(insertxAxis("Sequence Length(bp)", max_len + 1));
//    outhtml.append(insertyAxis("value"));
//    outhtml.append(insertSeriesBegin());
//    outhtml.append(insertSeriesData("line", LengthSequence, max_len + 1));
//    outhtml.append(insertSeriesEnd());
//    outhtml.append(insertOptionEnd());
//    outhtml.append(insertChartOption(LengthList));
//

    // Mean Quanlity
    outhtml.append(insertChart(MeanQuality1));
    //option
    outhtml.append(insertOptionBegin(MeanQuality1));
    outhtml.append(insertTitle("Mean Quanlity List(Before filtering)"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("Mean Sequence Quality(Phred Score)", 42));
    outhtml.append(insertyAxis("value"));
    outhtml.append(insertSeriesBegin());
    int64_t *qul_cnt = state1->GetQulCnt();
    for (int i = 0; i < 42; i++)
        tmp_double[i] = qul_cnt[i];
    outhtml.append(insertSeriesData("line", tmp_double, 42));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(MeanQuality1));


    outhtml.append(insertChart(MeanQuality2));
    //option
    outhtml.append(insertOptionBegin(MeanQuality2));
    outhtml.append(insertTitle("Mean Quanlity List(After filtering)"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("Mean Sequence Quality(Phred Score)", 42));
    outhtml.append(insertyAxis("value"));
    outhtml.append(insertSeriesBegin());
    qul_cnt = state1->GetQulCnt();
    for (int i = 0; i < 42; i++)
        tmp_double[i] = qul_cnt[i];
    outhtml.append(insertSeriesData("line", tmp_double, 42));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(MeanQuality2));




    //AGCT Content
    outhtml.append(insertChart(PositionContent1));
    //option
    outhtml.append(insertOptionBegin(PositionContent1));
    outhtml.append(insertTitle("AGCT Content(Before filtering)"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertLegend("\'A\',\'G\',\'C\',\'T\'"));
    outhtml.append(insertxAxis("Position in read(bp)", mx_len + 1));
    outhtml.append(insertyAxis("value", std::to_string(0), std::to_string(1)));
    outhtml.append(insertSeriesBegin());

    pos_cnt_ = state1->GetPosCnt();
    for (int i = 0; i < mx_len; i++) {
        int64_t sum_tot = 0.0;
        for (int j = 0; j < 8; j++) sum_tot += pos_cnt_[i * 8 + j];
        tmp_double[i] = 1.0 * pos_cnt_[i * 8 + 'A' & 0x07] / sum_tot;
    }
    outhtml.append(insertSeriesData("line", "A", tmp_double, mx_len));
    for (int i = 0; i < mx_len; i++) {
        int64_t sum_tot = 0.0;
        for (int j = 0; j < 8; j++) sum_tot += pos_cnt_[i * 8 + j];
        tmp_double[i] = 1.0 * pos_cnt_[i * 8 + 'G' & 0x07] / sum_tot;
    }
    outhtml.append(insertSeriesData("line", "G", tmp_double, mx_len));
    for (int i = 0; i < mx_len; i++) {
        int64_t sum_tot = 0.0;
        for (int j = 0; j < 8; j++) sum_tot += pos_cnt_[i * 8 + j];
        tmp_double[i] = 1.0 * pos_cnt_[i * 8 + 'C' & 0x07] / sum_tot;
    }
    outhtml.append(insertSeriesData("line", "C", tmp_double, mx_len));
    for (int i = 0; i < mx_len; i++) {
        int64_t sum_tot = 0.0;
        for (int j = 0; j < 8; j++) sum_tot += pos_cnt_[i * 8 + j];
        tmp_double[i] = 1.0 * pos_cnt_[i * 8 + 'T' & 0x07] / sum_tot;
    }
    outhtml.append(insertSeriesData("line", "T", tmp_double, mx_len));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(PositionContent1));


    outhtml.append(insertChart(PositionContent2));
    //option
    outhtml.append(insertOptionBegin(PositionContent2));
    outhtml.append(insertTitle("AGCT Content(After filtering)"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertLegend("\'A\',\'G\',\'C\',\'T\'"));
    outhtml.append(insertxAxis("Position in read(bp)", mx_len + 1));
    outhtml.append(insertyAxis("value", std::to_string(0), std::to_string(1)));
    outhtml.append(insertSeriesBegin());

    pos_cnt_ = state2->GetPosCnt();
    for (int i = 0; i < mx_len; i++) {
        int64_t sum_tot = 0.0;
        for (int j = 0; j < 8; j++) sum_tot += pos_cnt_[i * 8 + j];
        tmp_double[i] = 1.0 * pos_cnt_[i * 8 + 'A' & 0x07] / sum_tot;
    }
    outhtml.append(insertSeriesData("line", "A", tmp_double, mx_len));
    for (int i = 0; i < mx_len; i++) {
        int64_t sum_tot = 0.0;
        for (int j = 0; j < 8; j++) sum_tot += pos_cnt_[i * 8 + j];
        tmp_double[i] = 1.0 * pos_cnt_[i * 8 + 'G' & 0x07] / sum_tot;
    }
    outhtml.append(insertSeriesData("line", "G", tmp_double, mx_len));
    for (int i = 0; i < mx_len; i++) {
        int64_t sum_tot = 0.0;
        for (int j = 0; j < 8; j++) sum_tot += pos_cnt_[i * 8 + j];
        tmp_double[i] = 1.0 * pos_cnt_[i * 8 + 'C' & 0x07] / sum_tot;
    }
    outhtml.append(insertSeriesData("line", "C", tmp_double, mx_len));
    for (int i = 0; i < mx_len; i++) {
        int64_t sum_tot = 0.0;
        for (int j = 0; j < 8; j++) sum_tot += pos_cnt_[i * 8 + j];
        tmp_double[i] = 1.0 * pos_cnt_[i * 8 + 'T' & 0x07] / sum_tot;
    }
    outhtml.append(insertSeriesData("line", "T", tmp_double, mx_len));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(PositionContent2));



    //GC Content
    outhtml.append(insertChart(GCContent1));
    //option
    outhtml.append(insertOptionBegin(GCContent1));
    outhtml.append(insertTitle("Per Sequence GC content(Before filtering)"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("Mean GC content(%)", mx_len + 1, 1));
    outhtml.append(insertyAxis("value"));
    outhtml.append(insertSeriesBegin());
    int64_t *gc_cnt = state1->GetGcCnt();
    for (int i = 0; i < mx_len; i++)tmp_double[i] = gc_cnt[i];
    outhtml.append(insertSeriesData("line", tmp_double, mx_len));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(GCContent1));


    outhtml.append(insertChart(GCContent2));
    //option
    outhtml.append(insertOptionBegin(GCContent2));
    outhtml.append(insertTitle("Per Sequence GC content(After filtering)"));
    outhtml.append(insertTooltip());
    outhtml.append(insertDataZoom());
    outhtml.append(insertxAxis("Mean GC content(%)", mx_len + 1, 1));
    outhtml.append(insertyAxis("value"));
    outhtml.append(insertSeriesBegin());
    gc_cnt = state2->GetGcCnt();
    for (int i = 0; i < mx_len; i++)tmp_double[i] = gc_cnt[i];
    outhtml.append(insertSeriesData("line", tmp_double, mx_len));
    outhtml.append(insertSeriesEnd());
    outhtml.append(insertOptionEnd());
    outhtml.append(insertChartOption(GCContent2));
//
//    //N Content
//    outhtml.append(insertChart(NContent));
//    //option
//    outhtml.append(insertOptionBegin(NContent));
//    outhtml.append(insertTitle("Per Sequence N content"));
//    outhtml.append(insertTooltip());
//    outhtml.append(insertDataZoom());
//    outhtml.append(insertxAxis("position in read(%)", max_len));
//    outhtml.append(insertyAxis("value", std::to_string(0), std::to_string(1)));
//    outhtml.append(insertSeriesBegin());
//    for (int i = 0; i < max_len + 1; i++) {
//        double total = 0.0;
//        for (int j = 0; j < 8; j++) total += NumberList[i][j];
//        tmp_double[i] = NumberList[i]['N' & 0x07] / total;
//    }
//    outhtml.append(insertSeriesData("line", tmp_double, max_len));
//    outhtml.append(insertSeriesEnd());
//    outhtml.append(insertOptionEnd());
//    outhtml.append(insertChartOption(NContent));


    outhtml.append("</script>");
    outhtml.append("</html>");
    std::fstream fout = std::fstream("RabbitQC.html", std::ios::out | std::ios::binary);
    fout.write(outhtml.c_str(), outhtml.length());
    fout.close();

}