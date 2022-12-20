
#include <cstring>
#include <vector>
#include <unordered_map>
#include <unistd.h>

#include "tuple_spawn.hpp"
#include "slave.h"
#include "Reference.h"
#include "cmdinfo.h"
#include "state.h"
#include "tgsstate.h"
#include "threadinfo.h"
#define gcMax 100
extern "C" void* ldm_malloc(size_t size);
extern "C" void ldm_free(void *addr, size_t size);
struct qc_data {
    ThreadInfo **thread_info_;
    CmdInfo *cmd_info_;
    std::vector <neoReference> *data_;
    std::vector <neoReference> *pass_data_;
    int *cnt;
};

void updateTop5(int rlen, double avgQual, TGSStats *tgs_state) {
    int pos = 5;
    for (int i = 0; i < 5; i++) {
        if (tgs_state->mTop5QualReads[i].first < avgQual) {
            for (int j = 4; j > i; j--) {
                tgs_state->mTop5QualReads[j] = tgs_state->mTop5QualReads[j - 1];
            }
            tgs_state->mTop5QualReads[i] = {avgQual, rlen};
            break;
        }
    }
    for (int i = 0; i < 5; i++) {
        if (tgs_state->mTop5LengReads[i].first < rlen) {
            for (int j = 4; j > i; j--) {
                tgs_state->mTop5LengReads[j] = tgs_state->mTop5LengReads[j - 1];
            }
            tgs_state->mTop5LengReads[i] = {rlen, avgQual};
            break;
        }
    }
}

void tgsfunc(qc_data *para) {
    std::vector <neoReference> *data = para->data_;
    ThreadInfo *thread_info = para->thread_info_[_PEN];
    CmdInfo *cmd_info = para->cmd_info_;
    int data_num = data->size();
    int base_num = 0;
    // copy data to ldm
    int size_range_ = thread_info->TGS_state_->mMinlen >> 1;
    int64_t *head_qual_sum = (int64_t *)ldm_malloc(size_range_ * sizeof(int64_t));
    int64_t *tail_qual_sum = (int64_t *)ldm_malloc(size_range_ * sizeof(int64_t));
    dma_getn(thread_info->TGS_state_->head_qual_sum, head_qual_sum, size_range_);
    dma_getn(thread_info->TGS_state_->tail_qual_sum, tail_qual_sum, size_range_);

    int64_t *head_seq_pos_count[4];
    int64_t *tail_seq_pos_count[4];
    for(int i = 0; i < 4; i++){
        head_seq_pos_count[i] = (int64_t*)ldm_malloc(size_range_ * sizeof(int64_t));
        dma_getn(thread_info->TGS_state_->head_seq_pos_count[i], head_seq_pos_count[i], size_range_);
        tail_seq_pos_count[i] = (int64_t*)ldm_malloc(size_range_ * sizeof(int64_t));
        dma_getn(thread_info->TGS_state_->tail_seq_pos_count[i], tail_seq_pos_count[i], size_range_);
    }
    int64_t mBases51015Num[3];
    mBases51015Num[0] = thread_info->TGS_state_->mBases51015Num[0];
    mBases51015Num[1] = thread_info->TGS_state_->mBases51015Num[1];
    mBases51015Num[2] = thread_info->TGS_state_->mBases51015Num[2];
    for (int id = _PEN; id < data_num; id += 64) {
        auto ref = (*data)[id];
        const int rlen = ref.lseq;
        const char *seq = reinterpret_cast<const char *>(ref.base + ref.pseq);
        const char *quality = reinterpret_cast<const char *>(ref.base + ref.pqual);
        int size_range = size_range_;
        if(rlen < size_range) size_range = rlen;
        if(thread_info->TGS_state_->countLenNum < (1 << 20))
            thread_info->TGS_state_->mTotalReadsLen[thread_info->TGS_state_->countLenNum++] = rlen;
        thread_info->TGS_state_->mReadsNum++;
        thread_info->TGS_state_->mBasesNum += rlen;
        char c1, c2;
        int phredSub = 33;
        if (cmd_info->isPhred64_) phredSub = 64;
        int sumQual = 0;
        for (int i = 0; i < rlen; i++) {
            int qual = (quality[i] - phredSub);
            sumQual += qual;
            if (qual >= 5) mBases51015Num[0]++;
            if (qual >= 10) mBases51015Num[1]++;
            if (qual >= 15) mBases51015Num[2]++;
        }
        updateTop5(rlen, 1.0 * sumQual / rlen, thread_info->TGS_state_);
        for (int i = 0; i < size_range; ++i) {
            c1 = seq[i];
            head_qual_sum[i] += (quality[i] - phredSub);
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
        for (int i = rlen - 1; i >= rlen - size_range; --i) {
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
    dma_putn(thread_info->TGS_state_->head_qual_sum, head_qual_sum, size_range_);
    dma_putn(thread_info->TGS_state_->tail_qual_sum, tail_qual_sum, size_range_);
    ldm_free(head_qual_sum, size_range_ * sizeof(int64_t));
    ldm_free(tail_qual_sum, size_range_ * sizeof(int64_t));
    for(int i = 0; i < 4; i++){
        dma_putn(thread_info->TGS_state_->head_seq_pos_count[i], head_seq_pos_count[i], size_range_);
        dma_putn(thread_info->TGS_state_->tail_seq_pos_count[i], tail_seq_pos_count[i], size_range_);
        ldm_free(head_seq_pos_count[i], size_range_ * sizeof(int64_t));
        ldm_free(tail_seq_pos_count[i], size_range_ * sizeof(int64_t));
    }

    dma_putn(thread_info->TGS_state_->mBases51015Num, mBases51015Num, 3);
    //printf("=== %d base number %d\n", _PEN, base_num);
}



void StateInfo(State *state, neoReference &ref, CmdInfo *cmd_info_, int64_t *pos_cnt_, int64_t *pos_qul_, int64_t *len_cnt_, int64_t *gc_cnt_, int64_t *qul_cnt_, int64_t &q20bases_, int64_t &q30bases_, int64_t &tot_bases_, int64_t &gc_bases_, int &real_seq_len_, int64_t &lines_){
    int slen = ref.lseq;
    int qlen = ref.lqual;
    //ASSERT(slen == qlen);
    //if (slen > malloc_seq_len_) {
    //    ExtendBuffer(malloc_seq_len_, max(slen + 100, slen * 2));
    //}
    //TODO
    if(slen > real_seq_len_) real_seq_len_ = slen;
    len_cnt_[slen - 1]++;
    char *bases = reinterpret_cast<char *>(ref.base + ref.pseq);
    char *quals = reinterpret_cast<char *>(ref.base + ref.pqual);
    int flag = 4;
    int gc_cnt = 0;
    int qul_tot = 0;
    int kmer = 0;
    int phredSub = 33;
    if (cmd_info_->isPhred64_) phredSub = 64;
    for (int i = 0; i < slen; i++) {
        char b = bases[i] & 0x07;
        if (b == 3 || b == 7) gc_cnt++;
        int q = quals[i] - phredSub;
        if(q < 0) q = 0;
        qul_tot += q;
        if (q >= 30) {
            q20bases_++;
            q30bases_++;
        } else if (q >= 20) {
            q20bases_++;
        }
        pos_cnt_[i * 8 + b]++;
        pos_qul_[i * 8 + b] += q;
        //if (bases[i] == 'N') flag = 5;
        //int val = valAGCT[b];
        //kmer = ((kmer << 2) & 0x3FC) | val;
        //if (flag <= 0) state->kmer_[kmer]++;
        //flag--;
    }
    tot_bases_ += slen;
    gc_bases_ += gc_cnt;
    gc_cnt_[int(1.0 * gcMax * gc_cnt / slen)]++;
    qul_cnt_[int(1.0 * qul_tot / slen)]++;
    // do overrepresentation analysis for 1 of every 20 reads
    //if (do_over_represent_analyze_) {
    //    if (lines_ % over_representation_sampling_ == 0) {
    //        //StateORP(ref);
    //        //TODO
    //    }
    //}
    lines_++;
}

bool TrimSeq(neoReference &ref, int front, int tail, CmdInfo * cmd_info) {
    int new_seq_len = ref.lseq - front - tail;
    if (new_seq_len <= 0) return false;

    int w = cmd_info->cut_window_size_;
    int l = ref.lseq;
    const char *seq = reinterpret_cast<const char *>(ref.base + ref.pseq);
    const char *qualstr = reinterpret_cast<const char *>(ref.base + ref.pqual);
    int phredSub = 33;
    if (cmd_info->isPhred64_) phredSub = 64;
    // quality cutting forward
    if (cmd_info->trim_5end_) {
        int s = front;
        if (l - front - tail - w <= 0) {
            return 0;
        }

        int totalQual = 0;

        // preparing rolling
        for (int i = 0; i < w - 1; i++)
            totalQual += qualstr[s + i];

        for (s = front; s + w < l - tail; s++) {
            totalQual += qualstr[s + w - 1];
            // rolling
            if (s > front) {
                totalQual -= qualstr[s - 1];
            }
            // add 33 for phred33 transforming
            if ((double) totalQual / (double) w >= phredSub + cmd_info->cut_mean_quality_)
                break;
        }

        // the trimming in front is forwarded and rlen is recalculated
        if (s > 0)
            s = s + w - 1;
        while (s < l && seq[s] == 'N')
            s++;
        front = s;
        new_seq_len = l - front - tail;
    }

    // quality cutting backward
    if (cmd_info->trim_3end_) {
        if (l - front - tail - w <= 0) {
            return 0;
        }

        int totalQual = 0;
        int t = l - tail - 1;

        // preparing rolling
        for (int i = 0; i < w - 1; i++)
            totalQual += qualstr[t - i];

        for (t = l - tail - 1; t - w >= front; t--) {
            totalQual += qualstr[t - w + 1];
            // rolling
            if (t < l - tail - 1) {
                totalQual -= qualstr[t + 1];
            }
            // add 33 for phred33 transforming
            if ((double) totalQual / (double) w >= phredSub + cmd_info->cut_mean_quality_)
                break;
        }

        if (t < l - 1)
            t = t - w + 1;
        while (t >= 0 && seq[t] == 'N')
            t--;
        new_seq_len = t - front + 1;
    }

    if (new_seq_len <= 0 || front >= l - 1) {
        return 0;
    }


    ref.pseq += front;
    ref.pqual += front;
    ref.lseq = new_seq_len;
    ref.lqual = new_seq_len;

    return true;

}
int ReadFiltering(neoReference &ref, bool trim_res, bool isPhred64, CmdInfo * cmd_info) {
    if (!trim_res) return 2;
    int seq_len = ref.lseq;
    int qul_len = ref.lqual;
    //ASSERT(seq_len == qul_len);
    if (seq_len >= cmd_info->length_limit_) {
        return 3;
    }
    if (seq_len == 0 || seq_len < cmd_info->length_required_) {
        return 2;
    }
    int n_number = 0;
    int low_qual_number = 0;// < q15
    char *bases = reinterpret_cast<char *>(ref.base + ref.pseq);
    char *quals = reinterpret_cast<char *>(ref.base + ref.pqual);

    int phredSub = 33;
    if (isPhred64) phredSub = 64;
    for (int i = 0; i < seq_len; i++) {
        int q = std::max(0, quals[i] - phredSub);
        if (q < 15) {
            low_qual_number++;
        }
        if (bases[i] == 'N') {
            n_number++;
        }
    }
    if (low_qual_number > cmd_info->low_qual_perc_limit_ * seq_len / 100.0) return 4;
    if (n_number > cmd_info->n_number_limit_) return 1;
    return 0;
}



void ngsfunc(qc_data *para){
    std::vector <neoReference> *data = para->data_;
    std::vector <neoReference> *pass_data = para->pass_data_;
    ThreadInfo *thread_info = para->thread_info_[_PEN];
    CmdInfo *cmd_info_ = para->cmd_info_;
    int data_num = data->size();
    int pre_pos_seq_len = thread_info->pre_state1_->malloc_seq_len_ * 8;
    int64_t *pre_pos_cnt_ = (int64_t*)ldm_malloc(pre_pos_seq_len * sizeof(int64_t));
    int64_t *pre_pos_qul_ = (int64_t*)ldm_malloc(pre_pos_seq_len * sizeof(int64_t));
    dma_getn(thread_info->pre_state1_->pos_cnt_, pre_pos_cnt_, pre_pos_seq_len);
    dma_getn(thread_info->pre_state1_->pos_qul_, pre_pos_qul_, pre_pos_seq_len);

    int aft_pos_seq_len = thread_info->aft_state1_->malloc_seq_len_ * 8;
    int64_t *aft_pos_cnt_ = (int64_t*)ldm_malloc(aft_pos_seq_len * sizeof(int64_t));
    int64_t *aft_pos_qul_ = (int64_t*)ldm_malloc(aft_pos_seq_len * sizeof(int64_t));
    dma_getn(thread_info->aft_state1_->pos_cnt_, aft_pos_cnt_, aft_pos_seq_len);
    dma_getn(thread_info->aft_state1_->pos_qul_, aft_pos_qul_, aft_pos_seq_len);


    int pre_malloc_seq_len_ = thread_info->pre_state1_->malloc_seq_len_;
    int64_t *pre_len_cnt_ = (int64_t*)ldm_malloc(pre_malloc_seq_len_ * sizeof(int64_t));
    dma_getn(thread_info->pre_state1_->len_cnt_, pre_len_cnt_, pre_malloc_seq_len_);

    int aft_malloc_seq_len_ = thread_info->aft_state1_->malloc_seq_len_;
    int64_t *aft_len_cnt_ = (int64_t*)ldm_malloc(aft_malloc_seq_len_ * sizeof(int64_t));
    dma_getn(thread_info->aft_state1_->len_cnt_, aft_len_cnt_, aft_malloc_seq_len_);

    int64_t *pre_gc_cnt_ = (int64_t*)ldm_malloc((gcMax + 100) * sizeof(int64_t));
    dma_getn(thread_info->pre_state1_->gc_cnt_, pre_gc_cnt_, gcMax + 1);

    int64_t *aft_gc_cnt_ = (int64_t*)ldm_malloc((gcMax + 100) * sizeof(int64_t));
    dma_getn(thread_info->aft_state1_->gc_cnt_, aft_gc_cnt_, gcMax + 1);

    int pre_qul_range_ = thread_info->pre_state1_->qul_range_;
    int64_t *pre_qul_cnt_ = (int64_t*)ldm_malloc(pre_qul_range_ * sizeof(int64_t));
    dma_getn(thread_info->pre_state1_->qul_cnt_, pre_qul_cnt_, pre_qul_range_);

    int aft_qul_range_ = thread_info->aft_state1_->qul_range_;
    int64_t *aft_qul_cnt_ = (int64_t*)ldm_malloc(aft_qul_range_ * sizeof(int64_t));
    dma_getn(thread_info->aft_state1_->qul_cnt_, aft_qul_cnt_, aft_qul_range_);

    int64_t pre_q20bases_ = 0;
    int64_t pre_q30bases_ = 0;
    int64_t pre_tot_bases_ = 0;
    int64_t pre_gc_bases_ = 0;
    int pre_real_seq_len_ = 0;
    int64_t pre_lines_ = 0;

    int64_t aft_q20bases_ = 0;
    int64_t aft_q30bases_ = 0;
    int64_t aft_tot_bases_ = 0;
    int64_t aft_gc_bases_ = 0;
    int aft_real_seq_len_ = 0;
    int64_t aft_lines_ = 0;

    for(int id = _PEN; id < data_num; id += 64){
        auto item = (*data)[id];
        StateInfo(thread_info->pre_state1_, item, cmd_info_, pre_pos_cnt_, pre_pos_qul_, pre_len_cnt_, pre_gc_cnt_, pre_qul_cnt_, pre_q20bases_, pre_q30bases_, pre_tot_bases_, pre_gc_bases_, pre_real_seq_len_, pre_lines_); 
        if (cmd_info_->state_duplicate_) {
            //StateDup(item);
            //TODO
        }

        if (cmd_info_->add_umi_) {
            //umier_->ProcessSe(item);
            //TODO
        }
        bool trim_res = 0;
        trim_res = TrimSeq(item, cmd_info_->trim_front2_, cmd_info_->trim_tail1_, cmd_info_);

        if (trim_res && cmd_info_->trim_polyg_) {
            //PolyX::trimPolyG(item, cmd_info_->trim_poly_len_);
            //TODO
        }

        if (trim_res && cmd_info_->trim_polyx_) {
            //PolyX::trimPolyX(item, cmd_info_->trim_poly_len_);
            //TODO
        }

        if (trim_res && cmd_info_->trim_adapter_) {
            int res = 0;
            bool is_trimmed = false;
            if (cmd_info_->detect_adapter1_) {
                if (cmd_info_->print_what_trimmed_) {
                    //res = TrimAdapter(item, cmd_info_->adapter_seq1_, thread_info->aft_state1_->adapter_map_, cmd_info_->adapter_len_lim_, false);
                    //TODO
                } else {
                    //res = TrimAdapter(item, cmd_info_->adapter_seq1_, false);
                    //TODO
                }
                if (res) {
                    is_trimmed = true;
                    thread_info->aft_state1_->trim_adapter_++;
                    thread_info->aft_state1_->trim_adapter_bases_ += res;
                }
            }
            if (cmd_info_->adapter_from_fasta_.size() > 0) {
                if (cmd_info_->print_what_trimmed_) {
                    //res = TrimAdapters(item, cmd_info_->adapter_from_fasta_, thread_info->aft_state1_->adapter_map_, cmd_info_->adapter_len_lim_, false);
                    //TODO
                } else {
                    //res = TrimAdapters(item, cmd_info_->adapter_from_fasta_, false);
                    //TODO
                }
                if (res) {
                    if (!is_trimmed) thread_info->aft_state1_->trim_adapter_++;
                    thread_info->aft_state1_->trim_adapter_bases_ += res;
                }
            }
        }

        int filter_res = 0;
        filter_res = ReadFiltering(item, trim_res, cmd_info_->isPhred64_, cmd_info_);
        if (filter_res == 0) {
            StateInfo(thread_info->aft_state1_, item, cmd_info_, aft_pos_cnt_, aft_pos_qul_, aft_len_cnt_, aft_gc_cnt_, aft_qul_cnt_, aft_q20bases_, aft_q30bases_, aft_tot_bases_, aft_gc_bases_, aft_real_seq_len_, aft_lines_); 
            thread_info->aft_state1_->pass_reads_++;
            if (cmd_info_->write_data_) {
                (*pass_data)[id] = item;
            }
        } else {
            if (cmd_info_->write_data_) {
                item.lname = 0;
                item.lseq = 0;
                item.lstrand = 0;
                item.lqual = 0;
                (*pass_data)[id] = item;
            }
            if (filter_res == 1) {
                thread_info->aft_state1_->fail_N_++;
            } else if (filter_res == 2) {
                thread_info->aft_state1_->fail_short_++;
            } else if (filter_res == 3) {
                thread_info->aft_state1_->fail_long_++;
            } else if (filter_res == 4) {
                thread_info->aft_state1_->fail_lowq_++;
            }
        }
    }
    dma_putn(thread_info->pre_state1_->pos_cnt_, pre_pos_cnt_, pre_pos_seq_len);
    dma_putn(thread_info->pre_state1_->pos_qul_, pre_pos_qul_, pre_pos_seq_len);
    dma_putn(thread_info->aft_state1_->pos_cnt_, aft_pos_cnt_, aft_pos_seq_len);
    dma_putn(thread_info->aft_state1_->pos_qul_, aft_pos_qul_, aft_pos_seq_len);
    ldm_free(pre_pos_cnt_, pre_pos_seq_len * sizeof(int64_t));
    ldm_free(pre_pos_qul_, pre_pos_seq_len * sizeof(int64_t));
    ldm_free(aft_pos_cnt_, aft_pos_seq_len * sizeof(int64_t));
    ldm_free(aft_pos_qul_, aft_pos_seq_len * sizeof(int64_t));

    dma_putn(thread_info->pre_state1_->len_cnt_, pre_len_cnt_, pre_malloc_seq_len_);
    dma_putn(thread_info->aft_state1_->len_cnt_, aft_len_cnt_, aft_malloc_seq_len_);
    dma_putn(thread_info->pre_state1_->gc_cnt_, pre_gc_cnt_, gcMax + 1);
    dma_putn(thread_info->aft_state1_->gc_cnt_, aft_gc_cnt_, gcMax + 1);
    dma_putn(thread_info->pre_state1_->qul_cnt_, pre_qul_cnt_, pre_qul_range_);
    dma_putn(thread_info->aft_state1_->qul_cnt_, aft_qul_cnt_, aft_qul_range_);
    ldm_free(pre_len_cnt_, pre_malloc_seq_len_ * sizeof(int64_t));
    ldm_free(aft_len_cnt_, aft_malloc_seq_len_ * sizeof(int64_t));
    ldm_free(pre_gc_cnt_, (gcMax + 1) * sizeof(int64_t));
    ldm_free(aft_gc_cnt_, (gcMax + 1) * sizeof(int64_t));
    ldm_free(pre_qul_cnt_, pre_qul_range_ * sizeof(int64_t));
    ldm_free(aft_qul_cnt_, aft_qul_range_ * sizeof(int64_t));

    thread_info->pre_state1_->q20bases_ += pre_q20bases_;
    thread_info->pre_state1_->q30bases_ += pre_q30bases_;
    thread_info->pre_state1_->tot_bases_ += pre_tot_bases_;
    thread_info->pre_state1_->gc_bases_ += pre_gc_bases_;
    thread_info->pre_state1_->real_seq_len_ = pre_real_seq_len_;
    thread_info->pre_state1_->lines_ += pre_lines_;

    thread_info->aft_state1_->q20bases_ += aft_q20bases_;
    thread_info->aft_state1_->q30bases_ += aft_q30bases_;
    thread_info->aft_state1_->tot_bases_ += aft_tot_bases_;
    thread_info->aft_state1_->gc_bases_ += aft_gc_bases_;
    thread_info->aft_state1_->real_seq_len_ = aft_real_seq_len_;
    thread_info->aft_state1_->lines_ += aft_lines_;
}


ATHREAD_VISIBLE(tgsfunc);
//ATHREAD_VISIBLE(ngsfunc);



