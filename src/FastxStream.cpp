/*
   This file is a part of DSRC software distributed under GNU GPL 2 licence.
   The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc

Authors: Lucas Roguski and Sebastian Deorowicz

Version: 2.00

last modified by Zekun Yin 2020/5/18
*/

#include <iostream>
#include <string>

#include "FastxStream.h"
#include "Reference.h"

#ifdef Vec512

#include "immintrin.h"

#endif


#ifdef Vec256

#include "immintrin.h"

#endif

namespace rabbit {

    namespace fa {

        /**
          @brief Read the next chunk
          @return FastaChunk pointer if next chunk data has data, else return NULL
          */
        FastaChunk *FastaFileReader::readNextChunk() {
            FastaDataChunk *part = NULL;
            recordsPool.Acquire(part);
            FastaChunk *dataPart = new FastaChunk;
            dataPart->chunk = part;
            if (ReadNextChunk_(dataPart, this->seqInfos)) {
                return dataPart;
            } else {
                recordsPool.Release(part);
                return NULL;
            }
        }

        /**
         * @brief Read next listed chunk
         * @details this function make sure one FastaChunk(dataPart) contains at least
         * one whole sequence
         * @return FastaChunk pointer if next chunk data has data, else return NULL
         */
        FastaChunk *FastaFileReader::readNextChunkList() {
            FastaDataChunk *part = NULL;
            recordsPool.Acquire(part);
            FastaChunk *dataPart = new FastaChunk;
            dataPart->chunk = part;
            FastaDataChunk *current = part;
            bool continue_read = false;
            dataPart->start = this->totalSeqs;
            if (ReadNextFaChunk_(dataPart->chunk, this->seqInfos, continue_read)) {
                //FastaDataChunk *current = part;
                while (continue_read) {
                    FastaDataChunk *append = NULL;
                    recordsPool.Acquire(append);
                    if (ReadNextFaChunk_(append, this->seqInfos, continue_read)) {
                        current->next = append;
                        current = append;
                    } else {
                        recordsPool.Release(append);
                        break;
                    }
                }
                return dataPart;
            } else {
                recordsPool.Release(part);
                return NULL;
            }
        }

        bool FastaFileReader::ReadNextChunk_(FastaChunk *dataChunk_, SeqInfos &seqInfos) {
            FastaDataChunk *chunk_ = dataChunk_->chunk;
            if (Eof()) {
                chunk_->size = 0;
                return false;
            }

            // flush the data from previous incomplete chunk
            uchar *data = chunk_->data.Pointer();
            const uint64 cbufSize = chunk_->data.Size();
            chunk_->size = 0;
            int64 toRead = cbufSize - bufferSize;// buffersize: size left from last chunk

            if (bufferSize > 0) {
                std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
                chunk_->size = bufferSize;
                bufferSize = 0;
            }

            // read the next chunk
            int64 r = this->Read(data + chunk_->size, toRead);
            // std::cout << "r is :" << r << std::endl;
            // std::cout << "toRead: " << toRead << std::endl;

            if (r > 0) {
                if (r == toRead)// somewhere before end
                {
                    // dealing with halo region
                    uint64 chunkEnd = FindCutPos_(dataChunk_, data, cbufSize, mHalo, seqInfos);
                    chunk_->size = chunkEnd;// - 1; //1 char back from last '>'

                    if (usesCrlf) chunk_->size -= 1;

                    // copy tail to swapBuffer
                    int haloCount = 0;
                    uint64 tmpPtr = chunkEnd - 1;
                    if (mHalo > 0) {
                        while (true) {
                            if (data[tmpPtr] != '\n') {
                                haloCount++;
                                // std::cerr << (char)data[tmpPtr] << std::endl;
                                if (haloCount == mHalo) break;
                            }
                            tmpPtr--;
                        }
                        // std::cerr << "haloCount: " << haloCount << std::endl;
                    } else {
                        tmpPtr = chunkEnd;// nocopy
                    }
                    std::copy(data + tmpPtr, data + cbufSize, swapBuffer.Pointer());
                    bufferSize = cbufSize - tmpPtr;

                } else// at the end of file
                {
                    chunk_->size += r - 1;// skip the last EOF symbol
                    if (usesCrlf) chunk_->size -= 1;

                    // only for get seqsinfo
                    uint64 chunkEnd = FindCutPos_(dataChunk_, data, chunk_->size, mHalo, seqInfos);
                    // debug only
                    // std::string content((char*)data, chunk_->size);
                    // std::cout << "tail content ori: " << data << std::endl;
                    // std::cout << "tail content    : " << content << std::endl;
                    // end debug

                    eof = true;
                }
            } else {
                eof = true;
            }

            return true;
        }

        uint64
        FastaFileReader::FindCutPos_(FastaChunk *dataChunk_, uchar *data_, const uint64 size_, const uint64 halo_,
                                     SeqInfos &seqInfos) {
            int count = 0;
            uint64 pos_ = 0;
            uint64 cut_ = 0;     // cut_ point to next '>'
            uint64 lastSeq_ = 0; //-> the start of last sequences content
            uint64 lastName_ = 0;//-> the last '>'
            OneSeqInfo seqInfo;

            //dataChunk_ -> start and dataChunk_ -> end means start index and end index of current chunk
            if (data_[0] == '>')// start with '>'
            {
                dataChunk_->start = this->totalSeqs;
                while (pos_ < size_) {
                    if (data_[pos_] == '>') {
                        lastName_ = pos_;
                        if (FindEol(data_, pos_, size_))// find name
                        {
                            ++pos_;
                            lastSeq_ = pos_;

                            seqInfo.gid = this->totalSeqs;
                            seqInfos.push_back(seqInfo);

                            this->totalSeqs++;
                        } else {
                            cut_ = pos_;// find a cut: incomplete name
                            // std::cerr << "cut char: " << (char)data_[cut_] << std::endl;
                            break;
                        }
                    } else {
                        ++pos_;
                    }
                }

            } else {// start with {ACGT}
                dataChunk_->start = this->totalSeqs - 1;
                while (pos_ < size_) {
                    if (data_[pos_] == '>') {
                        lastName_ = pos_;
                        if (FindEol(data_, pos_, size_))// find name
                        {
                            ++pos_;
                            lastSeq_ = pos_;
                            seqInfo.gid = this->totalSeqs;
                            seqInfos.push_back(seqInfo);
                            this->totalSeqs++;
                        } else {
                            cut_ = pos_;// find a cut -> '>'
                            // std::cout << "cut char: " << (char)data_[cut_] << std::endl;
                            break;
                        }
                    } else {
                        ++pos_;
                    }
                }
            }

            // no tail cut
            // make sure that cbufSize > name_len + halo
            if (cut_ == 0) {
                uint64 lastSeqLen_ = size_ - lastSeq_;
                if (lastSeqLen_ < halo_) {
                    cut_ = lastName_;
                    this->totalSeqs--;
                }
            }

            dataChunk_->nseqs = this->totalSeqs - dataChunk_->start;
            dataChunk_->end = this->totalSeqs - 1;

            // if(cut_ != 0) std::cout << "cut: " << cut_ << std::endl;

            return cut_ ? cut_ : size_;
        }

        bool find_next_seq_start(uchar *data, uint64 size, uint64 &pos_) {
            if (pos_ == size - 1) return false;
            do {
                pos_++;
            } while (pos_ < size && data[pos_] != '>');

            return data[pos_] == '>' ? true : false;
        }

        bool FastaFileReader::ReadNextFaChunk_(FastaDataChunk *chunk_, SeqInfos &seqInfos, bool &continue_read) {
            if (Eof()) {
                chunk_->size = 0;
                return false;
            }

            // flush the data from previous incomplete chunk
            uchar *data = chunk_->data.Pointer();
            const uint64 cbufSize = chunk_->data.Size();
            chunk_->size = 0;
            int64 toRead = cbufSize - bufferSize;// buffersize: size left from last chunk

            if (bufferSize > 0) {
                std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
                chunk_->size = bufferSize;
                bufferSize = 0;
            }

            // read the next chunk
            int64 r = this->Read(data + chunk_->size, toRead);
            // std::cout << "r is :" << r << std::endl;
            // std::cout << "toRead: " << toRead << std::endl;

            if (r > 0) {
                if (r == toRead)// somewhere before end
                {
                    if (data[0] == '>') {
                        uint64 chunkEnd = 0, pos_ = 0;
                        //pos_ = chunk_->size;
                        while (find_next_seq_start(data, cbufSize, pos_)) {
                            chunkEnd = pos_;
                            this->totalSeqs++;
                        }
                        if (chunkEnd == 0) {// means this chunk is a big one, one chunk can not contain all
                            continue_read = true;
                            chunkEnd = cbufSize;
                        } else {
                            continue_read = false;
                        }
                        chunk_->size = chunkEnd - 1;
                        if (usesCrlf) chunk_->size -= 1;
                        std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
                        bufferSize = cbufSize - chunkEnd;
                    } else {// means this is a continue chunk for last sequence
                        uint64 chunkEnd = 0, pos_ = 0;
                        pos_ = chunk_->size;
                        if (find_next_seq_start(data, cbufSize, pos_)) {
                            chunkEnd = pos_;
                            this->totalSeqs++;
                        }
                        if (chunkEnd == 0) {// means this chunk is a big one, one chunk can not contain all
                            continue_read = true;
                            chunkEnd = cbufSize;
                        } else {
                            continue_read = false;
                        }
                        chunk_->size = chunkEnd - 1;
                        if (usesCrlf) chunk_->size -= 1;
                        std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
                        bufferSize = cbufSize - chunkEnd;
                    }
                } else {                  // at the end of file
                    chunk_->size += r - 1;// skip the last eof symbol
                    if (usesCrlf) chunk_->size -= 1;
                    eof = true;
                }
            } else {
                eof = true;
            }

            return true;
        }

    }// namespace fa

    namespace fq {

        int64 count_line(uchar *contenx, int64 read_bytes) {
            int64 count_n = 0;
#ifdef Vec512
            __m512i enter_con = _mm512_set1_epi8('\n');
            for (int i = 0; i < read_bytes; i += 64) {
                uint64_t tag = (1ll << std::min(64, int(read_bytes) - i)) - 1;
                __m512i item = _mm512_maskz_loadu_epi8(tag, contenx + i);
                __mmask64 mask = _mm512_cmp_epi8_mask(item, enter_con, _MM_CMPINT_EQ);
                count_n += _mm_popcnt_u64(mask);
            }
#elif Vec256
            int i = 0;
            __m256i enter_con = _mm256_set1_epi8('\n');

            for (; i + 32 <= read_bytes; i += 32) {
                __m256i item = _mm256_loadu_si256((__m256i const *)(contenx + i));
                __m256i res = _mm256_cmpeq_epi8(item, enter_con);
                unsigned mask = _mm256_movemask_epi8(res);
                count_n += _mm_popcnt_u32(mask);
            }
            for (; i < read_bytes; ++i) {
                if (contenx[i] == '\n') count_n++;
            }
#else

            for (int i = 0; i < read_bytes; ++i) {
                if (contenx[i] == '\n') count_n++;
            }
#endif

            return count_n;
        }

        void FastqFileReader::readChunk() {
            FastqDataChunk *part = NULL;

            recordsPool.Acquire(part);
            printf("fastqio: ready to into while\n");

            // while (!errorHandler.IsError() && fileReader.ReadNextChunk(part))
            while (ReadNextChunk_(part)) {
                ASSERT(part->size > 0);

                //printf("numParts is %d\n", numParts);
                numParts++;

                recordsPool.Release(part);
                recordsPool.Acquire(part);
            }

            ASSERT(part->size == 0);
            recordsPool.Release(part);// the last empty part

            // recordsQueue.SetCompleted();
        }

        /**
          @brief Read the next chunk
          @return FastqChunk pointer if next chunk data has data, else return NULL
          */

        FastqDataChunk *FastqFileReader::readNextChunkFromMem(int64 offset, int64 lim_size) {
            FastqDataChunk *part = NULL;
            recordsPool.Acquire(part);
            if (ReadNextChunkFromMem_(part, offset, lim_size)) {
                return part;
            } else {
                recordsPool.Release(part);
                return NULL;
            }
        }
        FastqDataChunk *FastqFileReader::readNextChunk(int64 offset, int64 lim_size) {
            FastqDataChunk *part = NULL;
            recordsPool.Acquire(part);
            if (ReadNextChunk_(part, offset, lim_size)) {
                return part;
            } else {
                recordsPool.Release(part);
                return NULL;
            }
        }

        FastqDataChunk *FastqFileReader::readNextInterChunk() {
            FastqDataChunk *part = NULL;
            recordsPool.Acquire(part);
            if (ReadNextInterChunk_(part)) {
                return part;
            } else {
                recordsPool.Release(part);
                return NULL;
            }
        }


        FastqDataChunk *FastqFileReader::readNextChunk(moodycamel::ReaderWriterQueue<std::pair<char *, int>> *q, std::atomic_int *d, std::pair<char *, int> &l) {
            FastqDataChunk *part = NULL;
            recordsPool.Acquire(part);
            if (ReadNextChunk_(part, q, d, l)) {
                return part;
            } else {
                recordsPool.Release(part);
                return NULL;
            }
        }

        FastqDataChunk *FastqFileReader::readNextPairChunkInterleaved() {
            FastqDataChunk *part = NULL;
            recordsPool.Acquire(part);
            uint64 lastChunkPos1;
            uint64 lastChunkPos2;
            //---------read chunk------------
            if (eof) {
                part->size = 0;
                // return false;
                recordsPool.Release(part);
                return NULL;
            }

            // flush the data from previous incomplete chunk
            uchar *data = part->data.Pointer();
            const uint64 cbufSize = part->data.Size();
            part->size = 0;
            int64 toRead;
            toRead = cbufSize - bufferSize;
            if (bufferSize > 0) {
                std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
                part->size = bufferSize;
                bufferSize = 0;
            }
            int64 r;
            r = mFqReader->Read(data + part->size, toRead);
            if (r > 0) {
                if (r == toRead) {
                    lastChunkPos1 = GetPreviousRecordPos_(data, cbufSize - 1, cbufSize);
                    lastChunkPos2 = GetPreviousRecordPos_(data, lastChunkPos1 - 1, cbufSize);

                } else {
                    // chunkEnd = r;
                    part->size += r - 1;
                    if (usesCrlf) part->size -= 1;
                    eof = true;
                }
            } else {
                eof = true;
                return NULL;
            }
            //------read chunk end------//

            if (!eof) {

                part->size = lastChunkPos1 - 1;
                if (usesCrlf) part->size -= 1;

                std::copy(data + lastChunkPos2, data + cbufSize, swapBuffer.Pointer());
                bufferSize = cbufSize - lastChunkPos2;
            }
            return part;
        }

        //add parallel in this fastq pe file read funciton


        FastqDataPairChunk *
        FastqFileReader::readNextPairChunkParallel(moodycamel::ReaderWriterQueue<std::pair<char *, int>> *q1,
                                                   moodycamel::ReaderWriterQueue<std::pair<char *, int>> *q2,
                                                   std::atomic_int *d1, std::atomic_int *d2,
                                                   std::pair<char *, int> &last1, std::pair<char *, int> &last2) {
            bool eof1 = false;
            bool eof2 = false;
            FastqDataPairChunk *pair = new FastqDataPairChunk;


            FastqDataChunk *leftPart = NULL;
            recordsPool.Acquire(leftPart);

            FastqDataChunk *rightPart = NULL;
            recordsPool.Acquire(rightPart);

            if (eof) {
                leftPart->size = 0;
                rightPart->size = 0;
                // return false;
                recordsPool.Release(leftPart);
                recordsPool.Release(rightPart);
                return NULL;
            }


            int64 left_line_count = 0;
            int64 chunkEnd = 0;
            uchar *data = NULL;
            uint64 cbufSize = 0;
            int64 r = 0;
            std::thread threadp1([&]() {
                data = leftPart->data.Pointer();
                cbufSize = leftPart->data.Size();
                leftPart->size = 0;
                int64 toRead = cbufSize - bufferSize;
                if (bufferSize > 0) {
                    std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
                    leftPart->size = bufferSize;
                    bufferSize = 0;
                }
                r = mFqReader->Read(data + leftPart->size, toRead, q1, d1, last1);
                if (r > 0) {
                    if (r == toRead) {
                        chunkEnd = cbufSize - GetNxtBuffSize;
                        chunkEnd = this->GetNextRecordPos_(data, chunkEnd, cbufSize);
                    } else {
                        leftPart->size += r - 1;
                        if (usesCrlf) leftPart->size -= 1;
                        eof1 = true;
                        chunkEnd = leftPart->size + 1;
                        cbufSize = leftPart->size + 1;
                    }
                } else {
                    eof1 = true;
                    chunkEnd = leftPart->size + 1;
                    cbufSize = leftPart->size + 1;
                    //if (leftPart->size == 0) {
                    //    return NULL;
                    //}
                }
            });
            //right part
            int64 right_line_count = 0;
            int64 chunkEnd_right = 0;
            uchar *data_right = NULL;
            uint64 cbufSize_right = 0;
            int64 r_right = 0;
            std::thread threadp2([&]() {
                data_right = rightPart->data.Pointer();
                cbufSize_right = rightPart->data.Size();
                rightPart->size = 0;
                int64 toRead_right = cbufSize_right - bufferSize2;
                if (bufferSize2 > 0) {
                    std::copy(swapBuffer2.Pointer(), swapBuffer2.Pointer() + bufferSize2, data_right);
                    rightPart->size = bufferSize2;
                    bufferSize2 = 0;
                }
                r_right = mFqReader2->Read(data_right + rightPart->size, toRead_right, q2, d2, last2);
                if (r_right > 0) {
                    if (r_right == toRead_right) {
                        chunkEnd_right = cbufSize_right - GetNxtBuffSize;
                        chunkEnd_right = this->GetNextRecordPos_(data_right, chunkEnd_right, cbufSize_right);
                    } else {
                        rightPart->size += r_right - 1;
                        if (usesCrlf) rightPart->size -= 1;
                        eof2 = true;
                        chunkEnd_right = rightPart->size + 1;
                        cbufSize_right = rightPart->size + 1;
                    }
                } else {
                    eof2 = true;
                    chunkEnd_right = rightPart->size + 1;
                    cbufSize_right = rightPart->size + 1;
                    //if (rightPart->size == 0) {
                    //    return NULL;
                    //}
                }
            });
            threadp1.join();

            threadp2.join();

            if (r == 0 && leftPart->size == 0) return NULL;
            if (r_right == 0 && rightPart->size == 0) return NULL;

            if (eof1 && eof2)
                eof = true;
            if (!eof) {
                left_line_count = count_line(data, chunkEnd);
                right_line_count = count_line(data_right, chunkEnd_right);
                int64 difference = left_line_count - right_line_count;
                if (difference > 0) {
                    while (chunkEnd >= 0) {
                        if (data[chunkEnd] == '\n') {
                            difference--;
                            if (difference == -1) {
                                chunkEnd++;
                                break;
                            }
                        }
                        chunkEnd--;
                    }
                } else if (difference < 0) {
                    while (chunkEnd_right >= 0) {
                        if (data_right[chunkEnd_right] == '\n') {
                            difference++;
                            if (difference == 1) {
                                chunkEnd_right++;
                                break;
                            }
                        }
                        chunkEnd_right--;
                    }
                }
                leftPart->size = chunkEnd - 1;
                if (usesCrlf) leftPart->size -= 1;
                std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
                bufferSize = cbufSize - chunkEnd;
                rightPart->size = chunkEnd_right - 1;
                if (usesCrlf) rightPart->size -= 1;
                std::copy(data_right + chunkEnd_right, data_right + cbufSize_right, swapBuffer2.Pointer());
                bufferSize2 = cbufSize_right - chunkEnd_right;
            }
            pair->left_part = leftPart;
            pair->right_part = rightPart;
            return pair;
        }


        FastqDataPairChunk *FastqFileReader::readNextPairChunkParallel() {
            bool eof1 = false;
            bool eof2 = false;
            FastqDataPairChunk *pair = new FastqDataPairChunk;


            FastqDataChunk *leftPart = NULL;
            recordsPool.Acquire(leftPart);

            FastqDataChunk *rightPart = NULL;
            recordsPool.Acquire(rightPart);

            if (eof) {
                leftPart->size = 0;
                rightPart->size = 0;
                // return false;
                recordsPool.Release(leftPart);
                recordsPool.Release(rightPart);
                return NULL;
            }


            int64 left_line_count = 0;
            int64 chunkEnd = 0;
            uchar *data = NULL;
            uint64 cbufSize = 0;
            int64 r = 0;
            std::thread threadp1([&]() {
                data = leftPart->data.Pointer();
                cbufSize = leftPart->data.Size();
                leftPart->size = 0;
                int64 toRead = cbufSize - bufferSize;
                if (bufferSize > 0) {
                    std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
                    leftPart->size = bufferSize;
                    bufferSize = 0;
                }
                r = mFqReader->Read(data + leftPart->size, toRead);
                if (r > 0) {
                    if (r == toRead) {
                        chunkEnd = cbufSize - GetNxtBuffSize;
                        chunkEnd = this->GetNextRecordPos_(data, chunkEnd, cbufSize);
                    } else {
                        leftPart->size += r - 1;
                        if (usesCrlf) leftPart->size -= 1;
                        eof1 = true;
                        chunkEnd = leftPart->size + 1;
                        cbufSize = leftPart->size + 1;
                    }
                } else {
                    eof1 = true;
                    chunkEnd = leftPart->size + 1;
                    cbufSize = leftPart->size + 1;
                    //if (leftPart->size == 0) {
                    //    return NULL;
                    //}
                }
            });
            //right part
            int64 right_line_count = 0;
            int64 chunkEnd_right = 0;
            uchar *data_right = NULL;
            uint64 cbufSize_right = 0;
            int64 r_right = 0;
            std::thread threadp2([&]() {
                data_right = rightPart->data.Pointer();
                cbufSize_right = rightPart->data.Size();
                rightPart->size = 0;
                int64 toRead_right = cbufSize_right - bufferSize2;
                if (bufferSize2 > 0) {
                    std::copy(swapBuffer2.Pointer(), swapBuffer2.Pointer() + bufferSize2, data_right);
                    rightPart->size = bufferSize2;
                    bufferSize2 = 0;
                }
                r_right = mFqReader2->Read(data_right + rightPart->size, toRead_right);
                if (r_right > 0) {
                    if (r_right == toRead_right) {
                        chunkEnd_right = cbufSize_right - GetNxtBuffSize;
                        chunkEnd_right = this->GetNextRecordPos_(data_right, chunkEnd_right, cbufSize_right);
                    } else {
                        rightPart->size += r_right - 1;
                        if (usesCrlf) rightPart->size -= 1;
                        eof2 = true;
                        chunkEnd_right = rightPart->size + 1;
                        cbufSize_right = rightPart->size + 1;
                    }
                } else {
                    eof2 = true;
                    chunkEnd_right = rightPart->size + 1;
                    cbufSize_right = rightPart->size + 1;
                    //if (rightPart->size == 0) {
                    //    return NULL;
                    //}
                }
            });
            threadp1.join();
            threadp2.join();
            if (r == 0 && leftPart->size == 0) return NULL;
            if (r_right == 0 && rightPart->size == 0) return NULL;

            if (eof1 && eof2) eof = true;
            if (!eof) {
                left_line_count = count_line(data, chunkEnd);
                right_line_count = count_line(data_right, chunkEnd_right);
                int64 difference = left_line_count - right_line_count;
                if (difference > 0) {
                    while (chunkEnd >= 0) {
                        if (data[chunkEnd] == '\n') {
                            difference--;
                            if (difference == -1) {
                                chunkEnd++;
                                break;
                            }
                        }
                        chunkEnd--;
                    }
                } else if (difference < 0) {
                    while (chunkEnd_right >= 0) {
                        if (data_right[chunkEnd_right] == '\n') {
                            difference++;
                            if (difference == 1) {
                                chunkEnd_right++;
                                break;
                            }
                        }
                        chunkEnd_right--;
                    }
                }
                leftPart->size = chunkEnd - 1;
                if (usesCrlf) leftPart->size -= 1;
                std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
                bufferSize = cbufSize - chunkEnd;
                rightPart->size = chunkEnd_right - 1;
                if (usesCrlf) rightPart->size -= 1;
                std::copy(data_right + chunkEnd_right, data_right + cbufSize_right, swapBuffer2.Pointer());
                bufferSize2 = cbufSize_right - chunkEnd_right;
            }
            pair->left_part = leftPart;
            pair->right_part = rightPart;
            return pair;
        }


        /**
  @brief Read the next paired chunk
  @return FastqDataPairChunk pointer if next chunk data has data, else return NULL
  */

        FastqDataPairChunk *FastqFileReader::readNextPairChunkFromMem(int64 offset, int64 lim_size) {
            //qwer

            //printf("offset limsize : %lld %lld\n", offset, lim_size);
            
            static int64 now_size = 0;
            static int64 now_size2 = 0;
            bool eof1 = false;
            bool eof2 = false;
            FastqDataPairChunk *pair = new FastqDataPairChunk;

            FastqDataChunk *leftPart = NULL;
            recordsPool.Acquire(leftPart);

            FastqDataChunk *rightPart = NULL;
            recordsPool.Acquire(rightPart);

            int64 left_line_count = 0;
            int64 right_line_count = 0;
            int64 chunkEnd = 0;
            int64 chunkEnd_right = 0;

            //---------read left chunk------------
            if (eof) {
                leftPart->size = 0;
                rightPart->size = 0;
                // return false;
                recordsPool.Release(leftPart);
                recordsPool.Release(rightPart);
                return NULL;
            }

            // flush the data from previous incomplete chunk
            uchar *data = leftPart->data.Pointer();
            uint64 cbufSize = leftPart->data.Size();
            leftPart->size = 0;
            int64 toRead;
            toRead = cbufSize - bufferSize;
            if (bufferSize > 0) {
                std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
                leftPart->size = bufferSize;
                bufferSize = 0;
            }
            if(now_size + toRead > lim_size) {
                toRead = lim_size - now_size;
            }
            int64 r;
            //if(offset == -1 || offset == -2) r = mFqReader->Read(data + leftPart->size, toRead);
            if(offset == -1 || offset == -2) r = ReadFromMem(data + leftPart->size, toRead);
            else {
                //r = mFqReader->ReadSeek(data + leftPart->size, toRead, offset);
                r = ReadSeekFromMem(data + leftPart->size, toRead, offset);
            }
            if(offset != -2) now_size += r;
            if (r > 0) {
                if (r == toRead && now_size < lim_size) {
                    chunkEnd = cbufSize - GetNxtBuffSize;
                    chunkEnd = GetNextRecordPos_(data, chunkEnd, cbufSize);
                } else {
                    leftPart->size += r - 1;
                    if (usesCrlf) leftPart->size -= 1;
                    eof1 = true;
                    chunkEnd = leftPart->size + 1;
                    cbufSize = leftPart->size + 1;
                }
            } else {
                eof1 = true;
                chunkEnd = leftPart->size + 1;
                cbufSize = leftPart->size + 1;
                if (leftPart->size == 0) {
                    return NULL;
                }
            }
            //------read left chunk end------//

            //-----------------read right chunk---------------------//
            uchar *data_right = rightPart->data.Pointer();
            uint64 cbufSize_right = rightPart->data.Size();
            rightPart->size = 0;
            toRead = cbufSize_right - bufferSize2;
            if (bufferSize2 > 0) {
                std::copy(swapBuffer2.Pointer(), swapBuffer2.Pointer() + bufferSize2, data_right);
                rightPart->size = bufferSize2;
                bufferSize2 = 0;
            }
            if(now_size2 + toRead > lim_size) {
                toRead = lim_size - now_size2;
            }
            //if(offset == -1 || offset == -2) r = mFqReader2->Read(data_right + rightPart->size, toRead);
            if(offset == -1 || offset == -2) r = ReadFromMem2(data_right + rightPart->size, toRead);
            else {
                //r = mFqReader2->ReadSeek(data_right + rightPart->size, toRead, offset);
                r = ReadSeekFromMem2(data_right + rightPart->size, toRead, offset);
            
            }
            if(offset != -2) now_size2 += r;

            if (r > 0) {
                if (r == toRead && now_size2 < lim_size) {
                    chunkEnd_right = cbufSize_right - GetNxtBuffSize;
                    chunkEnd_right = GetNextRecordPos_(data_right, chunkEnd_right, cbufSize_right);
                } else {
                    // chunkEnd_right += r;
                    rightPart->size += r - 1;
                    if (usesCrlf) rightPart->size -= 1;
                    eof2 = true;
                    chunkEnd_right = rightPart->size + 1;
                    cbufSize_right = rightPart->size + 1;
                }
            } else {
                eof2 = true;
                chunkEnd_right = rightPart->size + 1;
                cbufSize_right = rightPart->size + 1;
                if (rightPart->size == 0) {
                    return NULL;
                }
            }
            //--------------read right chunk end---------------------//
            if (eof1 && eof2) eof = true;
            if (!eof) {
                left_line_count = count_line(data, chunkEnd);
                right_line_count = count_line(data_right, chunkEnd_right);
                int64 difference = left_line_count - right_line_count;
                if(difference) printf("difference.......\n");
                if (difference > 0) {
                    while (chunkEnd >= 0) {
                        if (data[chunkEnd] == '\n') {
                            difference--;
                            if (difference == -1) {
                                chunkEnd++;
                                break;
                            }
                        }
                        chunkEnd--;
                    }
                } else if (difference < 0) {
                    while (chunkEnd_right >= 0) {
                        if (data_right[chunkEnd_right] == '\n') {
                            difference++;
                            if (difference == 1) {
                                chunkEnd_right++;
                                break;
                            }
                        }
                        chunkEnd_right--;
                    }
                }
                leftPart->size = chunkEnd - 1;
                if (usesCrlf) leftPart->size -= 1;
                std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
                bufferSize = cbufSize - chunkEnd;
                rightPart->size = chunkEnd_right - 1;
                if (usesCrlf) rightPart->size -= 1;
                std::copy(data_right + chunkEnd_right, data_right + cbufSize_right, swapBuffer2.Pointer());
                bufferSize2 = cbufSize_right - chunkEnd_right;
            }
            pair->left_part = leftPart;
            pair->right_part = rightPart;
            return pair;
        }

        FastqDataPairChunk *FastqFileReader::readNextPairChunk(int64 offset, int64 lim_size) {
            if(isZipped) lim_size = 1ll << 60;
            //qwer

            //printf("offset limsize : %lld %lld\n", offset, lim_size);
            
            static int64 now_size = 0;
            static int64 now_size2 = 0;
            bool eof1 = false;
            bool eof2 = false;
            FastqDataPairChunk *pair = new FastqDataPairChunk;

            FastqDataChunk *leftPart = NULL;
            recordsPool.Acquire(leftPart);

            FastqDataChunk *rightPart = NULL;
            recordsPool.Acquire(rightPart);

            int64 left_line_count = 0;
            int64 right_line_count = 0;
            int64 chunkEnd = 0;
            int64 chunkEnd_right = 0;

            //---------read left chunk------------
            if (eof) {
                leftPart->size = 0;
                rightPart->size = 0;
                // return false;
                recordsPool.Release(leftPart);
                recordsPool.Release(rightPart);
                return NULL;
            }

            // flush the data from previous incomplete chunk
            uchar *data = leftPart->data.Pointer();
            uint64 cbufSize = leftPart->data.Size();
            leftPart->size = 0;
            int64 toRead;
            toRead = cbufSize - bufferSize;
            if (bufferSize > 0) {
                std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
                leftPart->size = bufferSize;
                bufferSize = 0;
            }
            if(now_size + toRead > lim_size) {
                toRead = lim_size - now_size;
            }
            int64 r;
            if(offset == -1 || offset == -2) r = mFqReader->Read(data + leftPart->size, toRead);
            else {
                //std::cerr << "seek_read " << offset << std::endl;
                r = mFqReader->ReadSeek(data + leftPart->size, toRead, offset);
                //uchar *pos_tmp = data + leftPart->size;
                //for(int i = 0; i < 40; i++) std::cerr << pos_tmp[i];
                //std::cerr << std::endl;
            
            }
            if(offset != -2) now_size += r;
            //printf("=== res %lld, %lld\n", r, now_size);
            if (r > 0) {
                if (r == toRead && now_size < lim_size) {
                    chunkEnd = cbufSize - GetNxtBuffSize;
                    chunkEnd = GetNextRecordPos_(data, chunkEnd, cbufSize);
                } else {
                    // chunkEnd = r;
                    leftPart->size += r - 1;
                    if (usesCrlf) leftPart->size -= 1;
                    eof1 = true;
                    chunkEnd = leftPart->size + 1;
                    cbufSize = leftPart->size + 1;
                }
            } else {
                eof1 = true;
                chunkEnd = leftPart->size + 1;
                cbufSize = leftPart->size + 1;
                if (leftPart->size == 0) {
                    return NULL;
                }
            }
            //------read left chunk end------//

            //-----------------read right chunk---------------------//
            uchar *data_right = rightPart->data.Pointer();
            uint64 cbufSize_right = rightPart->data.Size();
            rightPart->size = 0;
            toRead = cbufSize_right - bufferSize2;
            if (bufferSize2 > 0) {
                std::copy(swapBuffer2.Pointer(), swapBuffer2.Pointer() + bufferSize2, data_right);
                rightPart->size = bufferSize2;
                bufferSize2 = 0;
            }
            if(now_size2 + toRead > lim_size) {
                toRead = lim_size - now_size2;
            }
            if(offset == -1 || offset == -2) r = mFqReader2->Read(data_right + rightPart->size, toRead);
            else {
                //std::cerr << "seek_read " << offset << std::endl;
                r = mFqReader2->ReadSeek(data_right + rightPart->size, toRead, offset);
                //uchar *pos_tmp = data_right + rightPart->size;
                //for(int i = 0; i < 40; i++) std::cerr << pos_tmp[i];
                //std::cerr << std::endl;
            
            }
            if(offset != -2) now_size2 += r;

            //printf("=== res %lld, %lld\n", r, now_size2);
            if (r > 0) {
                if (r == toRead && now_size2 < lim_size) {
                    chunkEnd_right = cbufSize_right - GetNxtBuffSize;
                    chunkEnd_right = GetNextRecordPos_(data_right, chunkEnd_right, cbufSize_right);
                } else {
                    // chunkEnd_right += r;
                    rightPart->size += r - 1;
                    if (usesCrlf) rightPart->size -= 1;
                    eof2 = true;
                    chunkEnd_right = rightPart->size + 1;
                    cbufSize_right = rightPart->size + 1;
                }
            } else {
                eof2 = true;
                chunkEnd_right = rightPart->size + 1;
                cbufSize_right = rightPart->size + 1;
                if (rightPart->size == 0) {
                    return NULL;
                }
            }
            //--------------read right chunk end---------------------//
            if (eof1 && eof2) eof = true;
            if (!eof) {
                left_line_count = count_line(data, chunkEnd);
                right_line_count = count_line(data_right, chunkEnd_right);
                int64 difference = left_line_count - right_line_count;
                if(difference) printf("difference.......\n");
                if (difference > 0) {
                    while (chunkEnd >= 0) {
                        if (data[chunkEnd] == '\n') {
                            difference--;
                            if (difference == -1) {
                                chunkEnd++;
                                break;
                            }
                        }
                        chunkEnd--;
                    }
                } else if (difference < 0) {
                    while (chunkEnd_right >= 0) {
                        if (data_right[chunkEnd_right] == '\n') {
                            difference++;
                            if (difference == 1) {
                                chunkEnd_right++;
                                break;
                            }
                        }
                        chunkEnd_right--;
                    }
                }
                leftPart->size = chunkEnd - 1;
                if (usesCrlf) leftPart->size -= 1;
                std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
                bufferSize = cbufSize - chunkEnd;
                rightPart->size = chunkEnd_right - 1;
                if (usesCrlf) rightPart->size -= 1;
                std::copy(data_right + chunkEnd_right, data_right + cbufSize_right, swapBuffer2.Pointer());
                bufferSize2 = cbufSize_right - chunkEnd_right;
            }
            pair->left_part = leftPart;
            pair->right_part = rightPart;
            return pair;
        }

        FastqDataPairChunk *
        FastqFileReader::readNextPairChunk(moodycamel::ReaderWriterQueue<std::pair<char *, int>> *q1,
                                           moodycamel::ReaderWriterQueue<std::pair<char *, int>> *q2,
                                           std::atomic_int *d1, std::atomic_int *d2,
                                           std::pair<char *, int> &last1, std::pair<char *, int> &last2) {
            bool eof1 = false;
            bool eof2 = false;
            FastqDataPairChunk *pair = new FastqDataPairChunk;

            FastqDataChunk *leftPart = NULL;
            recordsPool.Acquire(leftPart);// song: all in one datapool

            FastqDataChunk *rightPart = NULL;
            recordsPool.Acquire(rightPart);

            int64 left_line_count = 0;
            int64 right_line_count = 0;
            int64 chunkEnd = 0;
            int64 chunkEnd_right = 0;

            //---------read left chunk------------
            if (eof) {
                leftPart->size = 0;
                rightPart->size = 0;
                // return false;
                recordsPool.Release(leftPart);
                recordsPool.Release(rightPart);
                return NULL;
            }

            // flush the data from previous incomplete chunk
            uchar *data = leftPart->data.Pointer();
            uint64 cbufSize = leftPart->data.Size();
            leftPart->size = 0;
            int64 toRead;
            toRead = cbufSize - bufferSize;
            if (bufferSize > 0) {
                std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
                leftPart->size = bufferSize;
                bufferSize = 0;
            }
            int64 r;
            r = mFqReader->Read(data + leftPart->size, toRead, q1, d1, last1);
            if (r > 0) {
                if (r == toRead) {
                    chunkEnd = cbufSize - GetNxtBuffSize;
                    chunkEnd = GetNextRecordPos_(data, chunkEnd, cbufSize);
                } else {
                    leftPart->size += r - 1;
                    if (usesCrlf) leftPart->size -= 1;
                    eof1 = true;
                    chunkEnd = leftPart->size + 1;
                    cbufSize = leftPart->size + 1;
                }
            } else {
                eof1 = true;
                chunkEnd = leftPart->size + 1;
                cbufSize = leftPart->size + 1;
                if (leftPart->size == 0) {
                    return NULL;
                }
            }

            uchar *data_right = rightPart->data.Pointer();
            uint64 cbufSize_right = rightPart->data.Size();
            rightPart->size = 0;
            toRead = cbufSize_right - bufferSize2;
            if (bufferSize2 > 0) {
                std::copy(swapBuffer2.Pointer(), swapBuffer2.Pointer() + bufferSize2, data_right);
                rightPart->size = bufferSize2;
                bufferSize2 = 0;
            }
            r = mFqReader2->Read(data_right + rightPart->size, toRead, q2, d2, last2);
            if (r > 0) {
                if (r == toRead) {
                    chunkEnd_right = cbufSize_right - GetNxtBuffSize;
                    chunkEnd_right = GetNextRecordPos_(data_right, chunkEnd_right, cbufSize_right);
                } else {
                    rightPart->size += r - 1;
                    if (usesCrlf) rightPart->size -= 1;
                    eof2 = true;
                    chunkEnd_right = rightPart->size + 1;
                    cbufSize_right = rightPart->size + 1;
                }
            } else {
                eof2 = true;
                chunkEnd_right = rightPart->size + 1;
                cbufSize_right = rightPart->size + 1;
                if (rightPart->size == 0) {
                    return NULL;
                }
            }
            if (eof1 && eof2)
                eof = true;
            if (!eof) {
                left_line_count = count_line(data, chunkEnd);
                right_line_count = count_line(data_right, chunkEnd_right);
                int64 difference = left_line_count - right_line_count;
                if (difference > 0) {
                    while (chunkEnd >= 0) {
                        if (data[chunkEnd] == '\n') {
                            difference--;
                            if (difference == -1) {
                                chunkEnd++;
                                break;
                            }
                        }
                        chunkEnd--;
                    }
                } else if (difference < 0) {
                    while (chunkEnd_right >= 0) {
                        if (data_right[chunkEnd_right] == '\n') {
                            difference++;
                            if (difference == 1) {
                                chunkEnd_right++;
                                break;
                            }
                        }
                        chunkEnd_right--;
                    }
                }

                leftPart->size = chunkEnd - 1;
                if (usesCrlf) leftPart->size -= 1;
                std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
                bufferSize = cbufSize - chunkEnd;

                rightPart->size = chunkEnd_right - 1;
                if (usesCrlf) rightPart->size -= 1;
                std::copy(data_right + chunkEnd_right, data_right + cbufSize_right, swapBuffer2.Pointer());
                bufferSize2 = cbufSize_right - chunkEnd_right;
            }
            pair->left_part = leftPart;
            pair->right_part = rightPart;
            return pair;
        }

        void FastqFileReader::MemDataReader() {
            byte* tmpData = new byte[SwapBufferSize];
            while(true) {
                int64 r = mFqReader->Read(tmpData, SwapBufferSize);
                MemDataTotSize += r;
                if(r < SwapBufferSize) break;
            }
            //fprintf(stderr, "tot size is %lld\n", MemDataTotSize);
            MemData = new byte[MemDataTotSize];
            mFqReader->ReadSeek(MemData, MemDataTotSize, 0);
            //fprintf(stderr, "read to mem done\n");
        }


        void FastqFileReader::MemDataReader2() {
            byte* tmpData = new byte[SwapBufferSize];
            while(true) {
                int64 r = mFqReader2->Read(tmpData, SwapBufferSize);
                MemDataTotSize2 += r;
                if(r < SwapBufferSize) break;
            }
            //fprintf(stderr, "tot size2 is %lld\n", MemDataTotSize2);
            MemData2 = new byte[MemDataTotSize2];
            mFqReader2->ReadSeek(MemData2, MemDataTotSize2, 0);
            //fprintf(stderr, "read to mem2 done\n");
        }

        void FastqFileReader::ReleaseMemData() {
            delete MemData;
            
        }


        void FastqFileReader::ReleaseMemData2() {
            delete MemData2;
            
        }


        int64 FastqFileReader::ReadSeekFromMem(byte *memory_, uint64 size_, uint64 pos_) {
            MemDataNowPos = pos_;
            int64 lastDataSize = MemDataTotSize - MemDataNowPos;
            if(size_ > lastDataSize) {
               memcpy(memory_, MemData + MemDataNowPos, lastDataSize);
               MemDataNowPos += lastDataSize;
               MemDataReadFinish = 1;
               return lastDataSize;
            } else {
                memcpy(memory_, MemData + MemDataNowPos, size_);
                MemDataNowPos += size_;
                return size_;
            }

        }

        int64 FastqFileReader::ReadSeekFromMem2(byte *memory_, uint64 size_, uint64 pos_) {
            MemDataNowPos2 = pos_;
            int64 lastDataSize2 = MemDataTotSize2 - MemDataNowPos2;
            if(size_ > lastDataSize2) {
               memcpy(memory_, MemData2 + MemDataNowPos2, lastDataSize2);
               MemDataNowPos2 += lastDataSize2;
               MemDataReadFinish2 = 1;
               return lastDataSize2;
            } else {
                memcpy(memory_, MemData2 + MemDataNowPos2, size_);
                MemDataNowPos2 += size_;
                return size_;
            }

        }

        int64 FastqFileReader::ReadFromMem(byte *memory_, uint64 size_) {
            int64 lastDataSize = MemDataTotSize - MemDataNowPos;
            //fprintf(stderr, "read from mem %lld / %lld / %lld\n", MemDataTotSize, MemDataNowPos, size_);
            if(size_ > lastDataSize) {
               memcpy(memory_, MemData + MemDataNowPos, lastDataSize);
               MemDataNowPos += lastDataSize;
               MemDataReadFinish = 1;
               return lastDataSize;
            } else {
                memcpy(memory_, MemData + MemDataNowPos, size_);
                MemDataNowPos += size_;
                return size_;
            }

        }


        int64 FastqFileReader::ReadFromMem2(byte *memory_, uint64 size_) {
            int64 lastDataSize2 = MemDataTotSize2 - MemDataNowPos2;
            //fprintf(stderr, "read from mem2 %lld / %lld / %lld\n", MemDataTotSize2, MemDataNowPos2, size_);
            if(size_ > lastDataSize2) {
               memcpy(memory_, MemData2 + MemDataNowPos2, lastDataSize2);
               MemDataNowPos2 += lastDataSize2;
               MemDataReadFinish2 = 1;
               return lastDataSize2;
            } else {
                memcpy(memory_, MemData2 + MemDataNowPos2, size_);
                MemDataNowPos2 += size_;
                return size_;
            }

        }

        bool FastqFileReader::FinishReadFromMem() {
            return MemDataReadFinish;
        }
 
        bool FastqFileReader::FinishReadFromMem2() {
            return MemDataReadFinish2;
        }
            

        //tagg
          
        bool FastqFileReader::ReadNextChunkFromMem_(FastqDataChunk *chunk_, int64 offset, int64 lim_size) {
            
            static int64 now_size = 0;
            if ((FinishReadFromMem() || now_size == lim_size) && bufferSize == 0) {
                chunk_->size = 0;
                return false;
            }

            // flush the data from previous incomplete chunk
            uchar *data = chunk_->data.Pointer();
            uint64 cbufSize = chunk_->data.Size();
            chunk_->size = 0;
            int64 toRead = cbufSize - bufferSize;
            //---------------------------------
            if (bufferSize > 0) {
                std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
                chunk_->size = bufferSize;
                bufferSize = 0;
            }

            // read the next chunk
            if(now_size + toRead > lim_size) {
                toRead = lim_size - now_size;
            }
            //printf("read %lld %lld %lld\n", toRead, offset, lim_size);
            int64 r = 0;
            //if(offset == -1 || offset == -2) r = mFqReader->Read(data + chunk_->size, toRead);
            if(offset == -1 || offset == -2) r = ReadFromMem(data + chunk_->size, toRead);
            else {
                //std::cerr << "seek_read " << offset << " " << chunk_->size << std::endl;
                //r = mFqReader->ReadSeek(data + chunk_->size, toRead, offset);
                r = ReadSeekFromMem(data + chunk_->size, toRead, offset);
                //uchar *pos_tmp = data + chunk_->size;
                //for(int i = 0; i < 40; i++) std::cerr << pos_tmp[i];
                //std::cerr << std::endl;
            }
            if(offset != -2) now_size += r;
            //printf("=== res %lld, %lld\n", r, now_size);
            //if (!mFqReader->FinishRead() && now_size < lim_size) {
            if (!FinishReadFromMem() && now_size < lim_size) {
                cbufSize = r + chunk_->size;
                uint64 chunkEnd = cbufSize - (cbufSize < GetNxtBuffSize ? cbufSize : GetNxtBuffSize);
                chunkEnd = GetNextRecordPos_(data, chunkEnd, cbufSize);
                chunk_->size = chunkEnd - 1;
                if (usesCrlf) chunk_->size -= 1;

                std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
                bufferSize = cbufSize - chunkEnd;
            } else {                  // at the end of file
                chunk_->size += r - 1;// skip the last EOF symbol
                if (usesCrlf) chunk_->size -= 1;
                //mFqReader->setEof();
            }
            return true;
        }

 
        bool FastqFileReader::ReadNextChunk_(FastqDataChunk *chunk_, int64 offset, int64 lim_size) {
            if(isZipped) lim_size = 1ll << 60;
            
            static int64 now_size = 0;
            if ((mFqReader->FinishRead() || now_size == lim_size) && bufferSize == 0) {
                chunk_->size = 0;
                return false;
            }

            // flush the data from previous incomplete chunk
            uchar *data = chunk_->data.Pointer();
            uint64 cbufSize = chunk_->data.Size();
            chunk_->size = 0;
            int64 toRead = cbufSize - bufferSize;
            //---------------------------------
            if (bufferSize > 0) {
                std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
                chunk_->size = bufferSize;
                bufferSize = 0;
            }

            // read the next chunk
            if(now_size + toRead > lim_size) {
                toRead = lim_size - now_size;
            }
            int64 r = 0;
            if(offset == -1 || offset == -2) r = mFqReader->Read(data + chunk_->size, toRead);
            else {
                //std::cerr << "seek_read " << offset << " " << chunk_->size << std::endl;
                r = mFqReader->ReadSeek(data + chunk_->size, toRead, offset);
                //uchar *pos_tmp = data + chunk_->size;
                //for(int i = 0; i < 40; i++) std::cerr << pos_tmp[i];
                //std::cerr << std::endl;
            }
            if(offset != -2) now_size += r;
            //printf("=== res %lld, %lld\n", r, now_size);
            if (!mFqReader->FinishRead() && now_size < lim_size) {
                cbufSize = r + chunk_->size;
                uint64 chunkEnd = cbufSize - (cbufSize < GetNxtBuffSize ? cbufSize : GetNxtBuffSize);
                chunkEnd = GetNextRecordPos_(data, chunkEnd, cbufSize);
                chunk_->size = chunkEnd - 1;
                if (usesCrlf) chunk_->size -= 1;

                std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
                bufferSize = cbufSize - chunkEnd;
            } else {                  // at the end of file
                chunk_->size += r - 1;// skip the last EOF symbol
                if (usesCrlf) chunk_->size -= 1;
                mFqReader->setEof();
            }
            return true;
        }

        bool FastqFileReader::ReadNextChunk_(FastqDataChunk *chunk_, moodycamel::ReaderWriterQueue<std::pair<char *, int>> *q, std::atomic_int *d, std::pair<char *, int> &l) {
            if ((mFqReader->FinishRead() || mFqReader->Eof()) && bufferSize == 0) {
                chunk_->size = 0;
                return false;
            }

            // flush the data from previous incomplete chunk
            uchar *data = chunk_->data.Pointer();
            uint64 cbufSize = chunk_->data.Size();
            chunk_->size = 0;
            int64 toRead = cbufSize - bufferSize;
            //---------------------------------
            if (bufferSize > 0) {
                std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
                chunk_->size = bufferSize;
                bufferSize = 0;
            }

            // read the next chunk
            int64 r = mFqReader->Read(data + chunk_->size, toRead, q, d, l);
            // std::cout << "r is :" << r << std::endl;
            // std::cout << "r is :" << r << std::endl;

            // if (r > 0) {
            if (!mFqReader->FinishRead() && !mFqReader->Eof()) {
                cbufSize = r + chunk_->size;
                uint64 chunkEnd = cbufSize - (cbufSize < GetNxtBuffSize ? cbufSize : GetNxtBuffSize);
                chunkEnd = GetNextRecordPos_(data, chunkEnd, cbufSize);
                chunk_->size = chunkEnd - 1;
                if (usesCrlf) chunk_->size -= 1;

                std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
                bufferSize = cbufSize - chunkEnd;
            } else {                  // at the end of file
                chunk_->size += r - 1;// skip the last EOF symbol
                if (usesCrlf) chunk_->size -= 1;
                mFqReader->setEof();
            }
            if (r != toRead) {
                mFqReader->setEof();
            }
            //}else {
            //	mFqReader->setEof();
            //}
            return true;
        }


        bool FastqFileReader::ReadNextInterChunk_(FastqDataChunk *chunk_) {
            if (mFqReader->FinishRead() && bufferSize == 0) {
                chunk_->size = 0;
                return false;
            }

            uchar *data = chunk_->data.Pointer();
            uint64 cbufSize = chunk_->data.Size();
            uint64 chunkEnd = 0;
            chunk_->size = 0;
            int64 toRead = cbufSize - bufferSize;
            if (bufferSize > 0) {
                std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
                chunk_->size = bufferSize;
                bufferSize = 0;
            }

            int64 r = mFqReader->Read(data + chunk_->size, toRead);

            if (!mFqReader->FinishRead()) {
                cbufSize = r + chunk_->size;
                chunkEnd = cbufSize - (cbufSize < GetNxtBuffSize ? cbufSize : GetNxtBuffSize);
                chunkEnd = GetNextRecordPos_(data, chunkEnd, cbufSize);
            } else {                  // at the end of file
                chunk_->size += r - 1;// skip the last EOF symbol
                if (usesCrlf) chunk_->size -= 1;
                mFqReader->setEof();
            }
            if(!mFqReader->FinishRead()){
                int line_cnout = count_line(data, chunkEnd);
                ASSERT(line_cnout % 4 == 0);
                int item_count = line_cnout / 4;
                if(item_count & 1){
                    int difference = 4;
                    while (chunkEnd >= 0){
                        if (data[chunkEnd] == '\n'){
                            difference--;
                            if (difference == -1) {
                                chunkEnd++;
                                break;
                            }
                        }
                        chunkEnd--;
                    }
                    ASSERT(difference == -1);
                }
                chunk_->size = chunkEnd - 1;
                if (usesCrlf) chunk_->size -= 1;
                std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
                bufferSize = cbufSize - chunkEnd;
            }
            return true;
        }

        uint64 FastqFileReader::GetPreviousRecordPos_(uchar *data_, uint64 pos_, const uint64 size_) {
            int offset = 2;

            SkipToSol(data_, pos_, size_);
            if (usesCrlf) {
                offset = 3;
            }
            while (data_[pos_ + offset] != '@') {//+2
                // std::cout<<"pos_"<<pos_<<std::endl;
                SkipToSol(data_, pos_, size_);
            }
            //mark to check if the '@' is quality score
            uint64 pos0 = pos_ + offset;
            SkipToSol(data_, pos_, size_);
            if (data_[pos_ + offset] == '+') {
                // indicate that the '@' is quality score
                SkipToSol(data_, pos_, size_);
                SkipToSol(data_, pos_, size_);
                //it is name
                if (data_[pos_ + offset] != '@') {
                    std::cout << "core dump is " << data_[pos_ + offset] << std::endl;
                    return pos_ + offset;
                } else {
                    return pos_ + offset;
                }
            } else {
                return pos0;
            }
        }

        uint64 FastqFileReader::GetNextRecordPos_(uchar *data_, uint64 pos_, const uint64 size_) {
            SkipToEol(data_, pos_, size_);
            ++pos_;

            // find beginning of the next record
            while (data_[pos_] != '@') {
                SkipToEol(data_, pos_, size_);
                ++pos_;
            }
            uint64 pos0 = pos_;

            SkipToEol(data_, pos_, size_);
            ++pos_;

            if (data_[pos_] == '@')// previous one was a quality field
                return pos_;
            //-----[haoz:] is the following code necessary??-------------//
            SkipToEol(data_, pos_, size_);
            ++pos_;
            if (data_[pos_] != '+') std::cout << "core dump is pos: " << pos_ << " char: " << data_[pos_] << std::endl;
            ASSERT(data_[pos_] == '+');// pos0 was the start of tag
            return pos0;
        }

        // uint64 FastqFileReader::GetPrevRecordPos_(uchar *data_, uint64 pos_, const uint64 size_) {
        //	while(data_[pos_] != '@'){
        //		SkipToSol(data_, pos_, size_);
        //	}
        // }
        /**
 * @brief Skip to end of line
 * @param data_ base pointer
 * @param pos_ start position to skip
 * @param size_ data_ size
 */
        void FastqFileReader::SkipToEol(uchar *data_, uint64 &pos_, const uint64 size_) {
            // cerr << "pos: " << pos_ << " size: " << size_ << endl;
            ASSERT(pos_ < size_);

            while (data_[pos_] != '\n' && data_[pos_] != '\r' && pos_ < size_) ++pos_;

            if (data_[pos_] == '\r' && pos_ < size_) {
                if (data_[pos_ + 1] == '\n') {
                    usesCrlf = true;
                    ++pos_;
                }
            }
        }

        /**
 * @brief Skip to start of line
 * @param data_ base pointer
 * @param pos_ start position to skip
 * @param size_ data_ size
 */
        void FastqFileReader::SkipToSol(uchar *data_, uint64 &pos_, const uint64 size_) {
            ASSERT(pos_ < size_);
            if (data_[pos_] == '\n') {
                --pos_;
            }
            if (data_[pos_] == '\r') {
                usesCrlf = true;
                pos_--;
            }
            // find next '\n' or \n\r
            while (data_[pos_] != '\n' && data_[pos_] != '\r') {
                --pos_;
            }
            if (data_[pos_] == '\n') {
                --pos_;//?? is ++pos_ ?
            }
            if (data_[pos_] == '\r') {
                usesCrlf = true;
                pos_--;
            }
        }

    }// namespace fq

}// namespace rabbit
