/*
   This file is a part of DSRC software distributed under GNU GPL 2 licence.
   The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc

Authors: Lucas Roguski and Sebastian Deorowicz

Version: 2.00

last modified by Zekun Yin 2020/5/18
*/
#ifndef H_FASTQSTREAM
#define H_FASTQSTREAM

#include "Globals.h"

#include "Buffer.h"
#include "FastxChunk.h"
#include "FileReader.h"
#include "Reference.h"
#include "utils.h"
#include <iostream>
#include <string>
#include <zlib.h>//support gziped files, functional but inefficient

#if defined(_WIN32)
#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable : 4996)// D_SCL_SECURE
#pragma warning(disable : 4244)// conversion uint64 to uint32
//# pragma warning(disable : 4267)
#define FOPEN fopen
#define FDOPEN fdopen
#define FSEEK _fseeki64
#define FTELL _ftelli64
#define FCLOSE fclose
#elif __APPLE__// Apple by default suport 64 bit file operations (Darwin 10.5+)
#define FOPEN fopen
#define FDOPEN fdopen
#define FSEEK fseek
#define FTELL ftell
#define FCLOSE fclose
#else
#if !defined(_LARGEFILE_SOURCE)
#define _LARGEFILE_SOURCE
#if !defined(_LARGEFILE64_SOURCE)
#define _LARGEFILE64_SOURCE
#endif
#endif
#if defined(_FILE_OFFSET_BITS) && (_FILE_OFFSET_BITS != 64)
#undef _FILE_OFFSET_BITS
#endif
#if !defined(_FILE_OFFSET_BITS)
#define _FILE_OFFSET_BITS 64
#endif
#define FOPEN fopen64
#define FDOPEN fdopen
#define FSEEK fseeko64
#define FTELL ftello64
#define FCLOSE fclose
#endif

namespace rabbit {

    namespace fa {

        /*
         * @brief Fasta file reader
         * @details A class to read fasta (.fa) file
         */
        class FastaFileReader {
        private:
            /// Swap buffer size (16M default)
            static const uint32 SwapBufferSize = 1 << 26;// 16MB
            /// Fasta data data pool
            FastaDataPool &recordsPool;

        public:
            /*
             * @brief FastaFileReader Constructor
             * @param fileName_ Fasta file path
             * @param pool_ Data pool
             * @param halo size
             * @param isZippedNew if true, it will use gzopen to read fileName_
             */
            FastaFileReader(const std::string &fileName_, FastaDataPool &pool_, uint64 halo = 21,
                            bool isZippedNew = false)
                : swapBuffer(SwapBufferSize),
                  bufferSize(0),
                  eof(false),
                  usesCrlf(false),
                  totalSeqs(0),
                  mHalo(halo),
                  isZipped(isZippedNew),
                  recordsPool(pool_) {
                // if(ends_with(fileName_,".gz"))
                if (isZipped) {
                    mZipFile = gzopen(fileName_.c_str(), "r");
                    // isZipped=true;
                    gzrewind(mZipFile);

                } else {
                    mFile = FOPEN(fileName_.c_str(), "rb");
                    if (mFile == NULL) {
                        throw RioException(
                                ("Can not open file to read: " +
                                 fileName_)
                                        .c_str());//--------------need to change----------//
                    }
                }
            }

            /**
             * @brief FastaFileReader Constructor
             * @param fd Fasta file descriptor (if fasta file is opened)
             * @param pool_ Data pool
             * @param halo halo size
             * @param isZippedNew if true, it will use gzopen to read fileName_
             */
            FastaFileReader(int fd, FastaDataPool &pool_, uint64 halo = 21, bool isZippedNew = false)
                : swapBuffer(SwapBufferSize),
                  bufferSize(0),
                  eof(false),
                  usesCrlf(false),
                  totalSeqs(0),
                  mHalo(halo),
                  isZipped(isZippedNew),
                  recordsPool(pool_) {
                if (isZipped) {
                    mZipFile = gzdopen(fd, "r");
                    if (mZipFile == NULL) {
                        throw RioException("Can not open file to read!");//--------------need to change----------//
                    }
                    gzrewind(mZipFile);

                } else {
                    mFile = FDOPEN(fd, "rb");
                    if (mFile == NULL) {
                        throw RioException("Can not open file to read!");//--------------need to change----------//
                    }
                }
            }

            ~FastaFileReader() {
                // std::cerr << "totalSeqs: " << this->totalSeqs << std::endl;
                if (mFile != NULL || mZipFile != NULL) Close();
                // delete mFile;
                // delete mZipFile;
            }

            /// if it is end of file
            bool Eof() const { return eof; }

            FastaChunk *readNextChunk();

            FastaChunk *readNextChunkList();

            bool ReadNextChunk_(FastaChunk *chunk_, SeqInfos &seqInfos);

            bool ReadNextFaChunk_(FastaChunk *chunk_, SeqInfos &seqInfos);

            bool ReadNextFaChunk_(FastaDataChunk *dataChunk_, SeqInfos &seqInfos, bool &continue_read);

            /// close mFile
            void Close() {
                if (mFile != NULL) {
                    FCLOSE(mFile);
                    mFile = NULL;
                }

                if (mZipFile != NULL && isZipped) {
                    gzclose(mZipFile);
                    mZipFile = NULL;
                }
            }

            /**
             * @brief read data from file
             * @param memory_ pointer to store read file
             * @param size_ read size (byte)
             */
            int64 Read(byte *memory_, uint64 size_) {
                if (isZipped) {
                    int64 n = gzread(mZipFile, memory_, size_);
                    if (n == -1) std::cerr << "Error to read gzip file" << std::endl;
                    return n;
                } else {
                    int64 n = fread(memory_, 1, size_, mFile);
                    return n;
                }
            }

        private:
            core::Buffer swapBuffer;
            uint64 bufferSize;
            bool eof;
            bool usesCrlf;
            bool isZipped;

            FILE *mFile = NULL;
            gzFile mZipFile = NULL;

            uint64 mHalo;

        public:
            uint64 totalSeqs = 0;
            uint64 gid = 0;

            // added form fastareader
            SeqInfos seqInfos;
            uint32 numParts;

        private:
            uint64 FindCutPos_(FastaChunk *dataChunk_, uchar *data_, const uint64 size_, const uint64 halo_,
                               SeqInfos &seqInfos);

            /**
             * @brief Skip to end of line
             * @param data_ base pointer
             * @param pos_ start position to skip
             * @param size_ data_ size
             */
            void FastaSkipToEol(uchar *data_, uint64 &pos_, const uint64 size_) {
                // cout << "SkipToEol " << pos_ << " " << size_ << endl;
                ASSERT(pos_ <= size_);

                while (data_[pos_] != '\n' && data_[pos_] != '\r' && pos_ < size_) ++pos_;

                if (data_[pos_] == '\r' && pos_ < size_) {
                    if (data_[pos_ + 1] == '\n') {
                        usesCrlf = true;
                        ++pos_;
                    }
                }
            }

            void SkipToEol(uchar *data_, uint64 &pos_, const uint64 size_) {
                // cout << "SkipToEol " << pos_ << " " << size_ << endl;
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
            void SkipToSol(uchar *data_, uint64 &pos_, const uint64 size_) {
                ASSERT(pos_ < size_);
                if (data_[pos_] == '\n') {
                    --pos_;
                }
                if (data_[pos_] == '\r') {
                    usesCrlf = true;
                    pos_--;
                }
                //find EOL
                while (data_[pos_] != '\n' && data_[pos_] != '\r') {
                    --pos_;
                }
                if (data_[pos_] == '\n') {
                    --pos_;
                }
                if (data_[pos_] == '\r') {
                    usesCrlf = true;
                    pos_--;
                }
            }

            /**
             * @brief If contains '\n' in uchar buffer
             * @param data_ base pointer
             * @param pos_ start position to skip
             * @param size_ data_ size
             */
            bool FindEol(uchar *data_, uint64 &pos_, const uint64 size_) {
                bool found = false;
                uint64 pos0 = pos_;
                // TODO follow SkipToEol for both \n and \n\r
                while (pos0 < size_) {
                    if (data_[pos0] == '\n') {
                        pos_ = pos0;
                        found = true;
                        break;
                    } else {
                        pos0++;
                    }
                }
                return found;
            }
        };

    }// namespace fa

    namespace fq {
        class FastqFileReader {
        private:
            uint32 SwapBufferSize = 1 << 22;
            uint32 GetNxtBuffSize = 1 << 20;

        public:
            /**
             * @brief FastaFileReader Constructor
             * @param fileName_ Fastq file name
             * @param pool_ Data pool
             * @param fileName2_ the second file name if source file is pair-end sequence
             * @param isZippedNew if true, it will use gzopen to read fileName_ and fileName2_
             */
            FastqFileReader(const std::string &fileName_, FastqDataPool &pool_,
                            std::string fileName2_ = "", bool isZippedNew = false, uint32 mxLen_ = 1 << 20)
                : swapBuffer(SwapBufferSize),
                  swapBuffer2(SwapBufferSize),
                  bufferSize(0),
                  bufferSize2(0),
                  eof(false),
                  usesCrlf(false),
                  isZipped(isZippedNew),
                  numParts(0),
                  recordsPool(pool_) {
                GetNxtBuffSize = mxLen_;
                mFqReader = new FileReader(fileName_, isZipped);
                if (fileName2_ != "") {
                    mFqReader2 = new FileReader(fileName2_, isZipped);
                }
            }

            /**
             * @brief FastaFileReader Constructor
             * @param fileName_ Fastq file descriptor
             * @param pool_ Data pool
             * @param fileName2_ the second file descriptor if source file is pair-end sequence
             * @param isZippedNew if true, it will use gzopen to read fd and fd2
             */
            FastqFileReader(int fd, FastqDataPool &pool_, int fd2 = -1, bool isZippedNew = false)
                : swapBuffer(SwapBufferSize),
                  swapBuffer2(SwapBufferSize),
                  bufferSize(0),
                  bufferSize2(0),
                  eof(false),
                  usesCrlf(false),
                  isZipped(isZippedNew),
                  numParts(0),
                  recordsPool(pool_) {
                if (isZippedNew) {
                } else {
                    mFile = FDOPEN(fd, "rb");
                    if (fd != -1) {
                        mFile2 = FDOPEN(fd2, "rb");
                        if (mFile2 == NULL)
                            throw RioException("Can not open file to read: ");
                    }
                    if (mFile == NULL) {
                        throw RioException("Can not open file to read: ");
                    }
                }
            }

            ~FastqFileReader() {
                if (mFqReader != NULL) delete mFqReader;
                if (mFqReader2 != NULL) delete mFqReader2;
            }

            // added from fastxIO.h
            FastqDataChunk *readNextChunk();

            FastqDataChunk *readNextChunk(moodycamel::ReaderWriterQueue<std::pair<char *, int>> *q, atomic_int *d, pair<char *, int> &l);

            void readChunk();

            bool ReadNextChunk_(FastqDataChunk *chunk_);

            bool ReadNextChunk_(FastqDataChunk *chunk_, moodycamel::ReaderWriterQueue<std::pair<char *, int>> *q, atomic_int *d, pair<char *, int> &l);

            FastqDataPairChunk *readNextPairChunk();

            FastqDataPairChunk *readNextPairChunkParallel();

            FastqDataPairChunk *readNextPairChunkParallel(moodycamel::ReaderWriterQueue<std::pair<char *, int>> *q1,
                                                          moodycamel::ReaderWriterQueue<std::pair<char *, int>> *q2,
                                                          atomic_int *d1, atomic_int *d2,
                                                          pair<char *, int> &last1, pair<char *, int> &last2);


            FastqDataPairChunk *readNextPairChunk(moodycamel::ReaderWriterQueue<std::pair<char *, int>> *q1,
                                                  moodycamel::ReaderWriterQueue<std::pair<char *, int>> *q2,
                                                  atomic_int *d1, atomic_int *d2,
                                                  pair<char *, int> &last1, pair<char *, int> &last2);

            FastqDataChunk *readNextPairChunkInterleaved();


            bool ReadNextPairedChunk_(FastqDataChunk *chunk_);

            void Close() {
                if (mFile != NULL) {
                    FCLOSE(mFile);
                    mFile = NULL;
                }
                if (isZipped) {
                    if (mZipFile != NULL) {
                        gzclose(mZipFile);
                        mZipFile = NULL;
                    }
                }
            }

        private:
            core::Buffer swapBuffer;
            uint64 bufferSize;
            // just for pair-end usage
            core::Buffer swapBuffer2;
            uint64 bufferSize2;
            FILE *mFile2 = NULL;

            bool eof;
            bool usesCrlf;
            bool isZipped;

            FILE *mFile = NULL;
            gzFile mZipFile = NULL;

            FileReader *mFqReader = NULL;
            FileReader *mFqReader2 = NULL;

            // added from fastxIO.h
            FastqDataPool &recordsPool;
            uint32 numParts;

            uint64 lastOneReadPos;
            uint64 lastTwoReadPos;

            uint64 GetNextRecordPos_(uchar *data_, uint64 pos_, const uint64 size_);

            uint64 GetPreviousRecordPos_(uchar *data_, uint64 pos_, const uint64 size_);

            void SkipToSol(uchar *data_, uint64 &pos_, const uint64 size_);

            void SkipToEol(uchar *data_, uint64 &pos_, const uint64 size_);
        };

    }// namespace fq

}// namespace rabbit

#endif// H_FASTQSTREAM
