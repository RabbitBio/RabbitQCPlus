#include "Buffer.h"
#include "FastxChunk.h"
#include "Reference.h"
#include "utils.h"
#include <cstring>
#include <iostream>
#include <string>
#include <cstdio>

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


extern "C" {
#include <athread.h>
#include <pthread.h>
    void slave_decompressfunc();
}

#include <fstream>
#include "globalMutex.h"



namespace rabbit {

    class FileReader {
    private:
        static const uint32 IGZIP_IN_BUF_SIZE = 1 << 22;// 4M gziped file onece fetch
        static const uint32 GZIP_HEADER_BYTES_REQ = 1 << 16;

    public:
        FileReader(const std::string &fileName_, bool isZipped) {
            if (ends_with(fileName_, ".gz") || isZipped) {
#ifdef USE_LIBDEFLATE
                fprintf(stderr, "use libdeflate_slave\n");
                std::string index_file_name = fileName_ + ".swidx";
                fprintf(stderr, "index_file_name %s\n", index_file_name.c_str());
                iff_idx.open(index_file_name);
                if(!iff_idx.is_open()) {
                    fprintf(stderr, "%s do not exits\n", index_file_name.c_str());
                    exit(0);
                }
                for(int i = 0; i < 64; i++) {
                    in_buffer[i] = new char[BLOCK_SIZE];
                    out_buffer[i] = new char[BLOCK_SIZE];
                }
                to_read_buffer = new char[(64 + 8) * BLOCK_SIZE];
                buffer_tot_size = 0;
                buffer_now_pos = 0;
                input_file = fopen(fileName_.c_str(), "rb");
                if (input_file == NULL) {
                    perror("Failed to open input file");
                    exit(0);
                }
#else 
                mZipFile = gzopen(fileName_.c_str(), "r");
                gzrewind(mZipFile);
#endif
                this->isZipped = true;
            } else {
                mFile = FOPEN(fileName_.c_str(), "rb");
                if (fileName_ != "") {
                    mFile = FOPEN(fileName_.c_str(), "rb");
                    if (mFile == NULL)
                        throw RioException("Can not open file to read: ");
                }
                if (mFile == NULL) {
                    throw RioException(
                            ("Can not open file to read: " + fileName_).c_str());
                }
            }
            fprintf(stderr, "filereader init done\n");
        }

        FileReader(int fd, bool isZipped = false) {
            if (isZipped) {
                fprintf(stderr, "FileReader fd is todo ...\n");
                exit(0);
            } else {
                mFile = FDOPEN(fd, "rb");
                if (fd != -1) {
                    mFile = FDOPEN(fd, "rb");
                    if (mFile == NULL)
                        throw RioException("Can not open file to read: ");
                }
                if (mFile == NULL) {
                    throw RioException("Can not open file to read: ");
                }
            }
        }

        void DecompressMore() {
            size_t length_to_copy = buffer_tot_size - buffer_now_pos;
            if(length_to_copy > buffer_now_pos) {
                fprintf(stderr, "warning: memmove may GG\n");
            }
            std::memmove(to_read_buffer, to_read_buffer + buffer_now_pos, length_to_copy);
            buffer_tot_size = length_to_copy;
            buffer_now_pos = 0;
        
            int ok = 1;
            Para paras[64];
            size_t to_read = 0;
            size_t in_size = 0;
            size_t out_size[64] = {0};
            for(int i = 0; i < 64; i++) {
                iff_idx >> to_read;
                if(iff_idx.eof()) ok = 0;
                in_size = fread(in_buffer[i], 1, to_read, input_file);
                if (in_size == 0) ok = 0;
                //fprintf(stderr, "-- %d %d\n", to_read, in_size);
                paras[i].in_buffer = in_buffer[i];
                paras[i].out_buffer = out_buffer[i];
                paras[i].in_size = in_size;
                paras[i].out_size = &(out_size[i]);
            }
            //if(ok == 0) fprintf(stderr, "%d %d\n", iff_idx.eof(), feof(input_file));

            {
                std::lock_guard<std::mutex> guard(globalMutex);
                __real_athread_spawn((void *)slave_decompressfunc, paras, 1);
                athread_join();
            }

            //for(int i = 0; i < 64; i++) {
            //    if(in_size == 0) continue;
            //    size_t out_size_tmp = -1;
            //    libdeflate_decompressor* decompressor = libdeflate_alloc_decompressor();
            //    libdeflate_result result = libdeflate_gzip_decompress(decompressor, paras[i].in_buffer, paras[i].in_size, paras[i].out_buffer, BLOCK_SIZE, &out_size_tmp);
            //    if (result != LIBDEFLATE_SUCCESS) {
            //        fprintf(stderr, "Decompression failed\n");
            //    }
            //    *(paras[i].out_size) = out_size_tmp;
            //    libdeflate_free_decompressor(decompressor);
            //}
            for(int i = 0; i < 64; i++) {
                if (out_size[i] <= 0) continue;
                //fprintf(stderr, "====111 %d %d\n", out_size[i], buffer_tot_size);
                memcpy(to_read_buffer + buffer_tot_size, paras[i].out_buffer, out_size[i]);
                //fprintf(stderr, "====222 %d\n", out_size[i]);
                buffer_tot_size += out_size[i];
                if(buffer_tot_size > (64 + 8) * BLOCK_SIZE) {
                    fprintf(stderr, "buffer_tot_size is out of size\n");
                    exit(0);
                }
            }
        }


        int64 Read(byte *memory_, uint64 size_) {
            if (isZipped) {
#ifdef USE_LIBDEFLATE
                if (buffer_now_pos + size_ > buffer_tot_size) {
                    DecompressMore();
                    //fprintf(stderr, "decom more %lld %lld\n", buffer_now_pos, buffer_tot_size);
                    if(buffer_now_pos + size_ > buffer_tot_size) {
                        fprintf(stderr, "after decom still not big enough, read done\n");
                        size_ = buffer_tot_size - buffer_now_pos;
                    }
                }
                memcpy(memory_, to_read_buffer + buffer_now_pos, size_);
                buffer_now_pos += size_;
                return size_;
#else
                int64 n = gzread(mZipFile, memory_, size_);
                if (n == -1) {
                    int errNum;
                    const char* errorMsg = gzerror(mZipFile, &errNum);
                    if (errNum) {
                        std::cerr << "Error to read gzip file: " << errorMsg << std::endl;
                    }
                }
                return n;
#endif
            } else {
                int64 n = fread(memory_, 1, size_, mFile);
                return n;
            }
        }


        int64 ReadSeek(byte *memory_, uint64 size_, uint64 pos_) {
            if (isZipped) {
                if(pos_ == 0) {
                    return Read(memory_, size_);
                } else {
                    //TODO
                    fprintf(stderr, "TODO seek read of compressed file...\n");
                    exit(0);
                }
           } else {
                fseek(mFile, pos_, SEEK_SET);
                int64 n = fread(memory_, 1, size_, mFile);
                return n;
            }
        }

        /**
         *
         * @param buf
         * @param len
         * @param Q
         * @param done
         * @param L
         * @param num
         * @return
         */
        int64 Read(byte *buf, uint64 len, moodycamel::ReaderWriterQueue<std::pair<char *, int>>

                                                  *Q,
                   std::atomic_int *done, std::pair<char *, int> &L) {
            std::pair<char *, int> now;
            int64 ret;
            int64 got = 0;
            //                printf("now producer get data from pugz queue\n");
            if (L.second > 0) {
                if (L.second >= len) {
                    memcpy(buf, L.first, len);
                    char *tmp = new char[L.second - len];
                    memcpy(tmp, L.first + len, L.second - len);
                    memcpy(L.first, tmp, L.second - len);
                    delete[] tmp;
                    L.second = L.second - len;
                    ret = len;
                    buf += ret;
                    len -= ret;
                    got += ret;
                    return got;
                } else {
                    memcpy(buf, L.first, L.second);
                    ret = L.second;
                    L.second = 0;
                    buf += ret;
                    len -= ret;
                    got += ret;
                }
            }
            bool overWhile = false;
            while (len > 0) {
                while (Q->try_dequeue(now) == 0) {
                    if (Q->size_approx() == 0 && *done == 1) {
                        ret = 0;
                        overWhile = true;
                        break;
                    }
                    usleep(100);
                    //printf("producer waiting pugz...\n");
                }
                if (overWhile) {
                    ret = 0;
                    break;
                }
                //printf("get some data %d\n", now.second);
                if (now.second <= len) {
                    memcpy(buf, now.first, now.second);
                    delete[] now.first;
                    ret = now.second;
                } else {
                    int move_last = now.second - len;
                    memcpy(buf, now.first, len);
                    memcpy(L.first, now.first + len, move_last);
                    L.second = move_last;
                    delete[] now.first;
                    ret = len;
                }

                buf += ret;
                len -= ret;
                got += ret;
            }

            return got;
        }

        /// True means no need to call Read() function
        bool FinishRead() {
            if (isZipped) {
#ifdef USE_LIBDEFLATE
                //if(iff_idx.eof() != feof(input_file)) {
                //    fprintf(stderr, "iff_idx.eof() != feof(input_file)\n");
                //    exit(0);
                //}
                return iff_idx.eof() || feof(input_file);
#else
                return gzeof(mZipFile);
#endif
            } else {
                return feof(mFile);
            }
        }

        bool Eof() const {
            return eof;
            //if(eof) return eof;
            //return feof(mFile);
        }

        void setEof() {
            eof = true;
        }

        ~FileReader() {
            if (mIgInbuf != NULL) delete mIgInbuf;
            if (mFile != NULL) {
                fclose(mFile);
                mFile = NULL;
            }
            if (mZipFile != NULL) {
                gzclose(mZipFile);
            }
#ifdef USE_LIBDEFLATE
            if(isZipped) {
                if (input_file != NULL) {
                    fclose(input_file);
                    input_file = NULL;
                }
                if (iff_idx.is_open()) {
                    iff_idx.close();
                }
                for(int i = 0; i < 64; i++) {
                    delete[] in_buffer[i];
                    delete[] out_buffer[i];
                }
                delete[] to_read_buffer;
            }
#endif
        }

    private:
        char* in_buffer[64];
        char* out_buffer[64];
        char* to_read_buffer;
        size_t buffer_tot_size;
        size_t buffer_now_pos;
        FILE *input_file;


        FILE *mFile = NULL;
        gzFile mZipFile = NULL;
        //igzip usage
        unsigned char *mIgInbuf = NULL;
        bool isZipped = false;
        bool eof = false;
        std::ifstream iff_idx;
    };

}//namespace rabbit
