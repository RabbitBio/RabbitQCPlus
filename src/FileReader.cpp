#include "FileReader.h"

extern "C" {
#include <athread.h>
#include <pthread.h>
void slave_decompressfunc();
}
#define USE_LIBDEFLATE
namespace rabbit {

    FileReader::FileReader(const std::string &fileName_, bool isZipped, int64 startPos, int64 endPos, bool inMem) {
        read_in_mem = inMem;
        read_times = 0;
        start_pos = 0;
        end_pos = 0;
        if (ends_with(fileName_, ".gz") || isZipped) {
#ifdef USE_LIBDEFLATE
            start_line = startPos;
            end_line = endPos;
            std::string index_file_name = fileName_ + ".swidx";
            int block_size;
            iff_idx.open(index_file_name);
            if (!iff_idx.is_open()) {
                fprintf(stderr, "%s do not exits\n", index_file_name.c_str());
                exit(0);
            }
            while (iff_idx >> block_size) {
                block_sizes.push_back(block_size);
            }
            if(start_line == end_line) {
                end_line = block_sizes.size();
            }
            for(int i = 0; i < start_line; i++) start_pos += block_sizes[i];
            for(int i = 0; i < end_line; i++) end_pos += block_sizes[i];
            fprintf(stderr, "FileReader zip [%lld %lld] [%d %d]\n", start_pos, end_pos, start_line, end_line);
            now_block = start_line;
            block_sizes.resize(end_line);
            iff_idx_end = 0;
            for (int i = 0; i < 64; i++) {
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
            int64 file_offset = 0;
            for (int i = 0; i < start_line; i++) {
                file_offset += block_sizes[i];
            }
            fseek(input_file, file_offset, SEEK_SET);
            if (read_in_mem) {
                MemDataTotSize = end_pos - start_pos;
                fprintf(stderr, "read_in_mem %s %lld\n", fileName_.c_str(), MemDataTotSize);
                MemData = new byte[MemDataTotSize];
                int64 n = fread(MemData, 1, MemDataTotSize, input_file);
                fprintf(stderr, "read %lld bytes\n", n);
                if (n != MemDataTotSize) {
                    fprintf(stderr, "n != MemDataTotSize\n");
                    exit(0);
                }
            }

#else
            mZipFile = gzopen(fileName_.c_str(), "r");
            gzrewind(mZipFile);
            if (read_in_mem) {
                fprintf(stderr, "zlib not support read in mem\n");
                exit(0);
            }
#endif
            this->isZipped = true;
        } else {
            start_pos = startPos;
            end_pos = endPos;
            if(start_pos == end_pos) {
                std::ifstream gFile;
                gFile.open(fileName_.c_str());
                gFile.seekg(0, std::ios_base::end);
                end_pos = gFile.tellg();
                gFile.close();
            }
            fprintf(stderr, "FileReader [%lld %lld]\n", start_pos, end_pos);
            has_read = 0;
            total_read = end_pos - start_pos;
            mFile = FOPEN(fileName_.c_str(), "rb");
            if (mFile == NULL) {
                throw RioException(("Can not open file to read: " + fileName_).c_str());
            }
            fseek(mFile, start_pos, SEEK_SET);
            if (read_in_mem) {
                fprintf(stderr, "read_in_mem %s\n", fileName_.c_str());
                MemDataTotSize = end_pos - start_pos;
                MemData = new byte[MemDataTotSize];
                int64 n = fread(MemData, 1, MemDataTotSize, mFile);
                fprintf(stderr, "read %lld bytes\n", n);
                if (n != MemDataTotSize) {
                    fprintf(stderr, "n != MemDataTotSize\n");
                    exit(0);
                }
            }

        }
        //fprintf(stderr, "filereader init done\n");
    }

    FileReader::FileReader(int fd, bool isZipped) {
        fprintf(stderr, "FileReader fd is todo ...\n");
        exit(0);
    }

    FileReader::~FileReader() {
        if (mIgInbuf != NULL) delete mIgInbuf;
        if (read_in_mem) delete[] MemData;
        if (mFile != NULL) {
            fclose(mFile);
            mFile = NULL;
        }
        if (mZipFile != NULL) {
            gzclose(mZipFile);
        }
#ifdef USE_LIBDEFLATE
        if (isZipped) {
            if (input_file != NULL) {
                fclose(input_file);
                input_file = NULL;
            }
            if (iff_idx.is_open()) {
                iff_idx.close();
            }
            for (int i = 0; i < 64; i++) {
                delete[] in_buffer[i];
                delete[] out_buffer[i];
            }
            delete[] to_read_buffer;
        }
#endif
    }

    void FileReader::DecompressMore() {
        size_t length_to_copy = buffer_tot_size - buffer_now_pos;
        if (length_to_copy > buffer_now_pos) {
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
        for (int i = 0; i < 64; i++) {
            if (now_block == end_line) {
                to_read = 0;
                iff_idx_end = 1;
            } else {
                to_read = block_sizes[now_block++];
            }
//            fprintf(stderr, "to_read %zu, now_block %d\n", to_read, now_block);
            if(read_in_mem) {
                int64 lastDataSize = MemDataTotSize - MemDataNowPos;
                if (to_read > lastDataSize) {
                    memcpy(in_buffer[i], MemData + MemDataNowPos, lastDataSize);
                    MemDataNowPos += lastDataSize;
                    MemDataReadFinish = 1;
                    in_size = lastDataSize;
                } else {
                    memcpy(in_buffer[i], MemData + MemDataNowPos, to_read);
                    MemDataNowPos += to_read;
                    in_size = to_read;
                }
            } else {
                in_size = fread(in_buffer[i], 1, to_read, input_file);
            }
            if (in_size == 0) ok = 0;
            paras[i].in_buffer = in_buffer[i];
            paras[i].out_buffer = out_buffer[i];
            paras[i].in_size = in_size;
            paras[i].out_size = &(out_size[i]);
        }

        {
            std::lock_guard <std::mutex> guard(globalMutex);
            __real_athread_spawn((void *) slave_decompressfunc, paras, 1);
            athread_join();
        }

        for (int i = 0; i < 64; i++) {
            if (out_size[i] <= 0) continue;
            memcpy(to_read_buffer + buffer_tot_size, paras[i].out_buffer, out_size[i]);
            buffer_tot_size += out_size[i];
            if (buffer_tot_size > (64 + 8) * BLOCK_SIZE) {
                fprintf(stderr, "buffer_tot_size is out of size\n");
                exit(0);
            }
        }
    }


    int64 FileReader::Read(byte *memory_, uint64 size_) {
        read_times++;
//        if(read_in_mem) {
//            if(isZipped) {
//                if(read_times > 64) usleep(6000);
//            }
//            else {
//                if(read_times > 64) usleep(10000);
//            }
//        }

        if (isZipped) {
#ifdef USE_LIBDEFLATE
            if (buffer_now_pos + size_ > buffer_tot_size) {
                DecompressMore();
                if (buffer_now_pos + size_ > buffer_tot_size) {
                    fprintf(stderr, "after decom still not big enough, read done, -- %d\n", iff_idx_end);
                    size_ = buffer_tot_size - buffer_now_pos;
                }
            }
            memcpy(memory_, to_read_buffer + buffer_now_pos, size_);
            buffer_now_pos += size_;
            return size_;
#else
            if(read_in_mem) {
                fprintf(stderr, "zlib not support read in mem\n");
                exit(0);
            }
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
            if (read_in_mem) {
                int64 lastDataSize = MemDataTotSize - MemDataNowPos;
                if (size_ > lastDataSize) {
                    memcpy(memory_, MemData + MemDataNowPos, lastDataSize);
                    MemDataNowPos += lastDataSize;
                    MemDataReadFinish = 1;
                    return lastDataSize;
                } else {
                    memcpy(memory_, MemData + MemDataNowPos, size_);
                    MemDataNowPos += size_;
                    return size_;
                }
            } else {
                if(has_read + size_ > total_read) size_ = total_read - has_read;
                int64 n = fread(memory_, 1, size_, mFile);
                //fprintf(stderr, "normal read %lld\n", n);
                has_read += n;
                return n;
            }
        }
    }

    bool FileReader::FinishRead() {
        if (isZipped) {
#ifdef USE_LIBDEFLATE
            if(read_in_mem) return (iff_idx_end || MemDataReadFinish) && (buffer_now_pos == buffer_tot_size);
            else return (iff_idx_end || feof(input_file)) && (buffer_now_pos == buffer_tot_size);
#else
            return gzeof(mZipFile);
#endif
        } else {
            //fprintf(stderr, " finish ? %lld %lld\n", total_read, has_read);
            if (read_in_mem) return MemDataReadFinish;
            else return total_read == has_read;
        }
    }

    bool FileReader::Eof() {
        return eof;
    }

    void FileReader::setEof() {
        eof = true;
    }

}
