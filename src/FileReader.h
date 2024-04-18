#include "Buffer.h"
#include "FastxChunk.h"
#include "Reference.h"
#include "utils.h"
#include <cstring>
#include <iostream>
#include <string>
#include <cstdio>
#include <vector>

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


#include <fstream>
#include "globalMutex.h"


namespace rabbit {

    class FileReader {


    public:
        FileReader(const std::string &fileName_, bool isZipped, int64 startPos = 0, int64 endPos = 0, bool inMem = 0);


        FileReader(int fd, bool isZipped = false);


        void DecompressMore();


        int64 Read(byte *memory_, uint64 size_);


        int64 ReadSeek(byte *memory_, uint64 size_, uint64 pos_);


        bool FinishRead();

        bool Eof();

        void setEof();

        ~FileReader();

    private:
        static const uint32 IGZIP_IN_BUF_SIZE = 1 << 22;// 4M gziped file onece fetch
        static const uint32 GZIP_HEADER_BYTES_REQ = 1 << 16;
        byte* MemData = nullptr;
        int64 MemDataTotSize = 0;
        int64 MemDataNowPos = 0;
        bool MemDataReadFinish = false;
        bool read_in_mem = 0;

        char* in_buffer[64];
        char* out_buffer[64];
        char* to_read_buffer;
        std::vector<int> block_sizes;
        size_t buffer_tot_size;
        size_t buffer_now_pos;
        FILE *input_file;
        int now_block;


        FILE *mFile = NULL;
        gzFile mZipFile = NULL;
        //igzip usage
        unsigned char *mIgInbuf = NULL;
        bool isZipped = false;
        bool eof = false;
        std::ifstream iff_idx;
        bool iff_idx_end;
        int start_line;
        int end_line;
        int64 start_pos;
        int64 end_pos;
        int64 total_read;
        int64 has_read;
        int read_times;
    };

}//namespace rabbit
