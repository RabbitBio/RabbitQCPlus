/*
   This file is a part of DSRC software distributed under GNU GPL 2 licence.
   The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc

Authors: Lucas Roguski and Sebastian Deorowicz

Version: 2.00
*/

#ifndef H_GLOBALS
#define H_GLOBALS

#ifndef NDEBUG
#define DEBUG 1
#endif

// Visual Studio warning supression
//
#if defined(_WIN32)
#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable : 4996)// D_SCL_SECURE
#pragma warning(disable : 4244)// conversion uint64 to uint32
#pragma warning(disable : 4267)
#pragma warning(disable : 4800)// conversion byte to bool
#endif

// assertions
//
#if defined(DEBUG) || defined(_DEBUG)

#include "assert.h"

#define ASSERT(x) assert(x)
#else
#define ASSERT(x) (void) (x)
#endif
#define UMI_LOC_NONE 0
#define UMI_LOC_INDEX1 1
#define UMI_LOC_INDEX2 2
#define UMI_LOC_READ1 3
#define UMI_LOC_READ2 4
#define UMI_LOC_PER_INDEX 5
#define UMI_LOC_PER_READ 6

#include "atomicops.h"
#include "concurrentqueue.h"
#include "readerwriterqueue.h"
#include "util.h"
#include <atomic>
#include <stdexcept>
#include <string>
#include <sys/time.h>

#define slave_num 64

inline double GetTime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_sec + (double) tv.tv_usec / 1000000;
}

static int valAGCT[8] = {-1, 0, -1, 2, 1, -1, -1, 3};
static int valAGCT2[8] = {-1, 1, -1, 3, 2, -1, 5, 4};


inline bool overRepPassed(std::string &seq, int64_t count, int s) {
    switch (seq.length()) {
        case 10:
            return s * count > 500;
        case 20:
            return s * count > 200;
        case 40:
            return s * count > 100;
        case 100:
            return s * count > 50;
        default:
            return s * count > 20;
    }
}


inline std::string PaseFileName(std::string file_name) {
    int n_size = file_name.length();
    int suffix_pos = n_size;
    int prefix_pos = 0;
    for (int i = 0; i < n_size; i++) {
        if (file_name[i] == '/')
            prefix_pos = i + 1;
    }
    int p1 = n_size;
    int p2 = n_size;
    int p3 = n_size;
    if (file_name.find(".fastq") != std::string::npos) p1 = file_name.find(".fastq");
    if (file_name.find(".fq") != std::string::npos) p2 = file_name.find(".fq");
    if (file_name.find(".sra") != std::string::npos) p3 = file_name.find(".sra");
    suffix_pos = std::min(suffix_pos, p1);
    suffix_pos = std::min(suffix_pos, p2);
    suffix_pos = std::min(suffix_pos, p3);
    if (prefix_pos >= suffix_pos) return file_name;
    return file_name.substr(prefix_pos, suffix_pos - prefix_pos);
}


namespace rabbit {

    // basic types
    //
    typedef char int8;
    typedef unsigned char uchar, byte, uint8;
    typedef short int int16;
    typedef unsigned short int uint16;
    typedef int int32;
    typedef unsigned int uint32;
    typedef long long int64;
    typedef unsigned long long uint64;


    static char reMap[123] = {'0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
                              '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
                              '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
                              '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', 'T', 'B', 'G', 'D', 'E', 'F', 'C',
                              'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'A', 'U', 'V', 'W', 'X', 'Y',
                              'Z', '0', '0', '0', '0', '0', '0', 'T', 'b', 'G', 'd', 'e', 'f', 'C', 'h', 'i', 'j', 'k',
                              'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 'A', 'u', 'v', 'w', 'x', 'y', 'z'};

    // exception class
    //
    class RioException : public std::exception {
        std::string message;

    public:
        RioException(const char *msg_) : message(msg_) {}

        RioException(const std::string &msg_) : message(msg_) {}

        ~RioException() throw() {}

        const char *what() const throw()// for std::exception interface
        {
            return message.c_str();
        }
    };

}// namespace rabbit

#endif// H_GLOBALS
