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
#pragma warning(disable : 4996)  // D_SCL_SECURE
#pragma warning(disable : 4244)  // conversion uint64 to uint32
#pragma warning(disable : 4267)
#pragma warning(disable : 4800)  // conversion byte to bool
#endif

// assertions
//
#if defined(DEBUG) || defined(_DEBUG)

#include "assert.h"

#define ASSERT(x) assert(x)
#else
#define ASSERT(x) (void)(x)
#endif
#define UMI_LOC_NONE 0
#define UMI_LOC_INDEX1 1
#define UMI_LOC_INDEX2 2
#define UMI_LOC_READ1 3
#define UMI_LOC_READ2 4
#define UMI_LOC_PER_INDEX 5
#define UMI_LOC_PER_READ 6

#include <string>
#include <stdexcept>
#include <sys/time.h>


inline std::string replace(const std::string &str, const std::string &src, const std::string &dest) {
    std::string ret;

    std::string::size_type pos_begin = 0;
    std::string::size_type pos = str.find(src);
    while (pos != std::string::npos) {
        ret.append(str.data() + pos_begin, pos - pos_begin);
        ret += dest;
        pos_begin = pos + 1;
        pos = str.find(src, pos_begin);
    }
    if (pos_begin < str.length()) {
        ret.append(str.begin() + pos_begin, str.end());
    }
    return ret;
}

inline double GetTime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_sec + (double) tv.tv_usec / 1000000;
}

static int valAGCT[8] = {-1, 0, -1, 2, 1, -1, -1, 3};
static int valAGCT2[8] = {-1, 0, -1, 2, 1, -1, 4, 3};


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

        const char *what() const throw()  // for std::exception interface
        {
            return message.c_str();
        }
    };

}  // namespace rabbit

#endif  // H_GLOBALS
