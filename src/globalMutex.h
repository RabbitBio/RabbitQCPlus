#ifndef GLOBALMUTEX_H
#define GLOBALMUTEX_H

#define USE_LIBDEFLATE

#include <mutex>

#define BLOCK_SIZE (1 << 22) 
extern std::mutex globalMutex;
struct Para {
    char* in_buffer;
    char* out_buffer;
    size_t *out_size;
    size_t in_size;
    int level;
};



#endif
