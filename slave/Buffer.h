/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc

  Authors: Lucas Roguski and Sebastian Deorowicz

  Version: 2.00
  Modified by Zekun Yin
  zekun.yin@mail.sdu.edu.cn
*/

#ifndef BUFFER_H
#define BUFFER_H

#include <algorithm>
#include "Globals.h"
//#include "utils.h"

namespace rabbit {

    namespace core {

#define USE_64BIT_MEMORY 1

/**
 * @brief: buffer to store chunk data
 */
        class Buffer {
        public:
            Buffer(uint64 size_) {
            }

            /// return the pointer of buffer
            byte *Pointer() const { return (byte *) buffer; }

#if (USE_64BIT_MEMORY)

            uint64 *Pointer64() const { return buffer; }

#endif



        private:
            Buffer(const Buffer &) {}

            Buffer &operator=(const Buffer &) { return *this; }

#if (USE_64BIT_MEMORY)
            uint64 *buffer;
#else
            byte *buffer;
#endif
            uint64 size;
        };

/**
 * @brief: rabbitio chunk data wapper
 */
        struct DataChunk {
            /// default swap buffer size
            static const uint64 DefaultBufferSize = 1 << 20;  // 1 << 22

            /// chunk data
            Buffer data;
            /// chunk size
            uint64 size;
            /// list to matain all sequence chunk in one part
            DataChunk *next = NULL;  //[haoz:]added for all seqchunk in one part

            DataChunk(const uint64 bufferSize_ = DefaultBufferSize) : data(bufferSize_), size(0) {}

            void Reset() {
                size = 0;
                next = NULL;
            }
        };

        struct DataPairChunk {
            DataChunk *left_part;
            DataChunk *right_part;
        };

    }  // namespace core

}  // namespace rabbit

#endif  // BUFFER_H
