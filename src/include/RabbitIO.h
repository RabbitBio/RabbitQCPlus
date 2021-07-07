#include <algorithm>
#include "Globals.h"
#include "utils.h"

namespace rabbit {

namespace core {

#define USE_64BIT_MEMORY 1

/**
 * @brief: buffer to store chunk data
 */
class Buffer {
 public:
  Buffer(uint64 size_) {
    ASSERT(size_ != 0);

#if (USE_64BIT_MEMORY)
    uint64 size64 = size_ / 8;
    if (size64 * 8 < size_) size64 += 1;
    buffer = new uint64[size64];
#else
    buffer = new byte[size_];
#endif
    size = size_;
  }

  ~Buffer() { delete buffer; }

  uint64 Size() const { return size; }

	/// return the pointer of buffer
  byte *Pointer() const { return (byte *)buffer; }

#if (USE_64BIT_MEMORY)
  uint64 *Pointer64() const { return buffer; }
#endif
	/// resize the buffer
  void Extend(uint64 size_, bool copy_ = false) {
#if (USE_64BIT_MEMORY)
    uint64 size64 = size / 8;
    if (size64 * 8 < size) size64 += 1;

    uint64 newSize64 = size_ / 8;
    if (newSize64 * 8 < size_) newSize64 += 1;

    if (size > size_) return;

    if (size64 == newSize64) {
      size = size_;
      return;
    }

    uint64 *p = new uint64[newSize64];

    if (copy_) std::copy(buffer, buffer + size64, p);
#else
    if (size > size_) return;

    byte *p = new byte[size_];

    if (copy_) std::copy(buffer, buffer + size, p);

#endif
    delete[] buffer;

    buffer = p;
    size = size_;
  }

  void Swap(Buffer &b) {
    TSwap(b.buffer, buffer);
    TSwap(b.size, size);
  }

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
  DataChunk *next = NULL;  //list structre, for all seqchunk in one part

  DataChunk(const uint64 bufferSize_ = DefaultBufferSize) : data(bufferSize_), size(0) {}
  void Reset() { size = 0; next = NULL;}
};

struct DataPairChunk {
  DataChunk *left_part;
  DataChunk *right_part;
};

//[file:s] <util.h>
template <uint32 _TBitNum>
class TBitMask {
 public:
  static const uint64 Value = ((uint64)1 << (_TBitNum - 1)) | TBitMask<_TBitNum - 1>::Value;
};

template <>
class TBitMask<0> {
 public:
  static const uint64 Value = 0;
};

template <uint32 _TNum>
struct TLog2 {
  static const uint32 Value = TLog2<(_TNum >> 1)>::Value + 1;
};

template <>
struct TLog2<1> {
  static const uint32 Value = 0;
};

template <>
struct TLog2<0> {};

template <typename _T>
inline void TSwap(_T &x_, _T &y_) {
  _T tmp = x_;
  x_ = y_;
  y_ = tmp;
}

template <typename _T>
inline void TFree(_T *p_) {
  if (p_ != (_T *)0) delete p_, p_ = (_T *)0;
}

inline uint32 to_string(uchar *str, uint32 value) {
  uint32 digits;
  uint32 power = 1;

  if (value == 0) {
    str[0] = '0';
    return 1;
  }

  for (digits = 0; digits < 10; ++digits) {
    if (value < power) break;
    power *= 10;
  }

  power /= 10;
  for (uint32 i = 0; power; ++i, power /= 10) {
    int32 d = value / power;
    str[i] = (uchar)('0' + d);
    value -= d * power;
  }

  return digits;
}

inline bool extend_string(uchar *&str, uint32 &size) {
  uint32 new_size = size * 2;
  uchar *p = new uchar[new_size + 1];

  if (!p) return false;

  std::copy(str, str + size, p);
  size = new_size;
  delete[] str;
  str = p;

  return true;
}

inline bool extend_string_to(uchar *&str, uint32 &size, uint32 new_size) {
  if (new_size <= size) return true;

  uchar *p = new uchar[new_size + 1];

  if (!p) return false;

  std::copy(str, str + size, p);

  size = new_size;
  delete[] str;
  str = p;

  return true;
}

inline bool ends_with(const std::string &str_, const std::string &suff_) {
  return str_.size() >= suff_.size() && str_.compare(str_.size() - suff_.size(), suff_.size(), suff_) == 0;
}

inline uint32 int_log(uint32 x, uint32 base) {
  uint32 r = 0;

  if (base == 0) return 1;
  if (base == 1) base++;

  for (uint32 tmp = base; tmp <= x; tmp *= base) ++r;

  return r;
}

inline uint32 to_num(const uchar *str, uint32 len) {
  uint32 r = 0;

  for (uint32 i = 0; i < len; ++i) r = r * 10 + (str[i] - '0');

  return r;
}

inline bool is_num(const uchar *str_, uint32 len_, uint32 &val_) {
  val_ = 0;
  uint32 i;
  for (i = 0; i < len_; ++i) {
    if (str_[i] < '0' || str_[i] > '9') break;
    val_ = val_ * 10 + (str_[i] - '0');
  }

  return i == len_ && (len_ == 1 || str_[0] != '0');
}

inline uint32 bit_length(uint64 x) {
  for (uint32 i = 0; i < 32; ++i) {
    if (x < (1ull << i)) return i;
  }
  return 64;
}
//[file:d] <util.h>

}  // namespace core

}  // namespace rabbit
