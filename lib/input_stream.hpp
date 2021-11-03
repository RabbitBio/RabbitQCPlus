#ifndef INPUT_STREAM_HPP
#define INPUT_STREAM_HPP

#include <cstring> // memcpy
#include <type_traits>

#include "memory.hpp"
#include "assert.hpp"
#include "gzip_constants.h"

/**
 * @brief Model an compressed gzip input stream
 * It can be read by dequeing n<32 bits at a time or as byte aligned uint16_t words

 * The state of the "input bitstream" consists of the following variables:
 *
 *	- in_next: pointer to the next unread byte in the input buffer
 *
 *	- in_end: pointer just past the end of the input buffer
 *
 *	- bitbuf: a word-sized variable containing bits that have been read from
 *		  the input buffer.  The buffered bits are right-aligned
 *		  (they're the low-order bits).
 *
 *	- bitsleft: number of bits in 'bitbuf' that are valid.
 *
 */
class InputStream
{
  public: // protected: //FIXME
    /*
     * The type for the bitbuffer variable ('bitbuf' described above).  For best
     * performance, this should have size equal to a machine word.
     *
     * 64-bit platforms have a significant advantage: they get a bigger bitbuffer
     * which they have to fill less often.
     */
    using bitbuf_t      = machine_word_t;
    using bitbuf_size_t = uint_fast32_t;

    /** State */
    const byte* restrict in_next; /// Read pointer
    span<const byte>     data;
    bitbuf_t             bitbuf        = bitbuf_t(0); /// Bit buffer
    bitbuf_size_t        bitsleft      = 0;           /// Number of valid bits in the bit buffer
    bitbuf_size_t        overrun_count = 0;

    /**
     * Number of bits the bitbuffer variable can hold.
     */
    static constexpr bitbuf_size_t bitbuf_length = 8 * sizeof(bitbuf_t);

    /**
     * The maximum number of bits that can be requested to be in the bitbuffer
     * variable.  This is the maximum value of 'n' that can be passed
     * ensure_bits(n).
     *
     * This not equal to BITBUF_NBITS because we never read less than one byte at a
     * time.  If the bitbuffer variable contains more than (BITBUF_NBITS - 8) bits,
     * then we can't read another byte without first consuming some bits.  So the
     * maximum count we can ensure is (BITBUF_NBITS - 7).
     */
    static constexpr bitbuf_size_t bitbuf_max_ensure = bitbuf_length - 7;

    /**
     * Does the bitbuffer variable currently contain at least 'n' bits?
     */
    bool have_bits(size_t n) const { return bitsleft >= n; }

    /**
     * Fill the bitbuffer variable by reading the next word from the input buffer.
     * This can be significantly faster than fill_bits_bytewise().  However, for
     * this to work correctly, the word must be interpreted in little-endian format.
     * In addition, the memory access may be unaligned.  Therefore, this method is
     * most efficient on little-endian architectures that support fast unaligned
     * access, such as x86 and x86_64.
     */
    hot_fun void fill_bits_wordwise()
    {
        assert(bitsleft <= bitbuf_length);
        bitbuf_t dst;
        memcpy(&dst, in_next, sizeof(bitbuf_t));
        bitbuf |= dst << bitsleft;
        in_next += (bitbuf_length - bitsleft) >> 3u;
        bitsleft += (bitbuf_length - bitsleft) & ~7u;
        assert(bitsleft <= bitbuf_length);
    }

    /**
     * Fill the bitbuffer variable, reading one byte at a time.
     *
     * Note: if we would overrun the input buffer, we just don't read anything,
     * leaving the bits as 0 but marking them as filled.  This makes the
     * implementation simpler because this removes the need to distinguish between
     * "real" overruns and overruns that occur because of our own lookahead during
     * Huffman decoding.  The disadvantage is that a "real" overrun can go
     * undetected, and libdeflate_deflate_decompress() may return a success status
     * rather than the expected failure status if one occurs.  However, this is
     * irrelevant because even if this specific case were to be handled "correctly",
     * one could easily come up with a different case where the compressed data
     * would be corrupted in such a way that fully retains its validity.  Users
     * should run a checksum against the uncompressed data if they wish to detect
     * corruptions.
     */
    cold_fun noinline_fun void fill_bits_bytewise()
    {
        do {
            if (likely(in_next != data.end()))
                bitbuf |= bitbuf_t(*in_next++) << bitsleft;
            else
                overrun_count++;
            bitsleft += 8;
        } while (bitsleft <= bitbuf_length - 8);
        assert(bitsleft <= bitbuf_length);
    }

  public:
    InputStream(const byte* in, size_t len)
      : in_next(in)
      , data(in, len)
    {}

    InputStream(const InputStream&) = default;

    /// states asignments for the same stream (for backtracking)
    InputStream& operator=(const InputStream& from)
    {
        assert(data == from.data);

        in_next       = from.in_next;
        bitbuf        = from.bitbuf;
        bitsleft      = from.bitsleft;
        overrun_count = from.overrun_count;
        return *this;
    }

    InputStream(InputStream&&) = default;
    InputStream& operator=(InputStream&&) = default;

    bool consume_header()
    {
        align_input();
        if (available() < GZIP_MIN_OVERHEAD) return false;

        const byte* p = in_next;

        /* ID1 */
        if (*p++ != GZIP_ID1) return false;
        /* ID2 */
        if (*p++ != GZIP_ID2) return false;
        /* CM */
        if (*p++ != GZIP_CM_DEFLATE) return false;
        auto flg = static_cast<uint8_t>(*p++);
        /* MTIME */
        p += 4;
        /* XFL */
        p += 1;
        /* OS */
        p += 1;

        if (bool(flg & GZIP_FRESERVED)) return false;

        /* Extra field */
        if (bool(flg & GZIP_FEXTRA)) {
            uint16_t xlen;
            memcpy(&xlen, p, sizeof(uint16_t));
            p += sizeof(uint16_t);

            if (data.end() - p < uint32_t(xlen) + GZIP_FOOTER_SIZE) return false;

            p += xlen;
        }

        /* Original file name (zero terminated) */
        if (bool(flg & GZIP_FNAME)) {
            while (*p++ != byte(0) && p != data.end())
                ;
            if (data.end() - p < GZIP_FOOTER_SIZE) return false;
        }

        /* File comment (zero terminated) */
        if (bool(flg & GZIP_FCOMMENT)) {
            while (*p++ != byte(0) && p != data.end())
                ;
            if (data.end() - p < GZIP_FOOTER_SIZE) return false;
        }

        /* CRC16 for gzip header */
        if (bool(flg & GZIP_FHCRC)) {
            p += 2;
            if (data.end() - p < GZIP_FOOTER_SIZE) return false;
        }

        in_next = p;
        return true;
    }

    ssize_t consume_footer()
    {
        align_input();
        if (reinterpret_cast<uintptr_t>(in_next) % alignof(uint32_t) == 0 || available() < GZIP_FOOTER_SIZE) return -1;

        static_assert(sizeof(uint32_t[2]) == GZIP_FOOTER_SIZE, "Gzip footer size doesn't match sizeof(uint32_t[2])");
        uint32_t footer[2];
        memcpy(footer, in_next, GZIP_FOOTER_SIZE);
        in_next += GZIP_FOOTER_SIZE;

        return footer[1]; // decompressed size modulo 4GB
    }

    size_t size() const { return data.size(); }

    /** Remaining available bytes
     * @note align_input() should be called first in order to get accurate readings (or use available_bits() / 8)
     */
    size_t available() const restrict
    {
        assert(data.includes(in_next) || in_next == data.end());
        return size_t(data.end() - in_next);
    }

    /// Remaining available bits
    size_t available_bits() const { return 8 * available() + bitsleft; }

    /** Position in the stream in bits
     */
    size_t position_bits() const { return 8 * size_t(in_next - data.begin()) - bitsleft; }

    bool set_position_bits(size_t bit_pos)
    {
        size_t pos_bytes = bit_pos >> 3u;
        size_t pos_bits  = bit_pos & 7;

        in_next = data.begin() + pos_bytes;
        assert(data.includes(in_next));

        if (unlikely(available() < sizeof(bitbuf_t))) return false;

        memcpy(&bitbuf, in_next, sizeof(bitbuf_t));
        in_next += +sizeof(bitbuf_t);
        bitbuf >>= pos_bits;
        bitsleft = bitbuf_length - pos_bits;

        assert(position_bits() == bit_pos);
        return true;
    }

    /// Seek forward
    void skip(size_t offset)
    {
        align_input();
        in_next += offset;
        assert(data.includes(in_next));
    }

    /**
     * Load more bits from the input buffer until the specified number of bits is
     * present in the bitbuffer variable.  'n' cannot be too large; see MAX_ENSURE
     * and CAN_ENSURE().
     */
    template<bitbuf_size_t n> hot_fun forceinline_fun bool ensure_bits()
    {
        static_assert(n <= bitbuf_max_ensure, "Bit buffer is too small");
        if (have_bits(n)) {
            return true;
        } else {
            if (likely(available() >= sizeof(bitbuf_t))) {
                fill_bits_wordwise();
                return true;
            } else if (available() >= 1) {
                fill_bits_bytewise();
                return true;
            } else {
                return false; // This is not acceptable overrun
            }
        }
    }

    /**
     * Return the next 'n' bits from the bitbuffer variable without removing them.
     */
    template<typename T = uint32_t> T bits(bitbuf_size_t n = 8 * sizeof(T)) const
    {
        assert(bitsleft >= n);
        return T(bitbuf & ((bitbuf_t(1) << n) - 1));
    }

    /**
     * Remove the next 'n' bits from the bitbuffer variable.
     */
    void remove_bits(bitbuf_size_t n)
    {
        assert(bitsleft >= n);
        bitbuf >>= n;
        bitsleft -= n;
    }

    /**
     * Remove and return the next 'n' bits from the bitbuffer variable.
     */
    template<typename T = uint32_t> T pop_bits(bitbuf_size_t n = 8 * sizeof(T))
    {
        assert(n <= 8 * sizeof(T));
        T tmp = bits<T>(n);
        remove_bits(n);
        return tmp;
    }

    /**
     * Align the input to the next byte boundary, discarding any remaining bits in
     * the current byte.
     *
     * Note that if the bitbuffer variable currently contains more than 8 bits, then
     * we must rewind 'in_next', effectively putting those bits back.  Only the bits
     * in what would be the "current" byte if we were reading one byte at a time can
     * be actually discarded.
     */
    void align_input()
    {
        assert(overrun_count <= (bitsleft >> 3));
        in_next -= (bitsleft >> 3) - overrun_count;
        // was:
        // in_next -= (bitsleft >> 3) - std::std::min(overrun_count, bitsleft >> 3);
        bitbuf   = 0;
        bitsleft = 0;
    }

    /**
     * Read a 16-bit value from the input.  This must have been preceded by a call
     * to ALIGN_INPUT(), and the caller must have already checked for overrun.
     */
    uint16_t pop_u16()
    {
        assert(available() >= sizeof(uint16_t));
        uint16_t tmp;
        memcpy(&tmp, in_next, sizeof(uint16_t));
        in_next += 2;
        return tmp;
    }

    /**
     * Copy n bytes to the ouput buffer. The input buffer must be aligned with a
     * call to align_input()
     */
    template<typename char_t> void copy(char_t* restrict out, size_t n) restrict
    {
        // This version support characters representation in output stream wider than bytes
        size_t nbytes = n * sizeof(char_t);
        assert(available() >= n * nbytes);
        memcpy(out, in_next, nbytes);
        in_next += nbytes;
    }

    /**
     * Checks that the lenght next bytes are ascii
     * (for checked copy()ies of uncompressed blocks)
     */
    bool check_ascii(size_t n, uint8_t min_value = '\t', uint8_t max_value = '~')
    {
        if (unlikely(n > available())) return false;

        for (size_t i = 0; i < n; i++) {
            byte c = in_next[i];
            if (c > byte(max_value) || c < byte(min_value)) return false;
        }
        return true;
    }
};

#endif // INPUT_STREAM_HPP
