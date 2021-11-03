/*
 * deflate_decompress.c - a decompressor for DEFLATE
 *
 * Originally public domain; changes after 2016-09-07 are copyrighted.
 *
 * Copyright 2016 Eric Biggers
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 *
 * ---------------------------------------------------------------------------
 *
 * This is a highly optimized DEFLATE decompressor.  When compiled with gcc on
 * x86_64, it decompresses data in about 52% of the time of zlib (48% if BMI2
 * instructions are available).  On other architectures it should still be
 * significantly faster than zlib, but the difference may be smaller.
 *
 * Why this is faster than zlib's implementation:
 *
 * - Word accesses rather than byte accesses when reading input
 * - Word accesses rather than byte accesses when copying matches
 * - Faster Huffman decoding combined with various DEFLATE-specific tricks
 * - Larger bitbuffer variable that doesn't need to be filled as often
 * - Other optimizations to remove unnecessary branches
 * - Only full-buffer decompression is supported, so the code doesn't need to
 *   support stopping and resuming decompression.
 * - On x86_64, compile a version of the decompression routine using BMI2
 *   instructions and use it automatically at runtime when supported.
 */

#include <climits>
#include <pmmintrin.h>
#include <tmmintrin.h>

#include <mutex>
#include <condition_variable>
#include <atomic>

#include "common/exceptions.hpp"
#include "memory.hpp"
#include "assert.hpp"

#include "input_stream.hpp"
#include "decompressor.hpp"

#include "libdeflatePugz.h"

#include "immintrin.h"
#include "emmintrin.h"

/** The main class for the deflate decompression
 * Holds the decompressor tables and the input stream, but not the deflate window or output buffers
 */
class DeflateParser {
public:
    DeflateParser(const InputStream &in_stream)
            : _in_stream(in_stream) {}

    enum class block_result : unsigned {
        SUCCESS = 0, // Success, yet many work remaining
        LAST_BLOCK = 1, // Last block had just been decoded
        CAUGHT_UP_DOWNSTREAM = 2, // Caught up downstream decoder
        FLUSH_FAIL = 3, // Not enough space in the buffer
        INVALID_BLOCK_TYPE,
        INVALID_DYNAMIC_HT,
        INVALID_UNCOMPRESSED_BLOCK,
        INVALID_LITERAL,
        INVALID_MATCH,
        TOO_MUCH_INPUT,
        NOT_ENOUGH_INPUT,
        INVALID_PARSE,
    };

    static const char *block_result_to_cstr(block_result result) {
        static constexpr const char *block_result_strings[] = {
                "SUCCESS",
                "LAST_BLOCK",
                "CAUGHT_UP_DOWNSTREAM",
                "FLUSH_FAIL",
                "INVALID_BLOCK_TYPE",
                "INVALID_DYNAMIC_HT",
                "INVALID_UNCOMPRESSED_BLOCK",
                "INVALID_LITERAL",
                "INVALID_MATCH",
                "TOO_MUCH_INPUT",
                "NOT_ENOUGH_INPUT",
                "INVALID_PARSE",
        };
        return block_result_strings[static_cast<unsigned>(result)];
    }

protected:
    template<typename Window, typename Sink, typename Might = ShouldSucceed>
    block_result do_block(Window &window, Sink &sink, const Might &might_tag = {}) {

        /* Starting to read the next block.  */
        if (unlikely(!_in_stream.ensure_bits<1 + 2 + 5 + 5 + 4>())) return block_result::NOT_ENOUGH_INPUT;

        /* BFINAL: 1 bit  */
        block_result success = _in_stream.pop_bits(1) ? block_result::LAST_BLOCK : block_result::SUCCESS;

        /* BTYPE: 2 bits  */
        const libdeflate_decompressor *cur_d;
        switch (_in_stream.pop_bits<uint8_t>(2)) {
            case DEFLATE_BLOCKTYPE_DYNAMIC_HUFFMAN:
                if (might_tag.fail_if(!prepare_dynamic(might_tag))) return block_result::INVALID_DYNAMIC_HT;
                cur_d = &_decompressor;
                break;

            case DEFLATE_BLOCKTYPE_UNCOMPRESSED:
                if (might_tag.fail_if(!do_uncompressed(window, might_tag)))
                    return block_result::INVALID_UNCOMPRESSED_BLOCK;

                return success;

            case DEFLATE_BLOCKTYPE_STATIC_HUFFMAN:
                cur_d = &static_decompressor;
                break;

            default:
                return block_result::INVALID_BLOCK_TYPE;
        }

        /* The main DEFLATE decode loop  */
        for (;;) {
            /* Decode a litlen symbol.  */
            _in_stream.ensure_bits<DEFLATE_MAX_LITLEN_CODEWORD_LEN>();
            // FIXME: entry should be const
            uint32_t entry = cur_d->u.litlen_decode_table[_in_stream.bits<uint16_t>(LITLEN_TABLEBITS)];
            if (entry & HUFFDEC_SUBTABLE_POINTER) {
                /* Litlen subtable required (uncommon case)  */
                _in_stream.remove_bits(LITLEN_TABLEBITS);
                entry = cur_d->u.litlen_decode_table[((entry >> HUFFDEC_RESULT_SHIFT) & 0xFFFF)
                                                     + _in_stream.bits(entry & HUFFDEC_LENGTH_MASK)];
            }
            _in_stream.remove_bits(entry & HUFFDEC_LENGTH_MASK);
            // PRINT_DEBUG("_in_stream position %x\n",_in_stream.in_next);
            if (entry & HUFFDEC_LITERAL) {
                /* Literal  */
                if (unlikely(window.available() == 0)) {
                    if (unlikely(window.flush(sink) == 0)) { return block_result::FLUSH_FAIL; }
                    assert(window.available() != 0);
                }

                if (might_tag.fail_if(!window.push(uint8_t(entry >> HUFFDEC_RESULT_SHIFT)))) {
                    return block_result::INVALID_LITERAL;
                }

                continue;
            }

            /* Match or end-of-block  */
            entry >>= HUFFDEC_RESULT_SHIFT;
            _in_stream.ensure_bits<InputStream::bitbuf_max_ensure>();

            /* Pop the extra length bits and add them to the length base to
             * produce the full length.  */
            const uint32_t length
                    =
                    (entry >> HUFFDEC_LENGTH_BASE_SHIFT) + _in_stream.pop_bits(entry & HUFFDEC_EXTRA_LENGTH_BITS_MASK);
            assert(length <= 258);

            /* The match destination must not end after the end of the
             * output buffer.  For efficiency, combine this check with the
             * end-of-block check.  We're using 0 for the special
             * end-of-block length, so subtract 1 and it turn it into
             * SIZE_MAX.  */
            // static_assert(HUFFDEC_END_OF_BLOCK_LENGTH == 0);
            if (unlikely(length - 1 >= window.available())) {
                if (likely(length == HUFFDEC_END_OF_BLOCK_LENGTH)) { // Block done
                    return success;
                } else { // Needs flushing
                    if (unlikely(window.flush(sink) < length)) { return block_result::FLUSH_FAIL; }
                    assert(length <= window.available());
                }
            }
            assert(length >= 3);

            // if we end up here, it means we're at a match

            /* Decode the match offset.  */
            entry = cur_d->offset_decode_table[_in_stream.bits<uint8_t>(OFFSET_TABLEBITS)];
            if (unlikely(entry & HUFFDEC_SUBTABLE_POINTER)) {
                /* Offset subtable required (uncommon case)  */
                _in_stream.remove_bits(OFFSET_TABLEBITS);
                entry = cur_d->offset_decode_table[((entry >> HUFFDEC_RESULT_SHIFT) & 0xFFFF)
                                                   + _in_stream.bits(entry & HUFFDEC_LENGTH_MASK)];
            }
            _in_stream.remove_bits(entry & HUFFDEC_LENGTH_MASK);
            entry >>= HUFFDEC_RESULT_SHIFT;

            /* Pop the extra offset bits and add them to the offset base to
             * produce the full offset.  */
            const uint32_t offset
                    =
                    (entry & HUFFDEC_OFFSET_BASE_MASK) + _in_stream.pop_bits(entry >> HUFFDEC_EXTRA_OFFSET_BITS_SHIFT);

            /* Copy the match: 'length' bytes at 'window_next - offset' to
             * 'window_next'.  */
            if (might_tag.fail_if(!window.copy_match(length, offset))) { return block_result::INVALID_MATCH; }
        }
    }

private:
    template<typename Might = ShouldSucceed>
    bool prepare_dynamic(const Might &might_tag = {}) {

        /* The order in which precode lengths are stored.  */
        static constexpr uint8_t deflate_precode_lens_permutation[DEFLATE_NUM_PRECODE_SYMS]
                = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};

        /* Read the codeword length counts.  */
        unsigned num_litlen_syms = _in_stream.pop_bits<unsigned>(5) + 257;
        unsigned num_offset_syms = _in_stream.pop_bits<unsigned>(5) + 1;
        const unsigned num_explicit_precode_lens = _in_stream.pop_bits<unsigned>(4) + 4;

        /* Read the precode codeword lengths.  */
        _in_stream.ensure_bits<DEFLATE_NUM_PRECODE_SYMS * 3>();

        for (unsigned i = 0; i < num_explicit_precode_lens; i++)
            _decompressor.u.precode_lens[deflate_precode_lens_permutation[i]] = _in_stream.pop_bits<len_t>(3);

        for (unsigned i = num_explicit_precode_lens; i < DEFLATE_NUM_PRECODE_SYMS; i++)
            _decompressor.u.precode_lens[deflate_precode_lens_permutation[i]] = 0;

        /* Build the decode table for the precode.  */
        if (might_tag.fail_if(!build_precode_decode_table(&_decompressor, might_tag))) return false;

        /* Expand the literal/length and offset codeword lengths.  */
        for (unsigned i = 0; i < num_litlen_syms + num_offset_syms;) {
            _in_stream.ensure_bits<DEFLATE_MAX_PRE_CODEWORD_LEN + 7>();

            /* (The code below assumes that the precode decode table
             * does not have any subtables.)  */
            // static_assert(PRECODE_TABLEBITS == DEFLATE_MAX_PRE_CODEWORD_LEN);

            /* Read the next precode symbol.  */
            const uint32_t entry
                    = _decompressor.u.l.precode_decode_table[_in_stream.bits<uint8_t>(DEFLATE_MAX_PRE_CODEWORD_LEN)];
            _in_stream.remove_bits(entry & HUFFDEC_LENGTH_MASK);
            const unsigned presym = entry >> HUFFDEC_RESULT_SHIFT;

            if (presym < 16) {
                /* Explicit codeword length  */
                _decompressor.u.l.lens[i++] = len_t(presym);
                continue;
            }

            /* Run-length encoded codeword lengths  */

            /* Note: we don't need verify that the repeat count
             * doesn't overflow the number of elements, since we
             * have enough extra spaces to allow for the worst-case
             * overflow (138 zeroes when only 1 length was
             * remaining).
             *
             * In the case of the small repeat counts (presyms 16
             * and 17), it is fastest to always write the maximum
             * number of entries.  That gets rid of branches that
             * would otherwise be required.
             *
             * It is not just because of the numerical order that
             * our checks go in the order 'presym < 16', 'presym ==
             * 16', and 'presym == 17'.  For typical data this is
             * ordered from most frequent to least frequent case.
             */
            if (presym == 16) {
                /* Repeat the previous length 3 - 6 times  */
                if (might_tag.fail_if(!(i != 0))) {
                    PRINT_DEBUG_DECODING("fail at (i!=0)\n");
                    return false;
                }
                const uint8_t rep_val = _decompressor.u.l.lens[i - 1];
                const unsigned rep_count = 3 + _in_stream.pop_bits<unsigned>(2);
                _decompressor.u.l.lens[i + 0] = rep_val;
                _decompressor.u.l.lens[i + 1] = rep_val;
                _decompressor.u.l.lens[i + 2] = rep_val;
                _decompressor.u.l.lens[i + 3] = rep_val;
                _decompressor.u.l.lens[i + 4] = rep_val;
                _decompressor.u.l.lens[i + 5] = rep_val;
                i += rep_count;
            } else if (presym == 17) {
                /* Repeat zero 3 - 10 times  */
                const unsigned rep_count = 3 + _in_stream.pop_bits<unsigned>(3);
                _decompressor.u.l.lens[i + 0] = 0;
                _decompressor.u.l.lens[i + 1] = 0;
                _decompressor.u.l.lens[i + 2] = 0;
                _decompressor.u.l.lens[i + 3] = 0;
                _decompressor.u.l.lens[i + 4] = 0;
                _decompressor.u.l.lens[i + 5] = 0;
                _decompressor.u.l.lens[i + 6] = 0;
                _decompressor.u.l.lens[i + 7] = 0;
                _decompressor.u.l.lens[i + 8] = 0;
                _decompressor.u.l.lens[i + 9] = 0;
                i += rep_count;
            } else {
                /* Repeat zero 11 - 138 times  */
                const unsigned rep_count = 11 + _in_stream.pop_bits<unsigned>(7);
                memset(&_decompressor.u.l.lens[i], 0, rep_count * sizeof(_decompressor.u.l.lens[i]));
                i += rep_count;
            }
        }

        if (!build_offset_decode_table(&_decompressor, num_litlen_syms, num_offset_syms, might_tag)) {
            PRINT_DEBUG_DECODING(
                    "fail at build_offset_decode_table(_decompressor, num_litlen_syms, num_offset_syms)\n");
            return false;
        }
        if (!build_litlen_decode_table(&_decompressor, num_litlen_syms, might_tag)) {
            PRINT_DEBUG_DECODING(
                    "fail at build_litlen_decode_table(_decompressor, num_litlen_syms, num_offset_syms)\n");
            return false;
        }

        return true;
    }

    template<typename OutWindow, typename Might = ShouldSucceed>
    inline bool do_uncompressed(OutWindow &out, const Might &might_tag = {}) {
        /* Uncompressed block: copy 'len' bytes literally from the input
         * buffer to the output buffer.  */
        _in_stream.align_input();

        if (unlikely(_in_stream.available() < 4)) {
            PRINT_DEBUG_DECODING("bad block, uncompressed check less than 4 bytes in input\n");
            return false;
        }

        uint16_t len = _in_stream.pop_u16();
        uint16_t nlen = _in_stream.pop_u16();

        if (might_tag.fail_if(len != uint16_t(~nlen))) {
            // PRINT_DEBUG("bad uncompressed block: len encoding check\n");
            return false;
        }

        if (might_tag.fail_if(len > _in_stream.available())) {
            PRINT_DEBUG_DECODING("bad uncompressed block: len bigger than input stream \n");
            return false;
        }

        if (might_tag.fail_if(!out.copy(_in_stream, len))) {
            PRINT_DEBUG_DECODING("bad uncompressed block: rejected by output window (non-ascii)\n");
            return false;
        }
        return true;
    }

protected:
    InputStream _in_stream;

private:
    static inline struct libdeflate_decompressor make_static_decompressor() {
        struct libdeflate_decompressor sd;
        prepare_static(&sd);
        return sd;
    }

    static const libdeflate_decompressor static_decompressor;
    struct libdeflate_decompressor _decompressor = {};
};

const libdeflate_decompressor DeflateParser::static_decompressor = DeflateParser::make_static_decompressor();

namespace details {

/// Unrolled loops through template recurssion
    template<std::size_t N, typename F, typename R>
    struct repeat_t {
        forceinline_fun static R
        call(F
        f)
        {
            R res = f();
            if (!bool(res))
                return res;
            else
                return repeat_t<N - 1, F, R>::call(f);
        }
    };

    template<std::size_t N, typename F>
    struct repeat_t<N, F, void> {
        forceinline_fun static void call(F f) {
            f();
            repeat_t<N - 1, F, void>::call(f);
        }
    };

    template<typename F, typename R>
    struct repeat_t<0, F, R> {
        forceinline_fun static R
        call(F
        f) { return f(); }
    };

    template<typename F>
    struct repeat_t<0, F, void> {
        forceinline_fun static void call(F f) { f(); }
    };

/// Repeats the operation N times or until the function returns false
    template<std::size_t N, typename FunctionType>
    forceinline_fun auto
    repeat(FunctionType function)

    -> decltype(function()) {
    static_assert(N > 0, "N must be non-zero");
    return
    repeat_t<N - 1, FunctionType, decltype(function())>::call(function);
}

using vec_t = __m128i;
static constexpr size_t vec_size = sizeof(vec_t);

/// Copy size bytes of cache lines without loading the destination in caches
static inline void *
stream_memcpy(void *restrict _dst, const void *restrict _src, size_t size) {
    static_assert(cache_line_size % vec_size == 0, "A integer number of stream_ty should fit in cache line");

    assert(reinterpret_cast<uintptr_t>(_dst) % cache_line_size == 0);
    assert(reinterpret_cast<uintptr_t>(_src) % cache_line_size == 0);
    assert(size % cache_line_size == 0);

    auto dst = reinterpret_cast<vec_t *>(_dst);
    auto src = reinterpret_cast<const vec_t *>(_src);
    auto dst_end = dst + size / vec_size;

    while (repeat<10>([&]() {
        repeat<cache_line_size / vec_size>([&]() { _mm_stream_si128(dst++, _mm_load_si128(src++)); });
        return dst < dst_end;
    }));

    return _dst;
}

/// Copy sizes bytes from _dst - offset to offset, if offset < size, bytes are repeated.
/// Size is rounded up to the next multiple of vec_size (16 bytes currently)
static inline void *
overlap_memcpy(void *_dst, size_t offset, size_t size) {
    static constexpr __v16qi suffle_arr[16] = {
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0,  0,  0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0,  0,  0},
            {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,  1,  0,  1,  0,  1},
            {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1,  2,  0,  1,  2,  0},
            {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2,  3,  0,  1,  2,  3},
            {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0,  1,  2,  3,  4,  0},
            {0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4,  5,  0,  1,  2,  3},
            {0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3,  4,  5,  6,  0,  1},
            {0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2,  3,  4,  5,  6,  7},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1,  2,  3,  4,  5,  6},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0,  1,  2,  3,  4,  5},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0,  1,  2,  3,  4},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0,  1,  2,  3},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 0,  1,  2},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0,  1},
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0},
    };

    static constexpr uint8_t delta[18] = {0, 0, 0, 2, 0, 4, 2, 5, 0, 2, 4, 6, 8, 10, 12, 14, 0, 1};

    uint8_t *dst = static_cast<uint8_t *>(_dst);
    const uint8_t *const dst_end = dst + size;
    const uint8_t *src = dst - offset;

    if (offset < vec_size) {
        vec_t repeats;
        memcpy(&repeats, src, vec_size);
        repeats = _mm_shuffle_epi8(repeats, reinterpret_cast<vec_t>(suffle_arr[offset]));
        memcpy(dst, &repeats, vec_size);

        src = dst - delta[offset];
        dst += vec_size;
        if (dst >= dst_end) return _dst;
    }

    do {
        memcpy(dst, src, vec_size);
        dst += vec_size;
        src += vec_size;
    } while (dst < dst_end);

    return _dst;
}

}

/** Context window for a deflate parse
 * Can be used with different symbols types, alllowing to handle symbolic backreferences
 * The window buffer is wrapped using virtual memory mapping
 */
template<typename _char_t = char, unsigned _context_bits = 15>
class Window {
public:
    // TODO: should we set these at runtime ? context_bits=15 should be fine for all compression levels
    // buffer_size could be adjusted on L3 cache size.
    // But: there is a runtime cost as we loose some optimizations
    static constexpr unsigned context_bits = _context_bits;

    using char_t = _char_t;
    // Ascii range accepted for literals
    static constexpr char_t max_value = char_t('~'), min_value = char_t('\t');

    using wsize_t = uint_fast32_t; /// Type for positive offset in buffer
    using wssize_t = int_fast32_t;  /// Type for signed offsets in buffer

    static constexpr wsize_t context_size = wsize_t(1) << context_bits;
    // We add one 4k page for the actual size of the buffer, allowing to copy more than requested
    static constexpr wsize_t buffer_size = context_size + ::details::page_size / sizeof(char_t);

    static constexpr wsize_t batch_copymatch_size = details::vec_size / sizeof(char_t);
    static_assert(details::vec_size % sizeof(char_t) == 0, "Fractional number of symbols in copy match batches");

    Window(Window &&) noexcept = default;

    Window &operator=(Window &&) noexcept = default;

    Window(const Window &) = delete;

    Window &operator=(const Window &) noexcept = delete;

    Window(const char *shm_name = nullptr)
            : _buffer(alloc_mirrored<char_t>(buffer_size, 3, shm_name)), next(_buffer.end() - buffer_size),
              waterline(_buffer.end() - batch_copymatch_size), _last_flush_end(next) {}

    void clear() {
        next = _buffer.end() - buffer_size;
        _last_flush_end = next;
        waterline = _buffer.end() - batch_copymatch_size;
    }

    wsize_t available() const {
        assert(next <= waterline);
        return wsize_t(waterline - next);
    }

    bool push(char_t c) {
        if (unlikely(char_t(c) > max_value || char_t(c) < min_value)) {
            PRINT_DEBUG("fail, unprintable literal unexpected in fastq\n");
            return false;
        }

        assert(available() >= 1); // Enforced by do_block
        *next++ = char_t(c);
        return true;
    }

    /* return true if it's a reasonable offset, otherwise false */
    bool copy_match(wsize_t length, wsize_t offset) {
        if (unlikely(offset > context_size)) {
            PRINT_DEBUG("fail, copy_match, offset %d\n", (int) offset);
            return false;
        }

        // Could not happen with the way offset and length are encoded
        assert(length >= 3);
        assert(length <= 258);
        assert(offset != 0);
        assert(offset <= context_size);
        assert(available() >= length);

        char_t *dst = next;
        next += length;
        assert(_buffer.includes(dst, dst + details::round_up<batch_copymatch_size>(length)));
        assert(_buffer.includes(dst - offset));

        details::overlap_memcpy(dst, offset * sizeof(char_t), length * sizeof(char_t));

        return true;
    }

    bool copy(InputStream &in, wsize_t length) {
        if (unlikely(!in.check_ascii(length, min_value, max_value))) {
            PRINT_DEBUG("fail, unprintable uncompressed block unexpected in fastq\n");
            return false;
        }
        assert(available() >= length);
        in.copy(next, length);
        next += length;
        return true;
    }

    span<char_t> flushable() {
        span<char_t> flushing = {_last_flush_end, next};
        assert(_buffer.includes(flushing));
        assert(flushing.size() <= buffer_size - batch_copymatch_size);
        return flushing;
    }

    /// Move the 32K context to the start of the buffer
    template<typename Sink>
    size_t flush(Sink &sink) {
        // Waterline is just a cache of this value
        assert(waterline == _last_flush_end + (buffer_size - batch_copymatch_size));

        span<char_t> flushing = flushable();
        size_t sz = sink(flushing);
        assert(sz <= flushing.size());
        _last_flush_end += sz;
        assert(flushing.bounds(_last_flush_end));

        if (_last_flush_end + buffer_size > _buffer.end()) {
            next -= buffer_size;
            _last_flush_end -= buffer_size;
        }
        waterline = _last_flush_end + (buffer_size - batch_copymatch_size);

        // Here why we need 3 contiguous mappings of the same buffer: all theses pointer must fit in the address space
        // of the buffer next - context_size > next > waterline + batch_copymatch_size
        assert(_buffer.includes({next - context_size, next}));
        assert(_buffer.includes({next, waterline + batch_copymatch_size}));
        return sz;
    }

    span<const char_t> current_context() const { return {next - context_size, next}; }

    span<char_t> current_context() { return {next - context_size, next}; }

protected:
    mmap_span<char_t> _buffer; // A buffer_size buffer mapped 3 time contiguously
    char_t *next;
    const char_t *waterline;
    char_t *_last_flush_end;
};

/** A do nothing window for checking the validity of the stream while searching for a block during random access.
 */
struct DummyWindow {
    static constexpr unsigned context_bits = 15;

    using char_t = uint8_t;
    static constexpr char_t max_value = char_t('~'), min_value = char_t('\t');

    using wsize_t = uint_fast32_t; /// Type for positive offset in buffer
    using wssize_t = int_fast32_t;  /// Type for signed offsets in buffer

    static constexpr wsize_t context_size = wsize_t(1) << context_bits;

    void clear() { _size = 0; }

    wsize_t size() const { return _size; }

    wsize_t available() const { return context_size; }

    bool push(char_t c) {
        _size++;
        if (c <= max_value && c >= min_value) {
            return true;
        } else {
            PRINT_DEBUG_DECODING("Non ascii char literal %u\n", unsigned(c));
            return false;
        }
    }

    /* return true if it's a reasonable offset, otherwise false */
    bool copy_match(wsize_t length, wsize_t offset) {
        assert(length >= 3);
        assert(length <= 258);
        assert(offset != 0);
        _size += length;
        if (offset <= context_size) {
            return true;
        } else {
            PRINT_DEBUG_DECODING("invalid LZ77 offset %lu\n", offset);
            return false;
        }
    }

    bool copy(InputStream &in, wsize_t length) {
        _size += length;
        if (in.check_ascii(length, min_value, max_value)) {
            return true;
        } else {
            PRINT_DEBUG_DECODING("Non ascii char in uncompressed block\n");
            return false;
        }
    }

    /// Move the 32K context to the start of the buffer
    size_t flush(DummyWindow &) { return context_size; }

    bool notify_end_block(InputStream &) const { return true; }

protected:
    size_t _size = 0;
};

/** Window that flush its content to a buffer
 * Cache lines are flushed with streaming non-temporal store using the write combining buffer
 */
template<typename char_t>
class SinkBuffer : public span<char_t> {
private:
    static constexpr size_t cache_line_size = details::cache_line_size;
    static_assert(cache_line_size % sizeof(char_t) == 0, "A integer number of char_t should fits in cache line");

public:
    SinkBuffer(span<char_t> buf)
            : span<char_t>(buf) {}

    size_t operator()(span<char_t> data) {
        const size_t sz = size_t(::details::round_down<cache_line_size>(data.end()) - data.begin());
        if (unlikely(sz > this->size())) return 0;

        char_t *dst = this->begin();
        this->pop_front(sz);

        details::stream_memcpy(dst, data.begin(), sz * sizeof(char_t));
        return sz;
    }

    /// Copy the remaining data at the end of decompression, returns the remaining buffer space aligned to the next
    /// cache line
    span<char_t> final_flush(Window<char_t> &window) { // We're leaving together
        span<char_t> last_flush = window.flushable();
        size_t sz = last_flush.size();
        if (unlikely(sz > this->size())) return {}; // Not enough space

        memcpy(this->begin(), last_flush.begin(), sz * sizeof(char_t));
        this->pop_front(sz);

        return {::details::round_up<cache_line_size>(this->begin()), this->end()};
    }

    // Return a pointer after the last written symbol in the buffer
    char_t *buf_ptr() const { return this->begin(); }
};

// Virtual base class for pugz consumers
class ConsumerInterface {
public:
    void set_chunk_idx(unsigned chunk_idx, bool is_last = false) {
        _chunk_idx = chunk_idx;
        _last_chunk = is_last;
    }

    void set_section_idx(unsigned section_idx) { _section_idx = section_idx; }

    unsigned chunk_idx() const { return _chunk_idx; }

    unsigned section_idx() const { return _section_idx; }

    bool is_last_chunk() const { return _last_chunk; }

    size_t operator()(span<const uint8_t> data) {
        flush(data, false);
        return data.size();
    }

    // Flushes the already resolved context in ~32KB step, last indicates if this is the last of the chunk
    virtual void flush(span<const uint8_t> data, bool last) = 0;

    virtual void flush(span<uint16_t> data16bits,
                       span<const uint8_t> lkt16bits,
                       span<uint8_t> data8bits,
                       span<const uint8_t> lkt8bits)
    = 0;

    virtual ~ConsumerInterface() {}

private:
    unsigned _chunk_idx = 0;
    unsigned _section_idx = 0;
    bool _last_chunk = false;
};

/** Compresses the 16bits back-references symbols into 8bits using a lookup-table
 */
template<typename NarrowWindow = Window<uint8_t>, typename WideWindow = Window<uint16_t>>
struct BackrefMultiplexer {

    using narrow_t = typename NarrowWindow::char_t;
    using wide_t = typename WideWindow::char_t;

    static constexpr narrow_t first_backref_symbol = NarrowWindow::max_value + 1;
    static constexpr narrow_t last_backref_symbol = std::numeric_limits<narrow_t>::max();
    static constexpr unsigned total_available_symbols = unsigned(last_backref_symbol) + 1;
    static constexpr size_t context_size = WideWindow::context_size;

    static_assert(WideWindow::context_size == context_size, "Both window should have the same context size");

    BackrefMultiplexer()
            : lkt8to16bits(make_unique_span<wide_t>(total_available_symbols)),
              lkt16bits2chr(make_unique_span<narrow_t>(first_backref_symbol + context_size)),
              lkt8bits2chr(make_unique_span<narrow_t>(total_available_symbols)) {
        for (narrow_t i = 0; i < first_backref_symbol; i++) {
            lkt8to16bits[i] = 0;
        }
        // Prepare the linear part of lookup table
        for (unsigned i = 0; i < first_backref_symbol; i++) {
            lkt16bits2chr[i] = narrow_t(i);
            lkt8bits2chr[i] = narrow_t(i);
        }
    }

    /* Given an input_context with backref encoded as max_value+1+offset on 16bits, try to write
     * a 8bits context in output_context with backref encoded in [max_value+1, 255] through a lookup table
     */
    bool compress_backref_symbols(const WideWindow &input_context, NarrowWindow &output_context) {
        assert(lkt8to16bits);

        narrow_t next_symbol = first_backref_symbol;

        narrow_t *output_p = output_context.current_context().begin();
        for (wide_t c_from : input_context.current_context()) {
            narrow_t c_to;
            if (c_from < first_backref_symbol) {
                c_to = narrow_t(c_from);                        // An in range (resolved) character
            } else {                                            // Or a backref indexing the initial (unknown) context
                c_from = wide_t(c_from - first_backref_symbol); // Get the back-ref offset
                // Linear scan looking for an already allocated backref symbol
                c_to = narrow_t(0);
                for (unsigned i = first_backref_symbol; i < next_symbol; i++) {
                    if (lkt8to16bits[i] == c_from) {
                        c_to = narrow_t(i);
                        assert(c_to != narrow_t(0));
                        break;
                    }
                }
                if (c_to == narrow_t(0)) { // Not found
                    // Try to allocate a new symbol
                    if (next_symbol == 0) { // wrapped arround at previous allocation
                        is_compressed = false;
                        return false;
                    }
                    c_to = narrow_t(next_symbol);
                    lkt8to16bits[next_symbol++] = c_from;
                }
                assert(c_to >= first_backref_symbol);
            }
            *output_p++ = c_to;
        }
        assert(output_p == output_context.current_context().end());

#ifndef NDEBUG // Checks compression
        auto *pcomp = output_context.current_context().begin();
        for (wide_t cin : input_context.current_context()) {
            assert(next_symbol == 0 || *pcomp < next_symbol);
            if (*pcomp < first_backref_symbol) {
                assert(*pcomp == cin);
            } else {
                assert(lkt8to16bits[*pcomp] == cin - first_backref_symbol);
            }
            pcomp++;
        }
        assert(pcomp == output_context.current_context().end());
#endif
        is_compressed = true;
        return true;
    }

    /** Given a context (ie. a lkt: offset -> char) and the compression lkt (8bits code -> offset 16bit) computed by
     * compress_backref_symbols, returns a lkt: 8bits code -> char
     */
    void compose_context(span<const narrow_t> context) {
        assert(context.size() == context_size);

        // Copy the context for the 16bits translation lkt
        narrow_t *lkt16bits2chr_p = &lkt16bits2chr[first_backref_symbol];
        assert(lkt16bits2chr.end() - lkt16bits2chr_p == context_size);
        memcpy(&lkt16bits2chr[first_backref_symbol], context.begin(), context_size * sizeof(narrow_t));

        if (is_compressed) {
            // Compose the compression lookup table to get the tranlation lookup table for 8bit compressed symbols
            for (unsigned i = first_backref_symbol; i < total_available_symbols; i++) {
                wide_t offset = lkt8to16bits[i];
                assert(offset < NarrowWindow::context_size);
                narrow_t chr = context[offset];
                assert(chr >= NarrowWindow::min_value && chr <= NarrowWindow::max_value);
                lkt8bits2chr[i] = chr;
            }
        }
    }

    unique_span<wide_t> lkt8to16bits;  // 8bits code -> offset 16bit
    unique_span<narrow_t> lkt16bits2chr; // 16bits code -> 8bit char = [\0, '~'] + context
    unique_span<narrow_t> lkt8bits2chr;  // 8bits code -> 8bit char = [\0, '~'] + context[lkt8to16bits]
    bool is_compressed = false;
};

class gzip_error : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;

    gzip_error(DeflateParser::block_result res)
            : runtime_error(DeflateParser::block_result_to_cstr(res)) {}
};

/// Monomorphic base for passing information accross threads
class DeflateThread : public DeflateParser {
public:
    static constexpr size_t unset_stop_pos = ~0UL;

    DeflateThread(const InputStream &input_stream, ConsumerInterface &consumer)
            : DeflateParser(input_stream), _consumer(consumer) {}

    // Return the context and the position of the next block in the stream
    std::pair<locked_span<uint8_t>, size_t> get_context() {
        auto lock = std::unique_lock<std::mutex>(_mut);
        while (_state == state_t::RUNNING)
            _cond.wait(lock);

        if (_state == state_t::BORROWED_CONTEXT) {
            span<uint8_t> context;
            std::swap(context, _context);
            // Signal that the context is used: the thread will be resumed once we release the lock
            _state = state_t::RUNNING;
            _cond.notify_all();
            PRINT_DEBUG("%p give context for block %lu\n", (void *) this, _stoped_at);
            return {locked_span<uint8_t>{context, std::move(lock)}, _stoped_at};
        } else {
            assert(_state == state_t::FAIL);
            return {locked_span<uint8_t>{}, unset_stop_pos};
        }
    }

    /// Set the position of the first synced block upstream, so that this thread stops before this block */
    void set_end_block(size_t synced_pos) {
        PRINT_DEBUG("%p set to stop after %lu\n", (void *) this, synced_pos);
        _stop_after.store(synced_pos, std::memory_order_release);
    }

    void set_initial_context(span<uint8_t> context = {}) {
        wait_for_context_borrow();
        _window.clear();
        if (context) { memcpy(_window.current_context().begin(), context.begin(), _window.current_context().size()); }
    }

    // Decompress classically (typically used at position 0) until a certain position
    void go(size_t position_bits = 0) {
        wait_for_context_borrow();

        _in_stream.set_position_bits(position_bits);

        _window.clear();
        PRINT_DEBUG("clear\n");

        auto res = this->decompress_loop(_window, _consumer, []() { return false; });
        PRINT_DEBUG("decompress_loop\n");

        if (res > block_result::CAUGHT_UP_DOWNSTREAM) { throw_gzip_error(res); }

        this->set_context(_window.current_context());
        PRINT_DEBUG("set_context\n");

        _consumer.flush(_window.flushable(), true);
        PRINT_DEBUG("flush\n");
    }

    ~DeflateThread() {
        wait_for_context_borrow(); // Someone might be reading from our windows, so wait before freeing it.
        PRINT_DEBUG("~DeflateThread()");
    }

protected:
    void wait_for_context_borrow() {
        auto lock = std::unique_lock<std::mutex>(_mut);
        bool was_borrowed = _state == state_t::BORROWED_CONTEXT;
        while (_state == state_t::BORROWED_CONTEXT)
            _cond.wait(lock);

        if (was_borrowed) PRINT_DEBUG("%p get context borrow back\n", (void *) this);
    }

    // Post context and wait for it to be consumed
    void set_context(span<uint8_t> ctx) {
#ifndef NDEBUG
        PRINT_DEBUG("%p set context\n", (void *) this);
        for (uint8_t c : ctx) {
            assert(c >= Window<uint8_t>::min_value && c <= Window<uint8_t>::max_value);
        }
#endif

        auto lock = std::unique_lock<std::mutex>(_mut);
        _context = ctx;
        _stoped_at = _in_stream.position_bits();
        assert(_state == state_t::RUNNING);
        _state = state_t::BORROWED_CONTEXT;
        _cond.notify_all();
    }

    // Enter failed state
    void fail() {
        PRINT_DEBUG("%p failed\n", (void *) this);
        auto lock = std::unique_lock<std::mutex>(_mut);
        while (_state == state_t::BORROWED_CONTEXT)
            _cond.wait(lock);

        _state = state_t::FAIL;
        _cond.notify_all();
    }

    template<typename T>
    void throw_gzip_error(T msg) {
        fail();
        throw gzip_error(msg);
    }

    size_t get_stop_pos() const { return _stop_after.load(std::memory_order_acquire); }

    template<typename Window, typename Sink, typename Predicate>
    flatten_fun block_result
    decompress_loop(Window
    & window,
    Sink &sink, Predicate
    && predicate)
    {
        for (;;) {
            if (unlikely(predicate())) return block_result::SUCCESS;
            size_t target_stop = get_stop_pos();

            if (unlikely(_in_stream.position_bits() >= target_stop)) {
                if (_in_stream.position_bits() != target_stop) {
                    PRINT_DEBUG("%p stoped at %lu, %li bits after expected\n",
                                (void *) this,
                                _in_stream.position_bits(),
                                _in_stream.position_bits() - target_stop);
                } else {
                    PRINT_DEBUG("%p stoped at %lu\n", (void *) this, _in_stream.position_bits());
                }
                _stop_after.store(unset_stop_pos, std::memory_order_relaxed);
                return block_result::CAUGHT_UP_DOWNSTREAM;
            }

            block_result res = do_block(window, sink, ShouldSucceed{});

            if (unlikely(res != block_result::SUCCESS)) { return res; }
        }
    }

protected:
    Window<uint8_t> _window = {};
    ConsumerInterface &_consumer;

private:
    /* Members for synchronization and communication */
    std::mutex _mut{};
    std::condition_variable _cond{};
    std::atomic<size_t> _stop_after = {unset_stop_pos}; // Where we should stop
    size_t _stoped_at = unset_stop_pos;   // Where we stopped
    span<uint8_t> _context = {};
    enum class state_t {
        RUNNING, BORROWED_CONTEXT, FAIL
    };
    state_t _state = state_t::RUNNING;
};

constexpr size_t DeflateThread::unset_stop_pos;

class DeflateThreadRandomAccess : public DeflateThread {

public:
    static constexpr size_t buffer_virtual_size = 512ull << 20;

    DeflateThreadRandomAccess(const InputStream &input_stream, ConsumerInterface &consumer)
            : DeflateThread(input_stream, consumer), buffer(alloc_huge<uint8_t>(buffer_virtual_size)) {}

    DeflateThreadRandomAccess(const DeflateThreadRandomAccess &) = delete;

    DeflateThreadRandomAccess &operator=(const DeflateThreadRandomAccess &) = delete;

    ~DeflateThreadRandomAccess() {
        wait_for_context_borrow();
        PRINT_DEBUG("~DeflateThreadRandomAccess\n");
    }

    void set_upstream(DeflateThread *up_stream) { _up_stream = up_stream; }

    // Finds a new block of decompressed size >= min_block_size bits
    // between positions [skip, skip+max_bits_skip] in the compressed stream
    size_t sync(size_t skip,
                const size_t max_bits_skip = size_t(1) << (3 + 20), // 1MiB
                const size_t min_block_size = 1 << 13                // 8KiB
    ) {
        _in_stream.set_position_bits(skip);

        DummyWindow dummy_win;

        size_t pos = skip;
        const size_t max_pos = pos + std::min(8 * _in_stream.size(), max_bits_skip);

        for (_in_stream.ensure_bits<1>(); pos < max_pos; pos++) {
            assert(pos == _in_stream.position_bits());

            if (_in_stream.bits<uint8_t>(1)) { // We don't expect to find a final block
                _in_stream.remove_bits(1);
                _in_stream.ensure_bits<1>();
                continue;
            }

            PRINT_DEBUG_DECODING("trying to decode huffman block at %lu\n", pos);

            block_result res = do_block(dummy_win, dummy_win, ShouldFail{});

            if (unlikely(res == block_result::SUCCESS && dummy_win.size() >= min_block_size)) {
                PRINT_DEBUG("%p Candidate block start at %lubits\n", (void *) this, pos);
                _in_stream.set_position_bits(pos);
                _up_stream->set_end_block(pos); // This is not even needed !
                return pos;
            }

            dummy_win.clear();
            if (unlikely(!_in_stream.set_position_bits(pos + 1))) { break; }
        }
        return 8 * _in_stream.size();
    }

    // Decompress a chunk starting at position "skipbits" (in bits) in the compressed stream
    // will guess (by calling sync()) the position of the next block
    bool go(size_t skipbits) {
        assert(_up_stream != nullptr);

        size_t sync_bitpos = sync(skipbits);

        // Get the bit position where the chunk stops. Previously it came from the thread handling the upstream chunk,
        // now it is set up deterministically from go()'s caller.
        // FIXME: With the new setup, this could be combined with the previous check and checked in sync()
        size_t stop_bitpos = get_stop_pos();
        if (stop_bitpos != unset_stop_pos && sync_bitpos >= stop_bitpos) {
            // FIXME: We found our first block after where we are supposed to stop; Could be due to the
            // file being multi-part (hence not yet supported by pugz)
            throw_gzip_error("Failed to find a gzip block during random access");
        }

        // Prepare the wide_window
        span<uint16_t> wide_buffer = buffer.reinterpret<uint16_t>();
        SinkBuffer<uint16_t> wide_sink = wide_buffer;
        wide_window.clear();
        uint16_t sym = wide_window.max_value + 1;
        for (auto &c : wide_window.current_context())
            c = sym++;

        // Decompress to 16bits buffer until there is a small enough number of back-references
        multiplexer.is_compressed = false;
        size_t block_count = 0;
        auto res = decompress_loop(wide_window, wide_sink, [&]() {
            block_count++;
            if (block_count <= 8 || block_count % 2 == 0) return false;

            wait_for_context_borrow(); // the narrow_window is mutated from this point
            _window.clear();
            return multiplexer.compress_backref_symbols(wide_window, _window);
        });

        // Seal the 16bit buffer, and get the remaining for the 8bit part
        span<uint8_t> narrow_buffer = wide_sink.final_flush(wide_window).reinterpret<uint8_t>();
        wide_buffer = {wide_buffer.begin(), wide_sink.begin()};

        if (res == block_result::SUCCESS) {
            // Decompress to the 8bit buffer
            SinkBuffer<uint8_t> narrow_sink = narrow_buffer;
            res = this->decompress_loop(_window, narrow_sink, []() { return false; });

            if (res == block_result::CAUGHT_UP_DOWNSTREAM || res == block_result::LAST_BLOCK) {
                // Seal the narrow buffer
                narrow_sink.final_flush(_window);
                narrow_buffer = {narrow_buffer.begin(), narrow_sink.begin()};

                // Get the context and prepare lookup table
                if (!prepare_lookup_table(sync_bitpos)) return false;

                // Translate the context for the next block
                for (auto &c : _window.current_context()) {
                    c = multiplexer.lkt8bits2chr[c];
                    assert(c >= _window.min_value && c <= _window.max_value);
                }
                this->set_context(_window.current_context()); // From now on, the window is borrowed

                _consumer.flush(wide_buffer, multiplexer.lkt16bits2chr, narrow_buffer, multiplexer.lkt8bits2chr);
            } else if (res == block_result::FLUSH_FAIL) {
                throw_gzip_error(res);
                // FIXME: buffer too small for input buffer_virtual_size = 512MiB for 32MiB of input => max compression
                // ratio of 16x. The simplest way to handle that is stop everything and retry with smaller
                // sections/chunks
            } else {
                throw_gzip_error(res);
                // At this point we have narrowed down the back-reference count to 126 and decoded more than 8 block, so
            }

        } else if (res == block_result::CAUGHT_UP_DOWNSTREAM || res == block_result::LAST_BLOCK) {
            wait_for_context_borrow(); // the narrow_window is mutated from this point
            _window.clear();

            // Get the context and prepare lookup table
            if (!prepare_lookup_table(sync_bitpos)) return false;

            // Translate the context for the next block
            auto *p = wide_window.current_context().begin();
            for (auto &c : _window.current_context()) {
                c = multiplexer.lkt16bits2chr[*p++];
                assert(c >= _window.min_value && c <= _window.max_value);
            }
            assert(p == wide_window.current_context().end());

            this->set_context(_window.current_context()); // From now on, the window is borrowed

            _consumer.flush(wide_buffer, multiplexer.lkt16bits2chr, {}, {});
        } else if (res == block_result::FLUSH_FAIL) {
            throw_gzip_error(res); // FIXME: buffer overflow, see above
        } else {
            throw_gzip_error(res);
            // FIXME: parse error in the first blocks (before compressing the backrefs to 8bits)
            // Maybe the synchronisation is a false positive and we should retry at sync_bitpos + 1
        }

        if (res == block_result::LAST_BLOCK) {
            PRINT_DEBUG("%p last block at %lu\n", (void *) this, _in_stream.position_bits());
        }
        return true;
    }

private:
    bool prepare_lookup_table(size_t sync_bitpos) {
        auto upstream_context = _up_stream->get_context();
        if (upstream_context.second == unset_stop_pos) {
            // Upstream decompressor failed
            fail();
            return false;
        }
        // Check if the context position we got match with our start position
        if (upstream_context.second != sync_bitpos) {
            // FIXME: redecompress from this position (with resolved context this time)
            throw_gzip_error("Got a context from invalid position");
        }
        multiplexer.compose_context(upstream_context.first);
        return true;
    }

    malloc_span<uint8_t> buffer;
    Window<uint16_t> wide_window = {};
    BackrefMultiplexer<Window<uint8_t>, Window<uint16_t>> multiplexer = {};
    DeflateThread *_up_stream = nullptr;
};

class ConsumerSync {
public:
    using lock_t = std::unique_lock<std::mutex>;

    void wait(ConsumerInterface &consumer) {
        lock_t lock{_mut};
        while (_section_idx != consumer.section_idx() || _chunk_idx != consumer.chunk_idx())
            _cond.wait(lock);
        lock.release();
    }

    void notify(ConsumerInterface &consumer) {
        if (not consumer.is_last_chunk()) {
            _chunk_idx++;
        } else {
            _chunk_idx = 0;
            _section_idx++;
        }
        _cond.notify_all();
        _mut.unlock();
    }

private:
    std::mutex _mut = {};
    std::condition_variable _cond = {};

    unsigned _section_idx = 0;
    unsigned _chunk_idx = 0;
};

template<typename Consumer>
class ConsumerWrapper : public ConsumerInterface {
public:
    ConsumerWrapper(Consumer &consumer, ConsumerSync *sync = nullptr)
            : _consumer(consumer), _sync(sync) {}

    ConsumerWrapper(const ConsumerWrapper &) noexcept = default;

    ConsumerWrapper &operator=(const ConsumerWrapper &) noexcept = default;

protected:
    // Flushes the already resolved context in ~32KB step, last indicates if this is the last of the chunk
    virtual void flush(span<const uint8_t> data, bool last) {
        if (not last) {
            if (unlikely(_sync != nullptr && _resolved_idx == 0)) _sync->wait(*this);

            _consumer(data);

            _resolved_idx++;
        } else {
            _consumer(data);

            _resolved_idx = 0;
            if (_sync != nullptr) _sync->notify(*this);
        }
    }

    virtual void flush(span<uint16_t> data16bits,
                       span<const uint8_t> lkt16bits,
                       span<uint8_t> data8bits,
                       span<const uint8_t> lkt8bits) {
        if (_sync != nullptr) _sync->wait(*this);

        slice_span(data16bits, 16 << 10, [&](span<uint16_t> slice) {
            uint8_t *s = reinterpret_cast<uint8_t *>(slice.begin());
            uint8_t *p = s;
            for (auto sym : slice)
                *p++ = lkt16bits[sym];
            _consumer(span<const uint8_t>(s, p));
        });

        slice_span(data8bits, 32 << 10, [&](span<uint8_t> slice) {
            for (auto &sym : slice)
                sym = lkt8bits[sym];
            _consumer(span<const uint8_t>(slice));
        });
        if (_sync != nullptr) _sync->notify(*this);
    }

private:
    template<typename T, typename F>
    static void slice_span(span<T> data, size_t n, F f) {
        T *start = data.begin();
        T *end = start + n;
        for (;; start = end, end += n) {
            if (end < data.end()) {
                f(span<T>(start, end));
            } else {
                f(span<T>(start, data.end()));
                return;
            }
        }
    }

    Consumer &_consumer;
    ConsumerSync *_sync = nullptr;
    unsigned _resolved_idx = 0;
};

struct OutputConsumer {

    //    void operator()(span<const uint8_t> data) const { write(STDOUT_FILENO, data.begin(), data.size()); }
    void operator()(span<const uint8_t> data) {
        if (*pDone == 1) {
            return;
        }
//        printf("ready to put data to memory111\n");
        char *utmp = new char[data.size()];
        memcpy(utmp, data.begin(), data.size());

//        printf("ready to put data to memory222\n");
//        int cntt = 0;
        while (P->try_enqueue({utmp, data.size()}) == 0) {
            if (*pDone == 1) {
                return;
            }
//            printf("pugz waitint\n");
//            cout << "pugz " << num << " wait producer" << " queue size " << P->size_approx() << endl;
            usleep(100);
//            cout << "sleep " << cntt++ << endl;
        }
//        if (cntt)cout << "sleep done" << endl;
//        if (P->size_approx() % 100 == 0) {
//            cout << " queue " << num << " size " << P->size_approx() << endl;
//        }
    }

    moodycamel::ReaderWriterQueue<std::pair<char *, int>> *P;
    int num;
    atomic_int *pDone;
};

struct LineCounter {
    void operator()(span<const uint8_t> data) {
        size_t count = 0;
        const uint8_t *p = data.begin();
        const uint8_t *e = data.end();

        for (;;) {
            p = static_cast<const uint8_t *>(memchr(p, static_cast<int>('\n'), size_t(e - p)));
            if (p != nullptr) {
                count += 1;
                p++;
            } else {
                break;
            }
        }

        lines.fetch_add(count);
    }

    ~LineCounter() { fprintf(stdout, "%lu\n", lines.load()); }

    std::atomic<size_t> lines = {0};
};

/* namespace */
