#ifndef GATBL_COMMON_HPP
#define GATBL_COMMON_HPP

#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <exception>

#define packed_layout __attribute__((__packed__))
#define weak_sym __attribute__((weak))
#define restrict __restrict__

#ifndef NOCOMPILER_HINTS
#    define noinline_fun __attribute__((noinline))
#    define forceinline_fun inline __attribute__((always_inline))
#    define flatten_fun __attribute__((flatten))
#    define pure_fun __attribute__((const))
#    define hot_fun __attribute__((hot))
#    define cold_fun __attribute__((cold))
#    define likely(expr) __builtin_expect(!!(expr), 1)
#    define unlikely(expr) __builtin_expect(!!(expr), 0)
#    define prefetchr(addr) __builtin_prefetch((addr), 0)
#    define prefetchw(addr) __builtin_prefetch((addr), 1)
#else
#    define noinline_fun
#    define forceinline_fun inline
#    define flatten_fun
#    define pure_fun
#    define hot_fun
#    define cold_fun
#    define likely(expr) (expr)
#    define unlikely(expr) (expr)
#    define prefetchr(addr) (expr)
#    define prefetchw(addr) (expr)
#endif

#ifndef __has_feature
#    define __has_feature(x) (0)
#endif

#ifdef __has_cpp_attribute
#    if __has_cpp_attribute(noreturn)
#        define noreturn_attr [[noreturn]]
#    endif
#    if __has_cpp_attribute(nodiscard)
#        define nodiscard_attr [[nodiscard]]
#    elif __has_cpp_attribute(gnu::warn_unused_result)
#        define nodiscard_attr [[gnu::warn_unused_result]]
#    endif
#endif

#ifndef noreturn_attr
#    define noreturn_attr __attribute__((noreturn))
#endif

#ifndef nodiscard_attr
#    define nodiscard_attr __attribute__((warn_unused_result))
#endif

#ifdef NDEBUG
#    define DEBUG 0
#    define gatbl_assert(...) (static_cast<void>(0))
#    define gatbl_assume(expr, ...) (likely((expr)) ? static_cast<void>(0) : __builtin_unreachable())
#else
#    define DEBUG 1

namespace gatbl {

/** Handler for assertion/assumption fails
 * Do not throw exception on purpose (directly terminate)
 */
noreturn_attr noinline_fun inline void
              abort_message(const char* msg...)
{
    va_list args;
    va_start(args, msg);
    std::vfprintf(stderr, msg, args);
    va_end(args);
    std::fflush(stderr);
    std::abort();
}

}

#    define __gatbl_sourceloc_fail(what, msg, ...) gatbl::abort_message("%s:%u %s\n\t" #what " failed: " #msg "\n", __FILE__, __LINE__, static_cast<const char*>(__PRETTY_FUNCTION__), ##__VA_ARGS__))
#    define gatbl_assert(expr, ...) (likely((expr)) ? static_cast<void>(0) : __gatbl_sourceloc_fail("Assertion '" #expr "'", ##__VA_ARGS__)
#    define gatbl_assume(expr, ...) (likely((expr)) ? static_cast<void>(0) : __gatbl_sourceloc_fail("Assumption '" #expr "'", ##__VA_ARGS__)

#endif

#ifdef __cpp_lib_byte
using byte = std::byte;
#else
enum class byte : unsigned char {};
#endif

/*
 * Word type of the target architecture.  Use 'size_t' instead of 'unsigned
 * long' to account for platforms such as Windows that use 32-bit 'unsigned
 * long' on 64-bit architectures.
 */
typedef size_t machine_word_t;

#if defined(PRINT_DEBUG) && PRINT_DEBUG
#    undef PRINT_DEBUG
#    define PRINT_DEBUG(...)                                                                                           \
        {                                                                                                              \
            fprintf(stderr, __VA_ARGS__);                                                                              \
            fflush(stderr);                                                                                            \
        }
#else
#    undef PRINT_DEBUG
#    define PRINT_DEBUG(...) static_cast<void>(0)
#endif

#if defined(PRINT_DEBUG_DECODING) && PRINT_DEBUG_DECODING
#    undef PRINT_DEBUG_DECODING
#    define PRINT_DEBUG_DECODING(...)                                                                                  \
        {                                                                                                              \
            fprintf(stderr, __VA_ARGS__);                                                                              \
            fflush(stderr);                                                                                            \
        }
#else
#    undef PRINT_DEBUG_DECODING
#    define PRINT_DEBUG_DECODING(...) static_cast<void>(0)
#endif

#endif // GATBL_COMMON_HPP
