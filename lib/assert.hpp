#ifndef ASSERT_HPP
#define ASSERT_HPP

#include <cstdio>

#include "common.hpp"

#ifdef NDEBUG
#    define assert(expr) (likely((expr)) ? static_cast<void>(0) : __builtin_unreachable())
#else
#    define assert(expr)                                                                                               \
        (likely((expr)) ? static_cast<void>(0) : __assert_fail(#expr, __FILE__, __LINE__, __PRETTY_FUNCTION__))
#endif

[[noreturn]] inline __attribute__((noinline)) void
__assert_fail(const char* assertion, const char* file, unsigned int line, const char* function) noexcept
{
    std::fprintf(stderr, "%s:%u: Assertion '%s' failed in '%s'.\n", file, line, assertion, function);
    std::fflush(stderr);
    std::abort();
}

/// Allows to bias branch weight based on template parameter
struct ShouldSucceed
{
    static constexpr bool succeed_if(bool p) { return likely(p); }
    static constexpr bool fail_if(bool p) { return unlikely(p); }
};
struct ShouldFail
{
    static constexpr bool succeed_if(bool p) { return unlikely(p); }
    static constexpr bool fail_if(bool p) { return likely(p); }
};

#endif // ASSERT_HPP
