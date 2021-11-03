#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <cstdlib>

#include <exception>
#include <system_error>

// FIXME: separate assert.h
#include <sstream>
#include <iostream>

#include "common/common.hpp"

namespace gatbl { namespace sys {
// The oxymoric "noinline inline" means that we want a weak non inlineable function
noreturn_attr noinline_fun inline cold_fun void
              throw_syserr(const char fmt[]...) // No variadic template, avoiding to generate too much code
{
    va_list     args;
    std::string _what;

    int errcode = 0;
    std::swap(errcode, errno);

    int size = vsnprintf(nullptr, 0, fmt, args);
    if (size < 0) { std::terminate(); }
    _what.resize(static_cast<size_t>(size));

    throw std::system_error(errcode, std::generic_category(), _what);
}

template<typename... Args>
forceinline_fun hot_fun unsigned
                check_ret(int ret, unsigned min, const char* what, Args&&... args)
{
    if (unlikely(ret < int(min))) { throw_syserr(what, std::forward<Args>(args)...); }
    return unsigned(ret);
}

template<typename... Args>
forceinline_fun hot_fun unsigned
                check_ret(int ret, const char* what, Args&&... args)
{
    return check_ret(ret, 0u, what, std::forward<Args>(args)...);
}

template<typename T, typename... Args>
forceinline_fun hot_fun T*
                        check_ptr(T* ret, const char* what, Args&&... args)
{
    if (unlikely(ret == nullptr)) { throw_syserr(what, std::forward<Args>(args)...); }
    return ret;
}

} // namespace sys
} // namespace gatbl

#endif // EXCEPTIONS_H
