#ifndef SYSTEM_H
#define SYSTEM_H

#include <cstdlib>

#include <exception>
#include <system_error>
#include <string>

namespace utils { namespace system { namespace except {

void
throw_syserr(const std::string& what)
{
    int errcode = errno;
    errno       = 0;
    switch (errcode) {
        case ENOMEM: throw std::bad_alloc();
        default: throw std::system_error(errcode, std::generic_category(), what);
    }
}

inline int
check_ret(int ret, const std::string& what)
{
    if (ret < 0) throw_syserr(what);
    return ret;
}

template<typename T>
inline T*
check_ptr(T* ret, const std::string& what)
{
    if (ret == nullptr) throw_syserr(what);
    return ret;
}
}}}

#endif // SYSTEM_H
