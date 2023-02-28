#ifndef CARE_DESERIALIZE_HPP
#define CARE_DESERIALIZE_HPP

#include <cstring>
#include <cstdint>
#include <fstream>
#include <string>

namespace care {

// De-serialization helpers

// The implementations are they way they are to be absolutely 100% standard-compliant no matter how read() is implemented.
// The char buffer will (usually) be optimized out.

template<typename T>
inline T& read_one(std::ifstream& is, T& v) {
    char tmp[sizeof(T)];
    is.read(tmp, sizeof(T));
    std::memcpy(&v, tmp, sizeof(T));
    return v;
}

template<typename T>
inline T read_one(std::ifstream& is) {
    T ret;
    read_one(is, ret);
    return ret;
}

template<typename T = std::uint64_t>
inline std::string read_str(std::ifstream& is) {
    std::string ret(read_one<T>(is), char());
    is.read(&ret[0], ret.size());
    return ret;
}

} // namespace care

#endif