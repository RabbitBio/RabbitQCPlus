#ifndef MYSTRINGVIEW_HPP
#define MYSTRINGVIEW_HPP

#include <hpc_helpers.cuh>

#include <string>

#if __cplusplus >= 201703L
#define STD_STRING_VIEW_CAN_BE_USED
#endif

#ifdef STD_STRING_VIEW_CAN_BE_USED

#include <string_view>

#endif

struct MyStringView{

    constexpr MyStringView() = default;

    HOSTDEVICEQUALIFIER
    constexpr MyStringView(const char* p, int l) noexcept : MyStringView{p,p+l}{}

    HOSTDEVICEQUALIFIER
    constexpr MyStringView(const char* pb, const char* pe) noexcept : ptrbegin{pb}, ptrend{pe}{}

    HOSTDEVICEQUALIFIER
    constexpr MyStringView(const MyStringView& rhs) noexcept : MyStringView{rhs.ptrbegin, rhs.ptrend}{}

    HOSTDEVICEQUALIFIER
    constexpr MyStringView operator=(const MyStringView& rhs) noexcept{
        ptrbegin = rhs.ptrbegin;
        ptrend = rhs.ptrend;
        return *this;
    }

    MyStringView(const std::string& s) : MyStringView{s.data(), int(s.size())}{}

    #ifdef STD_STRING_VIEW_CAN_BE_USED

    MyStringView(const std::string_view& sv) : MyStringView{sv.data(), int(sv.size())}{}

    constexpr MyStringView operator=(const std::string_view& rhs) noexcept{
        ptrbegin = rhs.data();
        ptrend = rhs.data() + rhs.size();
        return *this;
    }

    #endif

    HOSTDEVICEQUALIFIER
    constexpr const char* begin() const noexcept{
        return ptrbegin;
    }

    HOSTDEVICEQUALIFIER
    constexpr const char* end() const noexcept{
        return ptrend;
    }

    HOSTDEVICEQUALIFIER
    constexpr char operator[](int i) const noexcept{
        return data()[i];
    }

    HOSTDEVICEQUALIFIER
    constexpr char at(int i) const noexcept{
        assert(i < length());
        return data()[i];
    }

    HOSTDEVICEQUALIFIER
    constexpr char front() const noexcept{
        return operator[](0);
    }

    HOSTDEVICEQUALIFIER
    constexpr char back() const noexcept{
        return operator[](length() - 1);
    }

    HOSTDEVICEQUALIFIER
    constexpr const char* data() const noexcept{
        return ptrbegin;
    }

    HOSTDEVICEQUALIFIER
    constexpr int size() const noexcept{
        return ptrend - ptrbegin;
    }

    HOSTDEVICEQUALIFIER
    constexpr int length() const noexcept{
        return size();
    }

    HOSTDEVICEQUALIFIER
    constexpr bool empty() const noexcept{
        return length() == 0;
    }

    HOSTDEVICEQUALIFIER
    constexpr void remove_prefix(int n) noexcept{
        ptrbegin += n;
    }

    HOSTDEVICEQUALIFIER
    constexpr void remove_suffix(int n) noexcept{
        ptrend -= n;
    }

    const char* ptrbegin{};
    const char* ptrend{};
};







template<class T, class Difference_Type = int, class Stride = int>
struct StrideIterator {
public:

    using value_type = T;
    using pointer = T*;
    using reference = T&;
    using difference_type = Difference_Type;

public:
    StrideIterator() = default;
        
    HOSTDEVICEQUALIFIER
    StrideIterator(T* ptr, Stride stride_) : ptrbegin{ptr}, stride{stride_}{}
        
    StrideIterator(const StrideIterator& rhs) = default;

    StrideIterator& operator=(const StrideIterator& rhs) = default;
    
    HOSTDEVICEQUALIFIER
    operator bool() const
    {
        if(ptrbegin)
            return true;
        else
            return false;
    }
    
    HOSTDEVICEQUALIFIER
    bool operator==(const StrideIterator& rhs) const noexcept{
        return ptrbegin == rhs.ptrbegin && stride == rhs.stride;
    }
    
    HOSTDEVICEQUALIFIER
    bool operator!=(const StrideIterator& rhs) const{
        return !(operator==(rhs));
    }
    
    HOSTDEVICEQUALIFIER
    StrideIterator& operator+=(const difference_type& i) noexcept{
        ptrbegin += i * getStride();
        return *this;
    }
        
    HOSTDEVICEQUALIFIER
    StrideIterator& operator-=(const difference_type& i) noexcept{
        ptrbegin -= i * getStride();
        return *this;
    }
    
    HOSTDEVICEQUALIFIER
    StrideIterator& operator++() noexcept{
        operator+=(1);        
        return *this;
    }
    
    HOSTDEVICEQUALIFIER
    StrideIterator& operator--() noexcept{
        operator-=(1);
        return *this;
    }
    
    HOSTDEVICEQUALIFIER
    StrideIterator operator++(int) noexcept{
        auto temp(*this);
        operator+=(1);
        return temp;
    }
    
    HOSTDEVICEQUALIFIER
    StrideIterator operator--(int) noexcept{
        auto temp(*this);
        operator-=(1);
        return temp;
    }
    
    HOSTDEVICEQUALIFIER    
    StrideIterator operator+(const difference_type& i) const noexcept{
        auto temp(*this);
        temp.operator+=(i);
        return temp;
    }
    
    HOSTDEVICEQUALIFIER
    StrideIterator operator-(const difference_type& i) const noexcept{
        auto temp(*this);
        temp.operator-=(i);
        return temp;
    }
    
    HOSTDEVICEQUALIFIER
    difference_type operator-(const StrideIterator& rhs){
        assert(isSameStride(rhs));

        return (rhs.ptrbegin - ptrbegin) / getStride();
    }
    
    HOSTDEVICEQUALIFIER
    T& operator*(){
        return *ptrbegin;
    }
    
    HOSTDEVICEQUALIFIER
    const T& operator*() const{
        return *ptrbegin;
    }
    
    HOSTDEVICEQUALIFIER
    T* operator->(){
        return ptrbegin;
    }
    
    HOSTDEVICEQUALIFIER
    Stride getStride() const noexcept{
        return stride;
    }

private:
    HOSTDEVICEQUALIFIER
    void isSameStride(const StrideIterator& rhs) const noexcept{
        return getStride() == rhs.getStride();
    }

    Stride stride;
    T* ptrbegin;
};








struct MyStridedStringView{
public:
    using value_type = char;
    using size_type = int;
    using difference_type = size_t;
    using stride_type = int;
    using const_iterator = StrideIterator<const value_type, difference_type, stride_type>;
    using iterator = StrideIterator<value_type, difference_type, stride_type>;

public:

    constexpr MyStridedStringView() = default;

    HOSTDEVICEQUALIFIER
    constexpr MyStridedStringView(const char* p, int l, int o) noexcept : MyStridedStringView{p, p + l*o, o}{}

    HOSTDEVICEQUALIFIER
    constexpr MyStridedStringView(const char* pb, const char* pe, int o) noexcept : ptrbegin{pb}, ptrend{pe}, offset{o}{}

    HOSTDEVICEQUALIFIER
    constexpr MyStridedStringView(const MyStridedStringView& rhs) noexcept : MyStridedStringView{rhs.ptrbegin, rhs.ptrend, rhs.offset}{}

    HOSTDEVICEQUALIFIER
    constexpr MyStridedStringView operator=(const MyStridedStringView& rhs) noexcept{
        ptrbegin = rhs.ptrbegin;
        ptrend = rhs.ptrend;
        offset = rhs.offset;
        return *this;
    }

    MyStridedStringView(const std::string& s) : MyStridedStringView{s.data(), int(s.size()), 1}{}

    #ifdef STD_STRING_VIEW_CAN_BE_USED

    MyStridedStringView(const std::string_view& sv) : MyStridedStringView{sv.data(), int(sv.size()), 1}{}

    constexpr MyStridedStringView operator=(const std::string_view& rhs) noexcept{
        ptrbegin = rhs.data();
        ptrend = rhs.data() + rhs.size();
        offset = 1;
        return *this;
    }

    #endif

    HOSTDEVICEQUALIFIER
    const_iterator begin() noexcept{
        return const_iterator{ptrbegin, stride()};
    }

    HOSTDEVICEQUALIFIER
    const_iterator end() const noexcept{
        return const_iterator{ptrend, stride()};
    }

    HOSTDEVICEQUALIFIER
    constexpr value_type operator[](int i) const noexcept{
        return data()[i * offset];
    }

    HOSTDEVICEQUALIFIER
    constexpr value_type at(int i) const noexcept{
        assert(i < length());
        return operator[](i);
    }

    HOSTDEVICEQUALIFIER
    constexpr value_type front() const noexcept{
        return operator[](0);
    }

    HOSTDEVICEQUALIFIER
    constexpr value_type back() const noexcept{
        return operator[](length() - 1);
    }

    HOSTDEVICEQUALIFIER
    constexpr const value_type* data() const noexcept{
        return ptrbegin;
    }

    HOSTDEVICEQUALIFIER
    constexpr int size() const noexcept{
        return (ptrend - ptrbegin) / offset;
    }

    HOSTDEVICEQUALIFIER
    constexpr int length() const noexcept{
        return size();
    }

    HOSTDEVICEQUALIFIER
    constexpr bool empty() const noexcept{
        return length() == 0;
    }

    HOSTDEVICEQUALIFIER
    constexpr void remove_prefix(int n) noexcept{
        ptrbegin += n * offset;
    }

    HOSTDEVICEQUALIFIER
    constexpr void remove_suffix(int n) noexcept{
        ptrend -= n * offset;
    }

    HOSTDEVICEQUALIFIER
    constexpr int stride() const noexcept{
        return offset;
    }

    int offset{};
    const value_type* ptrbegin{};
    const value_type* ptrend{};
};


#ifdef STD_STRING_VIEW_CAN_BE_USED
#undef STD_STRING_VIEW_CAN_BE_USED
#endif

#endif