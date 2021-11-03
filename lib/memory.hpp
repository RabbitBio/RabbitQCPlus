#ifndef MEMORY_HPP
#define MEMORY_HPP

#include <cstdlib>
#include <cstdint>
#include <memory>
#include <mutex>

#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h>
#include <sys/mman.h>

#include "assert.hpp"
#include "common/exceptions.hpp"

namespace sys = gatbl::sys;

namespace details {

// FIXME: we should ask the system, but many uses expect compile time constants...
// Beside, I'm not aware of a linux x86_64 system that differ from this.
static constexpr size_t cache_line_size = 64;
static constexpr size_t page_size       = 1ull << 12;
static constexpr size_t huge_page_size  = 1ull << 21;

template<size_t n, typename T>
inline T
round_down(T i)
{
    return i & T(-n);
}

template<size_t n, typename T>
inline T
round_up(T i)
{
    return (i + (n - 1)) & T(-n);
}

template<size_t n, typename T>
inline T*
round_down(T* p)
{
    return reinterpret_cast<T*>(round_down<n>(reinterpret_cast<uintptr_t>(p)));
}

template<size_t n, typename T>
inline T*
round_up(T* p)
{
    return reinterpret_cast<T*>(round_up<n>(reinterpret_cast<uintptr_t>(p)));
}

template<typename T>
inline void
destruct(T* p)
{
    if (std::is_trivially_destructible<T>::value) return;
    p->~T();
}

template<typename T>
inline void
destruct(T* p, T* e)
{
    if (std::is_trivially_destructible<T>::value) return;
    for (; p < e; ++p)
        p->~T();
}

// Allows to call our deleters specialized for arrays without in-band size metadata
template<typename D, typename T>
inline auto
invoke_deleter(D& d, T* ptr, size_t sz) -> decltype(d(ptr, sz))
{
    return d(ptr, sz);
}

template<typename D, typename T>
inline auto
invoke_deleter(D& d, T* ptr, size_t) -> decltype(d(ptr))
{
    return d(ptr);
}

/// Create an anonymous file in a tmpfs
inline int
tmpshm(const char* shm_name_base)
{
    shm_name_base = shm_name_base == nullptr ? "/mirrored_buffer" : shm_name_base;
    for (unsigned retries = 0; retries < 1024; retries++) { // Try different file names
        std::string shm_name{shm_name_base};
        if (retries > 0) {
            shm_name += std::to_string(retries - 1);
            int ret = shm_open(shm_name.c_str(), O_RDWR | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
            if (ret > 0) {
                sys::check_ret(shm_unlink(shm_name.c_str()), "shm_unlink");
                return ret;
            } else if (errno != EEXIST) {
                sys::throw_syserr("shm_open");
            }
        }
    }
    sys::throw_syserr("shm_open");
}

} // namespace details

template<typename T> class span
{
  public:
    using element_type    = T;
    using pointer         = T*;
    using const_pointer   = T*;
    using iterator        = T*;
    using const_iterator  = T*; // Not const because we don't own the range.
    using reference       = T&;
    using const_reference = T&;

    span() noexcept
      : _begin(nullptr)
      , _end(nullptr)
    {}

    span(T* start, T* last) noexcept
      : _begin(start)
      , _end(last)
    {
        assert(_end >= _begin);
    }

    span(T* start, size_t n) noexcept
      : _begin(start)
      , _end(start + n)
    {}

    template<typename R,
             typename
             = typename std::enable_if<std::is_convertible<decltype(std::declval<R>().begin()), T*>::value>::type>
    span(R&& r) noexcept
      : _begin(r.begin())
      , _end(r.end())
    {
        assert(_end >= _begin);
    }

    span(const span&) noexcept = default;
    span(span&&) noexcept      = default;
    span& operator=(const span&) noexcept = default;
    span& operator=(span&&) noexcept = default;

    size_t size() const restrict
    {
        assert(_end >= _begin);
        return size_t(_end - _begin);
    }
    bool empty() const restrict { return size() == 0; }
         operator bool() const restrict { return size() != 0; }

    reference pop_front(size_t n = 1) restrict
    {
        _begin += n;
        assert(_end >= _begin);
        return *_begin;
    }

    reference pop_back(size_t n = 1) restrict
    {
        _end -= n;
        assert(_end >= _begin);
        return *_end;
    }

    iterator begin() const restrict { return _begin; }
    iterator end() const restrict { return _end; }

    bool                      bounds(const T* p) const restrict { return p >= _begin && p <= _end; }
    bool                      includes(const T* p) const restrict { return p >= _begin && p < _end; }
    bool                      includes(const T* b, const T* e) const restrict { return b >= _begin && e <= _end; }
    template<typename R> auto includes(const R& r) -> decltype(includes(r.begin, r.end)) const
    {
        return r.begin() >= _begin && r.end() <= _end;
    }

    reference operator[](size_t i) const restrict
    {
        assert(_begin + i < _end);
        return _begin[i];
    }

    span<T> sub_range(size_t size, size_t start = 0) const restrict
    {
        assert(size + start < this->size());
        return span(_begin + start, _begin + start + size);
    }

    template<typename U,
             typename = decltype(reinterpret_cast<U*>(std::declval<T*>()))> // SFINAE tests for constness consistency
    span<U> reinterpret() const restrict
    {
        assert((size() * sizeof(T)) % sizeof(U) == 0);
        assert(reinterpret_cast<uintptr_t>(_begin) % alignof(U) == 0);
        return {reinterpret_cast<U*>(_begin), reinterpret_cast<U*>(_end)};
    }

  private:
    T* _begin;
    T* _end;
};

/// Like std::unique_ptr but is also a range, or like span but with owning semantics
template<typename T, typename D = std::default_delete<T[]>> class unique_span : private D
{
  public:
    using element_type    = T;
    using pointer         = T*;
    using const_pointer   = const T*;
    using iterator        = T*;
    using const_iterator  = const T*;
    using reference       = T&;
    using const_reference = const T&;

    unique_span() noexcept
      : D()
      , _begin(nullptr)
      , _end(nullptr)
    {}

    template<typename _D = D, typename = typename std::enable_if<std::is_nothrow_constructible<D, _D&&>::value>::type>
    unique_span(T* start, T* last, _D&& d = {}) noexcept
      : D(std::forward<_D>(d))
      , _begin(start)
      , _end(last)
    {
        assert(last >= start);
    }

    template<typename _D = D, typename = typename std::enable_if<std::is_nothrow_constructible<D, _D&&>::value>::type>
    unique_span(T* start, size_t n, _D&& d = {}) noexcept
      : D(std::forward<_D>(d))
      , _begin(start)
      , _end(start + n)
    {}

    template<typename R,
             typename _D = D,
             typename
             = typename std::enable_if<std::is_nothrow_constructible<D, _D&&>::value
                                       && std::is_convertible<decltype(std::declval<R>().begin()), T*>::value>::type>
    unique_span(R& r, _D&& d = {}) noexcept
      : D(std::forward<_D>(d))
      , _begin(r.begin())
      , _end(r.end())
    {
        assert(_end >= _begin);
    }

    unique_span(unique_span&& from) noexcept
      : D(std::move(from))
      , _begin(nullptr)
      , _end(nullptr)
    {
        std::swap(_begin, from._begin);
        std::swap(_end, from._end);
    }

    unique_span& operator=(unique_span&& from) noexcept
    {
        std::swap(_begin, from._begin);
        std::swap(_end, from._end);
        *static_cast<D*>(this) = static_cast<D&&>(from);
        return *this;
    }

    unique_span(const unique_span&) = delete;
    unique_span& operator=(const unique_span&) = delete;

    ~unique_span()
    {
        if (_begin != nullptr) {
            details::invoke_deleter(*static_cast<D*>(this), _begin, size());
            _begin = nullptr;
        }
    }

    size_t size() const { return size_t(_end - _begin); }
    bool   empty() const { return size() == 0; }
           operator bool() const { return size() != 0; }

    iterator       begin() { return _begin; }
    iterator       end() { return _end; }
    const_iterator begin() const { return _begin; }
    const_iterator end() const { return _end; }

    bool includes(const T* p) const { return p >= _begin && p < _end; }
    bool includes(const T* b, const T* e) const { return b >= _begin && e <= _end; }
    bool includes(span<const T> other) const { return other.begin() >= _begin && other.end() <= _end; }

    reference operator[](size_t i)
    {
        assert(_begin + i < _end);
        return _begin[i];
    }

    const_reference operator[](size_t i) const
    {
        assert(_begin + i < _end);
        return _begin[i];
    }

    span<T> sub_range(size_t size, size_t start = 0)
    {
        assert(size + start < this->size());
        return {_begin + start, _begin + start + size};
    }

    span<const T> sub_range(size_t size, size_t start = 0) const
    {
        assert(size + start < this->size());
        return {_begin + start, _begin + start + size};
    }

    template<typename U> span<const U> reinterpret() const
    {
        assert((size() * sizeof(T)) % sizeof(U) == 0);
        assert(reinterpret_cast<uintptr_t>(_begin) % alignof(U) == 0);
        return {reinterpret_cast<const U*>(_begin), reinterpret_cast<const U*>(_end)};
    }

    template<typename U> span<U> reinterpret()
    {
        assert((size() * sizeof(T)) % sizeof(U) == 0);
        assert(reinterpret_cast<uintptr_t>(_begin) % alignof(U) == 0);
        return {reinterpret_cast<U*>(_begin), reinterpret_cast<U*>(_end)};
    }

    bool operator==(const span<T>& other) const { return other._begin == this->_begin && other._end == this->_end; }

  protected:
    T* _begin;
    T* _end;
};

template<typename T>
inline unique_span<T>
make_unique_span(size_t n)
{
    return {new T[n], n};
}

template<typename T> struct free_deleter
{
    void operator()(T* ptr)
    {
        details::destruct(ptr);
        free(ptr);
    }
};

template<typename T> struct free_deleter<T[]>
{
    void operator()(T* ptr, size_t size)
    {
        details::destruct(ptr, ptr + size);
        free(ptr);
    }
};

template<typename T> using malloc_span = unique_span<T, free_deleter<T[]>>;

#ifdef __linux__
// Bet that if the libc doesn't support these, the kernel does. Otherwise it's just a nonop.
#    ifndef MADV_HUGEPAGE
#        define MADV_HUGEPAGE 15
#    endif

#    ifndef MADV_DONTDUMP
#        define MADV_DONTDUMP 16
#    endif

inline void
madvise_huge(void* addr, size_t len, bool dontdump = true)
{
    // In general, regions with huge page contains data that should not be in coredumps
    if (dontdump) madvise(addr, len, MADV_DONTDUMP);

    addr = details::round_up<details::huge_page_size>(addr);
    len  = details::round_down<details::huge_page_size>(len);

    auto res = ::madvise(addr, len, MADV_HUGEPAGE);
    static_cast<void>(res);
    PRINT_DEBUG("madvise(%p, 0x%lx, MADV_HUGEPAGE)=%d\n", addr, len, res);
}

#else
inline void
madvise_huge(void* p, size_t size)
{}
#endif

template<typename T>
inline malloc_span<T>
alloc_huge(size_t n)
{

    const size_t bytes = details::round_up<details::huge_page_size>(n * sizeof(T));

    T* ptr;
    if (posix_memalign(reinterpret_cast<void**>(&ptr), details::huge_page_size, bytes) != 0) throw std::bad_alloc();

    madvise_huge(ptr, bytes);
    return malloc_span<T>(ptr, bytes / sizeof(T));
}

template<typename T> struct mmap_deleter
{
    void operator()(T* ptr)
    {
        details::destruct(ptr);
        sys::check_ret(munmap(ptr, sizeof(T)), "munmap");
    }
};

template<typename T> struct mmap_deleter<T[]>
{
    void operator()(T* ptr, size_t size)
    {
        details::destruct(ptr, ptr + size);
        sys::check_ret(munmap(ptr, sizeof(T)), "munmap");
    }
};

template<typename T> using mmap_span = unique_span<T, mmap_deleter<T[]>>;

template<typename T>
inline mmap_span<T>
alloc_mirrored(size_t size, size_t ncopies, const char* shm_name_base = nullptr)
{
    const size_t nbytes = sizeof(T) * size;

    // Create a backing file in tmpfs
    int fd = details::tmpshm(shm_name_base);
    sys::check_ret(ftruncate(fd, ssize_t(nbytes)), "ftrunctate");

    // Allocate the virtual contiguous address space
    T* ptr = reinterpret_cast<T*>(sys::check_ptr(
      mmap(nullptr, ncopies * nbytes, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0), "mmap (anonymous)"));

    // Map the file ncopies times
    for (size_t i = 0; i < ncopies; i++) {
        T* img_ptr = ptr + i * size;
        if (mmap(img_ptr, nbytes, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_FIXED, fd, 0) != img_ptr)
            sys::throw_syserr("mmap (shared 0");
    }

    sys::check_ret(close(fd), "close");

    return {ptr, ncopies * size};
}

template<typename Lockable = std::mutex> struct lock_releaser : private std::unique_lock<Lockable>
{
    using std::unique_lock<Lockable>::unique_lock;
    lock_releaser() noexcept = default;

    lock_releaser(std::unique_lock<Lockable>&& lock) noexcept
      : std::unique_lock<Lockable>(std::move(lock))
    {}

    template<typename T> void operator()(T*)
    {
        if (this->owns_lock()) this->unlock();
    }
};

template<typename T, typename Lockable = std::mutex> using locked_span = unique_span<T, lock_releaser<Lockable>>;

#endif // MEMORY_HPP
