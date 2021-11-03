/*
 * prog_util.c - utility functions for programs
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
 */

#include "../lib/memory.hpp"
#include "prog_util.h"

#include <errno.h>
#include <fcntl.h>
#include <stdarg.h>
#include <time.h>
#ifdef _WIN32
#    include <windows.h>
#else
#    include <unistd.h>
#    include <sys/mman.h>
#    include <sys/time.h>
#endif

#ifndef O_BINARY
#    define O_BINARY 0
#endif
#ifndef O_SEQUENTIAL
#    define O_SEQUENTIAL 0
#endif
#ifndef O_NOFOLLOW
#    define O_NOFOLLOW 0
#endif
#ifndef O_NONBLOCK
#    define O_NONBLOCK 0
#endif
#ifndef O_NOCTTY
#    define O_NOCTTY 0
#endif

/* The invocation name of the program (filename component only) */
const tchar* _program_invocation_name;

static void
do_msg(const char* format, bool with_errno, va_list va)
{
    int saved_errno = errno;

    fprintf(stderr, "%" TS ": ", program_invocation_name);
    vfprintf(stderr, format, va);
    if (with_errno)
        fprintf(stderr, ": %s\n", strerror(saved_errno));
    else
        fprintf(stderr, "\n");

    errno = saved_errno;
}

/* Print a message to standard error */
void
msg(const char* format, ...)
{
    va_list va;

    va_start(va, format);
    do_msg(format, false, va);
    va_end(va);
}

/* Print a message to standard error, including a description of errno */
void
msg_errno(const char* format, ...)
{
    va_list va;

    va_start(va, format);
    do_msg(format, true, va);
    va_end(va);
}

/*
 * Return the number of timer ticks that have elapsed since some unspecified
 * point fixed at the start of program execution
 */
uint64_t
timer_ticks(void)
{
#ifdef _WIN32
    LARGE_INTEGER count;
    QueryPerformanceCounter(&count);
    return count.QuadPart;
#elif defined(HAVE_CLOCK_GETTIME)
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (1000000000 * (uint64_t)ts.tv_sec) + ts.tv_nsec;
#else
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return (1000000 * uint64_t(tv.tv_sec)) + uint64_t(tv.tv_usec);
#endif
}

/*
 * Return the number of timer ticks per second
 */
static uint64_t
timer_frequency(void)
{
#ifdef _WIN32
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    return freq.QuadPart;
#elif defined(HAVE_CLOCK_GETTIME)
    return 1000000000;
#else
    return 1000000;
#endif
}

/*
 * Convert a number of elapsed timer ticks to milliseconds
 */
uint64_t
timer_ticks_to_ms(uint64_t ticks)
{
    return ticks * 1000 / timer_frequency();
}

/*
 * Convert a byte count and a number of elapsed timer ticks to MB/s
 */
uint64_t
timer_MB_per_s(uint64_t bytes, uint64_t ticks)
{
    return bytes * timer_frequency() / ticks / 1000000;
}

/*
 * Retrieve a pointer to the filename component of the specified path.
 *
 * Note: this does not modify the path.  Therefore, it is not guaranteed to work
 * properly for directories, since a path to a directory might have trailing
 * slashes.
 */
const tchar*
get_filename(const tchar* path)
{
    const tchar* slash = tstrrchr(path, '/');
#ifdef _WIN32
    const tchar* backslash = tstrrchr(path, '\\');
    if (backslash != nullptr && (slash == nullptr || backslash > slash)) slash = backslash;
#endif
    if (slash != nullptr) return slash + 1;
    return path;
}

/* Create a copy of 'path' surrounded by double quotes */
static tchar*
quote_path(const tchar* path)
{
    size_t len    = tstrlen(path);
    tchar* result = new tchar[1 + len + 1 + 1];

    if (result == nullptr) return nullptr;
    result[0] = '"';
    tmemcpy(&result[1], path, len);
    result[1 + len]     = '"';
    result[1 + len + 1] = '\0';
    return result;
}

/* Open a file for reading, or set up standard input for reading */
int
xopen_for_read(const tchar* path, bool symlink_ok, struct file_stream* strm)
{
    strm->mmap_token = nullptr;
    strm->mmap_mem   = nullptr;

    if (path == nullptr) {
        strm->is_standard_stream = true;
        strm->name               = T("standard input");
        strm->fd                 = STDIN_FILENO;
#ifdef _WIN32
        _setmode(strm->fd, O_BINARY);
#endif
        return 0;
    }

    strm->is_standard_stream = false;

    strm->name = quote_path(path);
    if (strm->name == nullptr) return -1;

    strm->fd = topen(path, O_RDONLY | O_BINARY | O_NONBLOCK | O_NOCTTY | (symlink_ok ? 0 : O_NOFOLLOW) | O_SEQUENTIAL);
    if (strm->fd < 0) {
        msg_errno("Can't open %" TS " for reading", strm->name);
        delete strm->name;
        return -1;
    }

#if defined(HAVE_POSIX_FADVISE) && (O_SEQUENTIAL == 0)
    // posix_fadvise(strm->fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif

    return 0;
}

/* Open a file for writing, or set up standard output for writing */
int
xopen_for_write(const tchar* path, bool overwrite, struct file_stream* strm)
{
    int ret = -1;

    strm->mmap_token = nullptr;
    strm->mmap_mem   = nullptr;

    if (path == nullptr) {
        strm->is_standard_stream = true;
        strm->name               = T("standard output");
        strm->fd                 = STDOUT_FILENO;
#ifdef _WIN32
        _setmode(strm->fd, O_BINARY);
#endif
        return 0;
    }

    strm->is_standard_stream = false;

    strm->name = quote_path(path);
    if (strm->name == nullptr) goto err;
retry:
    strm->fd = topen(path, O_WRONLY | O_BINARY | O_NOFOLLOW | O_CREAT | O_EXCL, 0644);
    if (strm->fd < 0) {
        if (errno != EEXIST) {
            msg_errno("Can't open %" TS " for writing", strm->name);
            goto err;
        }
        if (!overwrite) {
            if (!isatty(STDERR_FILENO) || !isatty(STDIN_FILENO)) {
                msg("%" TS " already exists; use -f to overwrite", strm->name);
                ret = -2; /* warning only */
                goto err;
            }
            fprintf(stderr,
                    "%" TS ": %" TS " already exists; "
                    "overwrite? (y/n) ",
                    program_invocation_name,
                    strm->name);
            if (getchar() != 'y') {
                msg("Not overwriting.");
                goto err;
            }
        }
        if (tunlink(path) != 0) {
            msg_errno("Unable to delete %" TS, strm->name);
            goto err;
        }
        goto retry;
    }

    return 0;

err:
    delete strm->name;
    return ret;
}
#include <vector>
/* Read the full contents of a file into memory */
static int
read_full_contents(struct file_stream* strm)
{
    size_t            filled   = 0;
    size_t            capacity = 4096;
    std::vector<char> buf;
    ssize_t           ret;

    buf.resize(capacity);
    do {
        if (filled == capacity) {
            if (capacity == SIZE_MAX) goto oom;
            capacity += std::min(SIZE_MAX - capacity, capacity);
            buf.resize(capacity);
        }
        ret = xread(strm, &buf[filled], capacity - filled);
        if (ret < 0) goto err;
        filled += size_t(ret);
    } while (ret != 0);

    strm->mmap_mem  = &buf[0];
    strm->mmap_size = filled;
    return 0;

err:
    return int(ret);
oom:
    msg("Out of memory!  %" TS " is too large to be processed by "
        "this program as currently implemented.",
        strm->name);
    ret = -1;
    goto err;
}

/* Map the contents of a file into memory */
int
map_file_contents(struct file_stream* strm, uint64_t size)
{
    if (size == 0) /* mmap isn't supported on empty files */
        return read_full_contents(strm);

    if (size > SIZE_MAX) {
        msg("%" TS " is too large to be processed by this program", strm->name);
        return -1;
    }
#ifdef _WIN32
    strm->mmap_token
      = CreateFileMapping((HANDLE)(intptr_t)_get_osfhandle(strm->fd), nullptr, PAGE_READONLY, 0, 0, nullptr);
    if (strm->mmap_token == nullptr) {
        DWORD err = GetLastError();
        if (err == ERROR_BAD_EXE_FORMAT) /* mmap unsupported */
            return read_full_contents(strm);
        msg("Unable create file mapping for %" TS ": Windows error %u", strm->name, (unsigned int)err);
        return -1;
    }

    strm->mmap_mem = MapViewOfFile((HANDLE)strm->mmap_token, FILE_MAP_READ, 0, 0, size);
    if (strm->mmap_mem == nullptr) {
        msg("Unable to map %" TS " into memory: Windows error %u", strm->name, (unsigned int)GetLastError());
        CloseHandle((HANDLE)strm->mmap_token);
        return -1;
    }
#else /* _WIN32 */
    strm->mmap_mem = mmap(nullptr, size, PROT_READ, MAP_SHARED, strm->fd, 0);
    if (strm->mmap_mem == MAP_FAILED) {
        strm->mmap_mem = nullptr;
        if (errno == ENODEV) /* mmap isn't supported on this file */
            return read_full_contents(strm);
        if (errno == ENOMEM) {
            msg("%" TS " is too large to be processed by this "
                "program",
                strm->name);
        } else {
            msg_errno("Unable to map %" TS " into memory", strm->name);
        }
        return -1;
    }

    madvise_huge(strm->mmap_mem, size);

    strm->mmap_token = strm; /* anything that's not nullptr */

#endif /* !_WIN32 */
    strm->mmap_size = size;
    return 0;
}

/*
 * Read from a file, returning the full count to indicate all bytes were read, a
 * short count (possibly 0) to indicate EOF, or -1 to indicate error.
 */
ssize_t
xread(struct file_stream* strm, void* buf, size_t count)
{
    char*  p          = static_cast<char*>(buf);
    size_t orig_count = count;

    while (count != 0) {
        ssize_t res = read(strm->fd, p, std::min(count, size_t(INT_MAX)));
        if (res == 0) break;
        if (res < 0) {
            if (errno == EAGAIN || errno == EINTR) continue;
            msg_errno("Error reading from %" TS, strm->name);
            return -1;
        }
        p += size_t(res);
        count -= size_t(res);
    }
    return ssize_t(orig_count) - ssize_t(count);
}

/* Write to a file, returning 0 if all bytes were written or -1 on error */
int
full_write(struct file_stream* strm, const void* buf, size_t count)
{
    const char* p = static_cast<const char*>(buf);

    while (count != 0) {
        ssize_t res = write(strm->fd, p, std::min(count, size_t(INT_MAX)));
        if (res <= 0) {
            msg_errno("Error writing to %" TS, strm->name);
            return -1;
        }
        p += size_t(res);
        count -= size_t(res);
    }
    return 0;
}

/* Close a file, returning 0 on success or -1 on error */
int
xclose(struct file_stream* strm)
{
    int ret = 0;

    if (!strm->is_standard_stream) {
        if (close(strm->fd) != 0) {
            msg_errno("Error closing %" TS, strm->name);
            ret = -1;
        }
        delete[] strm->name;
    }

    if (strm->mmap_token != nullptr) {
#ifdef _WIN32
        UnmapViewOfFile(strm->mmap_mem);
        CloseHandle((HANDLE)strm->mmap_token);
#else
        munmap(strm->mmap_mem, strm->mmap_size);
#endif
        strm->mmap_token = nullptr;
    } else {
        free(strm->mmap_mem);
    }
    strm->mmap_mem = nullptr;
    strm->fd       = -1;
    strm->name     = nullptr;
    return ret;
}
