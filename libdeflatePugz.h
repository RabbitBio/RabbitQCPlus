/*
 * libdeflate.h - public header for libdeflate
 */

#ifndef LIBDEFLATE_H
#define LIBDEFLATE_H

#define LIBDEFLATE_VERSION_MAJOR 0
#define LIBDEFLATE_VERSION_MINOR 8
#define LIBDEFLATE_VERSION_STRING "0.8"

#include <cstddef>
#include <cstdint>

/* ========================================================================== */
/*                             Decompression                                  */
/* ========================================================================== */

/*
 * Result of a call to libdeflate_deflate_decompress(),
 * libdeflate_zlib_decompress(), or libdeflate_gzip_decompress().
 */
enum libdeflate_result {
    /* Decompression was successful.  */
    LIBDEFLATE_SUCCESS = 0,

    /* Decompressed failed because the compressed data was invalid, corrupt,
     * or otherwise unsupported.  */
    LIBDEFLATE_BAD_DATA = 1,

    /* A nullptr 'actual_out_nbytes_ret' was provided, but the data would have
     * decompressed to fewer than 'out_nbytes_avail' bytes.  */
    LIBDEFLATE_SHORT_OUTPUT = 2,

    /* The data would have decompressed to more than 'out_nbytes_avail'
     * bytes.  */
    LIBDEFLATE_INSUFFICIENT_SPACE = 3,
};

#endif /* LIBDEFLATE_H */
