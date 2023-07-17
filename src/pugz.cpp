//
// Created by ylf9811 on 2021/11/3.
//

#include "pugz.h"

#include "../lib/gzip_decompress.hpp"
#include "prog_util.h"

using namespace std;

static int
stat_file(struct file_stream *in, stat_t *stbuf, bool allow_hard_links) {
    if (tfstat(in->fd, stbuf) != 0) {
        msg("%" TS ": unable to stat file", in->name);
        return -1;
    }

    if (!S_ISREG(stbuf->st_mode) && !in->is_standard_stream) {
        msg("%" TS " is %s -- skipping", in->name, S_ISDIR(stbuf->st_mode) ? "a directory" : "not a regular file");
        return -2;
    }

    if (stbuf->st_nlink > 1 && !allow_hard_links) {
        msg("%" TS " has multiple hard links -- skipping "
            "(use -f to process anyway)",
            in->name);
        return -2;
    }

    return 0;
}


void main_pugz(string in_name, int threads, moodycamel::ReaderWriterQueue<pair<char *, int>> *Q, atomic_int *producerDone) {
    struct file_stream in;
    stat_t stbuf;
    int ret;
    const byte *in_p;

    ret = xopen_for_read(in_name.c_str(), true, &in);
    if (ret != 0) {
        fprintf(stderr, "gg on xopen_for_read\n");
        exit(0);
    }

    ret = stat_file(&in, &stbuf, true);
    if (ret != 0) {
        fprintf(stderr, "gg on stat_file\n");
        exit(0);
    }
    /* TODO: need a streaming-friendly solution */
    ret = map_file_contents(&in, size_t(stbuf.st_size));

    if (ret != 0) {
        fprintf(stderr, "gg on map_file_contents\n");
        exit(0);
    }

    in_p = static_cast<const byte *>(in.mmap_mem);
    OutputConsumer output{};

    output.P = Q;
    output.pDone = producerDone;
    ConsumerSync sync{};
    libdeflate_gzip_decompress(in_p, in.mmap_size, threads, output, &sync);
    //    xclose(&in);
}
