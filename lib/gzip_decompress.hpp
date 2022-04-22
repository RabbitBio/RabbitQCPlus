/*
 * gzip_decompress.c - decompress with a gzip wrapper
 *
 * Originally public domain; changes after 2016-09-07 are copyrighted.
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

#include "gzip_constants.h"

#include "lib/libdeflatePugz.h"
#include <exception>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>

#include "deflate_decompress.hpp" //FIXME

template<typename Consumer>
static enum libdeflate_result
libdeflate_gzip_decompress(const byte *in, size_t in_nbytes, unsigned nthreads, Consumer &consumer,
                           ConsumerSync *sync) {
    // FIXME: handle header parsing inside DeflateThread*, allowing multimember gzip files
    InputStream in_stream2(in, in_nbytes);
    in_stream2.consume_header();
    InputStream in_stream(in_stream2.in_next, in_stream2.available());
    size_t in_size = in_stream2.available();
    //TODO GG when gz file is not big
    //nthreads = std::min(1 + unsigned(in_size >> 25), nthreads);
    nthreads = std::min(1 + unsigned(in_size >> 21), nthreads);

    PRINT_DEBUG("Using %u threads\n", nthreads);

    std::vector<std::thread> threads;
    std::vector<DeflateThread *> deflate_threads(nthreads);

    std::atomic<size_t> nready = {0};
    //std::vector<std::atomic<size_t>> nreadys;

    std::condition_variable ready;
    std::mutex ready_mtx;
    std::exception_ptr exception;

    threads.reserve(nthreads);

    // Sections of file decompressed sequentially
    size_t max_section_size = nthreads * (32ull << 20); // 32MB per thread
    size_t section_size = std::min(max_section_size, in_size);
    size_t n_sections = (in_size + section_size - 1) / section_size;
    if (nthreads == 1) n_sections = 1;
    section_size = in_size / n_sections;

    // Section are chunked to nethreads
    size_t chunk_size = section_size / nthreads;
    // The first thread is working with resolved context so its faster
    size_t first_chunk_size = chunk_size + (4UL << 20); // FIXME: ratio instead of delta
    if (nthreads > 1) {
        chunk_size = (nthreads * chunk_size - first_chunk_size) / (nthreads - 1);
        first_chunk_size = section_size - chunk_size * (nthreads - 1);
    }
    PRINT_DEBUG("in_size %zu\n", in_size);
    PRINT_DEBUG("n_sections %zu\n", n_sections);
    PRINT_DEBUG("section_size %zu\n", section_size);
    PRINT_DEBUG("nthreads %zu\n", nthreads);
    PRINT_DEBUG("chunk_size %zu\n", chunk_size);
    PRINT_DEBUG("first_chunk_size %zu\n", first_chunk_size);
    //nreadys.reserve(n_sections);
    std::atomic<size_t> *nreadys=new std::atomic<size_t>[n_sections];
    for (unsigned i = 0; i < n_sections; i++) {
        nreadys[i] = 0;
    }

    for (unsigned chunk_idx = 0; chunk_idx < nthreads; chunk_idx++) {
        if (chunk_idx == 0) {
            threads.emplace_back([&, chunk_idx]() {
                ConsumerWrapper<Consumer> consumer_wrapper{consumer, sync};
                consumer_wrapper.set_chunk_idx(0, nthreads == 1);
                DeflateThread deflate_thread(in_stream, consumer_wrapper);
                deflate_thread.threadNum = chunk_idx;

                PRINT_DEBUG("chunk 0 is %p\n", (void *) &deflate_thread);
                {
                    std::unique_lock<std::mutex> lock{ready_mtx};
                    deflate_threads[0] = &deflate_thread;
                    nready++;
                    ready.notify_all();
                    while (nready != nthreads)
                        ready.wait(lock);
                    PRINT_DEBUG("thread %lu init sync\n", chunk_idx);
                    //! thread sync
                    if (chunk_idx < nthreads - 1)
                        deflate_thread.set_downstream(deflate_threads[chunk_idx + 1]);
                }
                //usleep(100);

                try {
                    {
                        std::unique_lock<std::mutex> lock{ready_mtx};
                        nreadys[0]++;
                        ready.notify_all();
                        while (nreadys[0] != nthreads)
                            ready.wait(lock);
                        PRINT_DEBUG("thread %lu sync on section %lu\n", chunk_idx, 0);
                    }
                    //usleep(100);
                    const byte *last_unmapped = details::round_up<details::huge_page_size>(in);

                    // First chunk of first section: no context needed
                    if (chunk_idx == nthreads - 1)
                        deflate_thread.set_end_block(first_chunk_size * 8);
                    deflate_thread.go(0);

                    // First chunks of next sections get their contexts from the last chunk of the previous section
                    auto &prev_chunk = *deflate_threads[nthreads - 1];
                    for (unsigned section_idx = 1; section_idx < n_sections && nready == nthreads; section_idx++) {
                        {
                            std::unique_lock<std::mutex> lock{ready_mtx};
                            nreadys[section_idx]++;
                            ready.notify_all();
                            while (nreadys[section_idx] != nthreads)
                                ready.wait(lock);
                            PRINT_DEBUG("thread %lu sync on section %lu\n", chunk_idx, section_idx);
                        }
                        const size_t start = section_idx * section_size;
                        const size_t stop = start + first_chunk_size;

                        // Get the context and position of the first block of the section
                        size_t resume_bitpos;
                        { // Synchronization point
                            auto ctx = prev_chunk.get_context();
                            deflate_thread.set_initial_context(ctx.first);
                            resume_bitpos = ctx.second;
                        }

                        //! ?
                        // Unmmap the section previously decompressed (free RSS, usefull for large files)
                        const byte *unmap_end
                                = details::round_down<details::huge_page_size>(
                                        in_stream.data.begin() + resume_bitpos / 8);
                        sys::check_ret(munmap(const_cast<byte *>(last_unmapped), size_t(unmap_end - last_unmapped)),
                                       "munmap");
                        last_unmapped = unmap_end;
                        //! ?


                        PRINT_DEBUG("=====%p chunk 0 of section %u: [%lu<=%lu, %lu]\n",
                                    (void *) &deflate_thread,
                                    section_idx,
                                    start * 8,
                                    resume_bitpos,
                                    stop * 8);

                        assert(resume_bitpos >= start * 8);
                        assert(resume_bitpos < stop * 8);
                        //! ylf fix0
                        //! all add a condition before set block end, actually in sync(),
                        //! all stream a real block end from its down stream, except the last one
                        //! and, if set here, maybe overflow the real block value then get error
                        if (chunk_idx == nthreads - 1)
                            deflate_thread.set_end_block(stop * 8);
                        consumer_wrapper.set_section_idx(section_idx);
                        deflate_thread.go(resume_bitpos);
                    }

                } catch (...) {
                    std::unique_lock<std::mutex> lock{ready_mtx};
                    if (!exception) {
                        exception = std::current_exception();
                        nready = 0; // Stop the thread pool
                        for (unsigned i = 0; i < n_sections; i++)
                            nreadys[i] = 0;
                    }
                    return;
                }

                if (nthreads == 1) deflate_thread.get_context();
            });
        } else {
            threads.emplace_back([&, chunk_idx]() {
                ConsumerWrapper<Consumer> consumer_wrapper{consumer, sync};
                consumer_wrapper.set_chunk_idx(chunk_idx, chunk_idx == nthreads - 1);
                DeflateThreadRandomAccess deflate_thread{in_stream, consumer_wrapper};
                deflate_thread.threadNum = chunk_idx;
                PRINT_DEBUG("=====chunk %u is %p\n", chunk_idx, (void *) &deflate_thread);
                {
                    //! init all deflate_thread
                    std::unique_lock<std::mutex> lock{ready_mtx};
                    deflate_threads[chunk_idx] = &deflate_thread;
                    nready++;
                    ready.notify_all();

                    while (nready != nthreads)
                        ready.wait(lock);
                    //! thread sync
                    PRINT_DEBUG("thread %lu init sync\n", chunk_idx);
                    if (chunk_idx < nthreads - 1)
                        deflate_thread.set_downstream(deflate_threads[chunk_idx + 1]);
                    deflate_thread.set_upstream(deflate_threads[chunk_idx - 1]);
                }
                //usleep(100);

                try {
                    const size_t chunk_offset_start = first_chunk_size + chunk_size * (chunk_idx - 1);
                    const size_t chunk_offset_stop = first_chunk_size + chunk_size * chunk_idx;


                    for (unsigned section_idx = 0; section_idx < n_sections && nready == nthreads; section_idx++) {
                        {
                            std::unique_lock<std::mutex> lock{ready_mtx};
                            nreadys[section_idx]++;
                            ready.notify_all();
                            while (nreadys[section_idx] != nthreads)
                                ready.wait(lock);
                            PRINT_DEBUG("thread %lu sync on section %lu\n", chunk_idx, section_idx);
                        }
                        //usleep(100);
                        const size_t section_offset = section_idx * section_size;
                        const size_t start = section_offset + chunk_offset_start;
                        const size_t stop = section_offset + chunk_offset_stop;

                        //! roughly divided by size
                        PRINT_DEBUG("=====%p chunk %u of section %u [%lu, %lu]\n",
                                    (void *) &deflate_thread,
                                    chunk_idx,
                                    section_idx,
                                    start * 8,
                                    stop * 8);
                        //! set block end roughly
                        if (chunk_idx == nthreads - 1)
                            deflate_thread.set_end_block(stop * 8);

                        consumer_wrapper.set_section_idx(section_idx);
                        //! main process
                        if (!deflate_thread.go(start * 8)) return;
                        assert(chunk_idx != nthreads - 1 || stop == section_offset + section_size);
                    }

                    // No one will need the context from the last chunk of the last section
                    if (chunk_idx == nthreads - 1)
                        // Releasing it will enable the imediate destruction of the DeflateThread
                        deflate_thread.get_context();
                } catch (...) {
                    std::unique_lock<std::mutex> lock{ready_mtx};
                    if (!exception) {
                        exception = std::current_exception();
                        nready = 0; // Stop the thread pool
                        for (unsigned i = 0; i < n_sections; i++)
                            nreadys[i] = 0;
                    }
                    return;
                }
            });
        }
    }

    for (auto &thread: threads)
        thread.join();

    if (exception) { std::rethrow_exception(exception); }

    return LIBDEFLATE_SUCCESS;
}
