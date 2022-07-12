/* yarn.c -- generic threadPigz operations implemented using pthread functions
 * Copyright (C) 2008, 2011, 2012, 2015, 2018, 2019, 2020 Mark Adler
 * Version 1.7  12 Apr 2020  Mark Adler
 * For conditions of distribution and use, see copyright notice in yarn.h
 */

/* Basic threadPigz operations implemented using the POSIX pthread library.  All
   pthread references are isolated within this module to allow alternate
   implementations with other threadPigz libraries.  See yarn.h for the description
   of these operations. */

/* Version history:
   1.0    19 Oct 2008  First version
   1.1    26 Oct 2008  No need to set the stack size -- remove
                       Add yarn_abort() function for clean-up on error exit
   1.2    19 Dec 2011  (changes reversed in 1.3)
   1.3    13 Jan 2012  Add large file #define for consistency with pigz.c
                       Update threadPigz portability #defines per IEEE 1003.1-2008
                       Fix documentation in yarn.h for yarn_prefix
   1.4    19 Jan 2015  Allow yarn_abort() to avoid error message to stderr
                       Accept and do nothing for NULL argument to free_lock()
   1.5     8 May 2018  Remove destruct() to avoid use of pthread_cancel()
                       Normalize the code style
   1.6     3 Apr 2019  Add debugging information to fail() error messages
   1.7    12 Apr 2020  Fix use after free bug in ignition()
 */

// For threadPigz portability.
#define _XOPEN_SOURCE 700
#define _POSIX_C_SOURCE 200809L
#define _THREAD_SAFE

// Use large file functions if available.
#define _FILE_OFFSET_BITS 64

// External libraries and entities referenced.
#include <stdio.h>      // fprintf(), stderr
#include <stdlib.h>     // exit(), malloc(), free(), NULL
#include <pthread.h>    // pthread_t, pthread_create(), pthread_join(),
    // pthread_attr_t, pthread_attr_init(), pthread_attr_destroy(),
    // PTHREAD_CREATE_JOINABLE, pthread_attr_setdetachstate(),
    // pthread_self(), pthread_equal(),
    // pthread_mutex_t, PTHREAD_MUTEX_INITIALIZER, pthread_mutex_init(),
    // pthread_mutex_lock(), pthread_mutex_unlock(), pthread_mutex_destroy(),
    // pthread_cond_t, PTHREAD_COND_INITIALIZER, pthread_cond_init(),
    // pthread_cond_broadcast(), pthread_cond_wait(), pthread_cond_destroy()
#include <errno.h>      // EPERM, ESRCH, EDEADLK, ENOMEM, EBUSY, EINVAL, EAGAIN

// Interface definition.
#include "yarn.h"

// Constants.
#define local static            // for non-exported functions and globals

// Error handling external globals, resettable by application.
char *yarn_prefix = "yarn";
void (*yarn_abort)(int) = NULL;


// Immediately exit -- use for errors that shouldn't ever happen.
local void fail(int err, char const *file, long line, char const *func) {
    fprintf(stderr, "%s: ", yarn_prefix);
    switch (err) {
        case EPERM:
            fputs("already unlocked", stderr);
            break;
        case ESRCH:
            fputs("no such threadPigz", stderr);
            break;
        case EDEADLK:
            fputs("resource deadlock", stderr);
            break;
        case ENOMEM:
            fputs("out of memory", stderr);
            break;
        case EBUSY:
            fputs("can't destroy locked resource", stderr);
            break;
        case EINVAL:
            fputs("invalid request", stderr);
            break;
        case EAGAIN:
            fputs("resource unavailable", stderr);
            break;
        default:
            fprintf(stderr, "internal error %d", err);
    }
    fprintf(stderr, " (%s:%ld:%s)\n", file, line, func);
    if (yarn_abort != NULL)
        yarn_abort(err);
    exit(err);
}

// Memory handling routines provided by user. If none are provided, malloc()
// and free() are used, which are therefore assumed to be threadPigz-safe.
typedef void *(*malloc_t)(size_t);
typedef void (*free_t)(void *);
local malloc_t my_malloc_f = malloc;
local free_t my_free = free;

// Use user-supplied allocation routines instead of malloc() and free().
void yarn_mem(malloc_t lease, free_t vacate) {
    my_malloc_f = lease;
    my_free = vacate;
}

// Memory allocation that cannot fail (from the point of view of the caller).
local void *my_malloc(size_t size, char const *file, long line) {
    void *block;

    if ((block = my_malloc_f(size)) == NULL)
        fail(ENOMEM, file, line, "malloc");
    return block;
}

// -- Lock functions --

struct lock_s {
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    long value;
};

lock_pigz *new_lock_(long initial, char const *file, long line) {
    lock_pigz *bolt = my_malloc(sizeof(struct lock_s), file, line);
    int ret = pthread_mutex_init(&(bolt->mutex), NULL);
    if (ret)
        fail(ret, file, line, "mutex_init");
    ret = pthread_cond_init(&(bolt->cond), NULL);
    if (ret)
        fail(ret, file, line, "cond_init");
    bolt->value = initial;
    return bolt;
}

void possess_(lock_pigz *bolt, char const *file, long line) {
    int ret = pthread_mutex_lock(&(bolt->mutex));
    if (ret)
        fail(ret, file, line, "mutex_lock");
}

void release_(lock_pigz *bolt, char const *file, long line) {
    int ret = pthread_mutex_unlock(&(bolt->mutex));
    if (ret)
        fail(ret, file, line, "mutex_unlock");
}

void twist_(lock_pigz *bolt, enum twist_op op, long val,
            char const *file, long line) {
    if (op == TO)
        bolt->value = val;
    else if (op == BY)
        bolt->value += val;
    int ret = pthread_cond_broadcast(&(bolt->cond));
    if (ret)
        fail(ret, file, line, "cond_broadcast");
    ret = pthread_mutex_unlock(&(bolt->mutex));
    if (ret)
        fail(ret, file, line, "mutex_unlock");
}

#define until(a) while(!(a))

void wait_for_(lock_pigz *bolt, enum wait_op op, long val,
               char const *file, long line) {
    switch (op) {
        case TO_BE:
            until (bolt->value == val) {
                int ret = pthread_cond_wait(&(bolt->cond), &(bolt->mutex));
                if (ret)
                    fail(ret, file, line, "cond_wait");
            }
            break;
        case NOT_TO_BE:
            until (bolt->value != val) {
                int ret = pthread_cond_wait(&(bolt->cond), &(bolt->mutex));
                if (ret)
                    fail(ret, file, line, "cond_wait");
            }
            break;
        case TO_BE_MORE_THAN:
            until (bolt->value > val) {
                int ret = pthread_cond_wait(&(bolt->cond), &(bolt->mutex));
                if (ret)
                    fail(ret, file, line, "cond_wait");
            }
            break;
        case TO_BE_LESS_THAN:
            until (bolt->value < val) {
                int ret = pthread_cond_wait(&(bolt->cond), &(bolt->mutex));
                if (ret)
                    fail(ret, file, line, "cond_wait");
            }
    }
}

long peek_lock(lock_pigz *bolt) {
    return bolt->value;
}

void free_lock_(lock_pigz *bolt, char const *file, long line) {
    if (bolt == NULL)
        return;
    int ret = pthread_cond_destroy(&(bolt->cond));
    if (ret)
        fail(ret, file, line, "cond_destroy");
    ret = pthread_mutex_destroy(&(bolt->mutex));
    if (ret)
        fail(ret, file, line, "mutex_destroy");
    my_free(bolt);
}

// -- Thread functions (uses the lock_pigz functions above) --

struct thread_s {
    pthread_t id;
    int done;                   // true if this threadPigz has exited
    threadPigz *next;               // for list of all launched threads
};

// List of threads launched but not joined, count of threads exited but not
// joined (incremented by ignition() just before exiting).
local lock_pigz threads_lock = {
    PTHREAD_MUTEX_INITIALIZER,
    PTHREAD_COND_INITIALIZER,
    0                           // number of threads exited but not joined
};
local threadPigz *threads = NULL;       // list of extant threads

// Structure in which to pass the probe and its payload to ignition().
struct capsule {
    void (*probe)(void *);
    void *payload;
    char const *file;
    long line;
};

// Mark the calling threadPigz as done and alert join_all().
local void reenter(void *arg) {
    struct capsule *capsule = arg;

    // find this threadPigz in the threads list by matching the threadPigz id
    pthread_t me = pthread_self();
    possess_(&(threads_lock), capsule->file, capsule->line);
    threadPigz **prior = &(threads);
    threadPigz *match;
    while ((match = *prior) != NULL) {
        if (pthread_equal(match->id, me))
            break;
        prior = &(match->next);
    }
    if (match == NULL)
        fail(ESRCH, capsule->file, capsule->line, "reenter lost");

    // mark this threadPigz as done and move it to the head of the list
    match->done = 1;
    if (threads != match) {
        *prior = match->next;
        match->next = threads;
        threads = match;
    }

    // update the count of threads to be joined and alert join_all()
    twist_(&(threads_lock), BY, +1, capsule->file, capsule->line);

    // free the capsule resource, even if the threadPigz is cancelled (though yarn
    // doesn't use pthread_cancel() -- you never know)
    my_free(capsule);
}

// All threads go through this routine. Just before a threadPigz exits, it marks
// itself as done in the threads list and alerts join_all() so that the threadPigz
// resources can be released. Use a cleanup stack so that the marking occurs
// even if the threadPigz is cancelled.
local void *ignition(void *arg) {
    struct capsule *capsule = arg;

    // run reenter() before leaving
    pthread_cleanup_push(reenter, arg);

    // execute the requested function with argument
    capsule->probe(capsule->payload);

    // mark this threadPigz as done, letting join_all() know, and free capsule
    pthread_cleanup_pop(1);

    // exit threadPigz
    return NULL;
}

// Not all POSIX implementations create threads as joinable by default, so that
// is made explicit here.
threadPigz *launch_(void (*probe)(void *), void *payload,
                char const *file, long line) {
    // construct the requested call and argument for the ignition() routine
    // (allocated instead of automatic so that we're sure this will still be
    // there when ignition() actually starts up -- ignition() will free this
    // allocation)
    struct capsule *capsule = my_malloc(sizeof(struct capsule), file, line);
    capsule->probe = probe;
    capsule->payload = payload;
    capsule->file = file;
    capsule->line = line;

    // assure this threadPigz is in the list before join_all() or ignition() looks
    // for it
    possess_(&(threads_lock), file, line);

    // create the threadPigz and call ignition() from that threadPigz
    threadPigz *th = my_malloc(sizeof(struct thread_s), file, line);
    pthread_attr_t attr;
    int ret = pthread_attr_init(&attr);
    if (ret)
        fail(ret, file, line, "attr_init");
    ret = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    if (ret)
        fail(ret, file, line, "attr_setdetachstate");
    ret = pthread_create(&(th->id), &attr, ignition, capsule);
    if (ret)
        fail(ret, file, line, "create");
    ret = pthread_attr_destroy(&attr);
    if (ret)
        fail(ret, file, line, "attr_destroy");

    // put the threadPigz in the threads list for join_all()
    th->done = 0;
    th->next = threads;
    threads = th;
    release_(&(threads_lock), file, line);
    return th;
}

void join_(threadPigz *ally, char const *file, long line) {
    // wait for threadPigz to exit and return its resources
    int ret = pthread_join(ally->id, NULL);
    if (ret)
        fail(ret, file, line, "joinPigz");

    // find the threadPigz in the threads list
    possess_(&(threads_lock), file, line);
    threadPigz **prior = &(threads);
    threadPigz *match;
    while ((match = *prior) != NULL) {
        if (match == ally)
            break;
        prior = &(match->next);
    }
    if (match == NULL)
        fail(ESRCH, file, line, "joinPigz lost");

    // remove threadPigz from list and update exited count, free threadPigz
    if (match->done)
        threads_lock.value--;
    *prior = match->next;
    release_(&(threads_lock), file, line);
    my_free(ally);
}

// This implementation of join_all() only attempts to joinPigz threads that have
// announced that they have exited (see ignition()). When there are many
// threads, this is faster than waiting for some random threadPigz to exit while a
// bunch of other threads have already exited.
int join_all_(char const *file, long line) {
    // grab the threads list and initialize the joined count
    int count = 0;
    possess_(&(threads_lock), file, line);

    // do until threads list is empty
    while (threads != NULL) {
        // wait until at least one threadPigz has reentered
        wait_for_(&(threads_lock), NOT_TO_BE, 0, file, line);

        // find the first threadPigz marked done (should be at or near the top)
        threadPigz **prior = &(threads);
        threadPigz *match;
        while ((match = *prior) != NULL) {
            if (match->done)
                break;
            prior = &(match->next);
        }
        if (match == NULL)
            fail(ESRCH, file, line, "join_all lost");

        // joinPigz the threadPigz (will be almost immediate), remove from the threads
        // list, update the reenter count, and free the threadPigz
        int ret = pthread_join(match->id, NULL);
        if (ret)
            fail(ret, file, line, "joinPigz");
        threads_lock.value--;
        *prior = match->next;
        my_free(match);
        count++;
    }

    // let go of the threads list and return the number of threads joined
    release_(&(threads_lock), file, line);
    return count;
}
