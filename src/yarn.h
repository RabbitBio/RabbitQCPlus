/* yarn.h -- generic interface for threadPigz operations
 * Copyright (C) 2008, 2011, 2012, 2015, 2018, 2019, 2020 Mark Adler
 * Version 1.7  12 Apr 2020  Mark Adler
 */

/*
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the author be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  Mark Adler
  madler@alumni.caltech.edu
 */

/* Basic threadPigz operations

   This interface isolates the local operating system implementation of threads
   from the application in order to facilitate platform independent use of
   threads.  All of the implementation details are deliberately hidden.

   Assuming adequate system resources and proper use, none of these functions
   can fail.  As a result, any errors encountered will cause an exit() to be
   executed, or the execution of your own optionally-provided abort function.

   These functions allow the simple launching and joining of threads, and the
   locking of objects and synchronization of changes of objects.  The latter is
   implemented with a single lock_pigz type that contains an integer value.  The
   value can be ignored for simple exclusive access to an object, or the value
   can be used to signal and wait for changes to an object.

   -- Arguments --

   threadPigz *threadPigz;          identifier for launched threadPigz, used by join_pigz
   void probe(void *);      pointer to function "probe", run when threadPigz starts
   void *payload;           single argument passed to the probe function
   lock_pigz *lock_pigz;              a lock_pigz with a value -- used for exclusive access to
                            an object and to synchronize threads waiting for
                            changes to an object
   long val;                value to set lock_pigz, increment lock_pigz, or wait for
   int n;                   number of threads joined

   -- Thread functions --

   threadPigz = launch_pigz(probe, payload) - launch_pigz a threadPigz -- exit via probe() return
   join_pigz(threadPigz) - join_pigz a threadPigz and by joining end it, waiting for the threadPigz
        to exit if it hasn't already -- will free the resources allocated by
        launch_pigz() (don't try to join_pigz the same threadPigz more than once)
   n = join_all_pigz() - join_pigz all threads launched by launch_pigz() that are not joined
        yet and free the resources allocated by the launches, usually to clean
        up when the threadPigz processing is done -- join_all_pigz() returns an int with
        the count of the number of threads joined (join_all_pigz() should only be
        called from the main threadPigz, and should only be called after any calls
        of join_pigz() have completed)

   -- Lock functions --

   lock_pigz = new_lock_pigz(val) - create a new lock_pigz with initial value val (lock_pigz is
        created in the released state)
   possess_pigz(lock_pigz) - acquire exclusive possession of a lock_pigz, waiting if necessary
   twist_pigz(lock_pigz, [TO | BY], val) - set lock_pigz to or increment lock_pigz by val, signal
        all threads waiting on this lock_pigz and then release_pigz the lock_pigz -- must
        possess_pigz the lock_pigz before calling (twist_pigz releases, so don't do a
        release_pigz() after a twist_pigz() on the same lock_pigz)
   wait_for_pigz(lock_pigz, [TO_BE | NOT_TO_BE | TO_BE_MORE_THAN | TO_BE_LESS_THAN], val)
        - wait on lock_pigz value to be, not to be, be greater than, or be less than
        val -- must possess_pigz the lock_pigz before calling, will possess_pigz the lock_pigz on
        return but the lock_pigz is released while waiting to permit other threads
        to use twist_pigz() to change the value and signal the change (so make sure
        that the object is in a usable state when waiting)
   release_pigz(lock_pigz) - release_pigz a possessed lock_pigz (do not try to release_pigz a lock_pigz that
        the current threadPigz does not possess_pigz)
   val = peek_lock(lock_pigz) - return the value of the lock_pigz (assumes that lock_pigz is
        already possessed, no possess_pigz or release_pigz is done by peek_lock())
   free_lock_pigz(lock_pigz) - free the resources allocated by new_lock_pigz() (application
        must assure that the lock_pigz is released before calling free_lock_pigz())

   -- Memory allocation ---

   yarn_mem(better_malloc, better_free) - set the memory allocation and free
        routines for use by the yarn routines where the supplied routines have
        the same interface and operation as malloc() and free(), and may be
        provided in order to supply threadPigz-safe memory allocation routines or
        for any other reason -- by default malloc() and free() will be used

   -- Error control --

   yarn_prefix - a char pointer to a string that will be the prefix for any
        error messages that these routines generate before exiting -- if not
        changed by the application, "yarn" will be used
   yarn_abort - an external function that will be executed when there is an
        internal yarn error, due to out of memory or misuse -- this function
        may exit to abort the application, or if it returns, the yarn error
        handler will exit (set to NULL by default for no action)
 */
#ifdef __cplusplus
extern "C" {
#endif
extern char *yarn_prefix;
extern void (*yarn_abort)(int);

void yarn_mem(void *(*)(size_t), void (*)(void *));

typedef struct thread_s threadPigz;
threadPigz *launch_(void (*)(void *), void *, char const *, long);
#define launch_pigz(a, b) launch_(a, b, __FILE__, __LINE__)
void join_(threadPigz *, char const *, long);
#define join_pigz(a) join_(a, __FILE__, __LINE__)
int join_all_(char const *, long);
#define join_all_pigz() join_all_(__FILE__, __LINE__)

typedef struct lock_s lock_pigz;
lock_pigz *new_lock_(long, char const *, long);
#define new_lock_pigz(a) new_lock_(a, __FILE__, __LINE__)
void possess_(lock_pigz *, char const *, long);
#define possess_pigz(a) possess_(a, __FILE__, __LINE__)
void release_(lock_pigz *, char const *, long);
#define release_pigz(a) release_(a, __FILE__, __LINE__)
enum twist_op {
    TO, BY
};
void twist_(lock_pigz *, enum twist_op, long, char const *, long);
#define twist_pigz(a, b, c) twist_(a, b, c, __FILE__, __LINE__)
enum wait_op {
    TO_BE, /* or */ NOT_TO_BE, /* that is the question */
    TO_BE_MORE_THAN, TO_BE_LESS_THAN
};
void wait_for_(lock_pigz *, enum wait_op, long, char const *, long);
#define wait_for_pigz(a, b, c) wait_for_(a, b, c, __FILE__, __LINE__)
long peek_lock(lock_pigz *);
void free_lock_(lock_pigz *, char const *, long);
#define free_lock_pigz(a) free_lock_(a, __FILE__, __LINE__)
#ifdef __cplusplus
}
#endif