/*****************************************************************************
 *
 * AM utilities
 *
 * released under MIT license
 *
 * 2008-2018 André Müller
 *
 *****************************************************************************/

#ifndef AM_PARALLEL_ECXECUTOR_H_
#define AM_PARALLEL_ECXECUTOR_H_

#include <vector>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <functional>
#include <thread>
#include <iterator>
#include <type_traits>
#include <tuple>

#include "tuple_apply.h" //because we don't have std::apply yet :-(


namespace am {


/*****************************************************************************
 *
 * @brief    runs a sequence of tasks in parallel
 *
 *           A task in this context is anything that is callable and that
 *           either returns void or the return value can be discarded.
 *
 * @details  Internally a number of threads is always kept alive and fired up
 *           as needed.
 *
 * @tparam   TaskIterator: InputIterator over tasks
 *           the task type itself must be callable, so either be a function
 *           pointer or a type that implements an operator() that can be
 *           invoked with CallArgs... e.g. std::function<void(CallArgs...)>
 *
 * @tparam   CallArgs...:
 *
 *****************************************************************************/
template<class TaskIterator, class... CallArgs>
class parallel_executor
{
    //used to store call arguments
    using args_type = std::tuple<CallArgs...>;

    enum class status { idle = 0, busy = 1, terminate = 2 };


public:
    //-----------------------------------------------------
    using iterator  = TaskIterator;
    using task_type = std::decay_t<decltype(*std::declval<iterator>())>;


    //---------------------------------------------------------------
    explicit
    parallel_executor(std::size_t concurrency = std::thread::hardware_concurrency()):
        mutables_{},
        isDone_{}, wakeUp_{},
        last_{}, first_{last_}, numUnfinished_{0},
        status_{status::idle},
        args_{}, workers_{}
    {
        make_workers(concurrency);
    }

    //-----------------------------------------------------
    parallel_executor(const parallel_executor&) = delete;
    parallel_executor(parallel_executor&& ) = delete;


    //---------------------------------------------------------------
    parallel_executor& operator = (const parallel_executor&) = delete;
    parallel_executor& operator = (parallel_executor&&) = delete;


    //-----------------------------------------------------
    ~parallel_executor() {
        abort();
    }


    //---------------------------------------------------------------
    std::size_t
    concurrency() const noexcept {
        return workers_.size();
    }


    //---------------------------------------------------------------
    /**
     * @brief runs a range of tasks and
     *        blocks execution of calling thread until all tasks are complete
     */
    template<class Range, class... Args, class T = int, class =
        std::enable_if_t<(sizeof...(Args) == sizeof...(CallArgs)),T>>
    void operator () (Range& range, Args&&... args)
    {
        using std::begin;
        using std::end;

        operator()(iterator(begin(range)), iterator(end(range)),
                   std::forward<Args>(args)...);
    }

    //-----------------------------------------------------
    /**
     * @brief runs a range of tasks and
     *        blocks execution of calling thread until all tasks are complete
     */
    template<class... Args, class T = int, class = typename
        std::enable_if_t<(sizeof...(Args) == sizeof...(CallArgs)),T>>
    void
    operator () (iterator first, iterator last, Args&&... args)
    {
        if(running() || (first == last)) return;

        std::unique_lock<std::mutex> lock(mutables_);

        first_ = first;
        last_ = last;
        using std::distance;
        numUnfinished_ = distance(first,last);

        //store call arguments
        args_ = args_type{std::forward<Args>(args)...};

        //fire up all worker threads
        status_.store(status::busy);
        wakeUp_.notify_all();

        //block execution of local thread until all tasks are complete
        while(running()) isDone_.wait(lock);
    }


    //---------------------------------------------------------------
    void
    stop() {
        if(status_ != status::busy) return;

        std::lock_guard<std::mutex> lock(mutables_);

        //mark all tasks as finished
        last_ = first_;

        //signal all threads to suspend after finishing their current task
        status_.store(status::idle);
        wakeUp_.notify_all();
    }


    //---------------------------------------------------------------
    bool
    running() const noexcept {
        return status_.load() == status::busy;
    }


private:

    //---------------------------------------------------------------
    void
    make_workers(std::size_t n)
    {
        if(n < 1) n = 1;

        workers_.reserve(n);

        //construct worker threads and suspend execution
        for(std::size_t i = 0; i < n; ++i) {
            workers_.emplace_back(std::thread( [&] {
                while(!abort_requested()) {
                    wait_until_needed();
                    if(!abort_requested())
                        try_run_next_task();
                }
            }));
        }
    }


    //---------------------------------------------------------------
    bool
    abort_requested() const noexcept {
        return (status_.load() == status::terminate);
    }

    //-----------------------------------------------------
    /// @brief this will render the object non-runnable
    void
    abort() {
        using std::next;

        {
            std::lock_guard<std::mutex> lock(mutables_);

            //mark all tasks as finished
            last_ = next(first_);
            numUnfinished_ = 0;

            //signal all threads to terminate immediatly
            status_.store(status::terminate);
            wakeUp_.notify_all();
        }

        //wait for all threads to be finished
        for(auto& w : workers_) w.join();
    }


    //---------------------------------------------------------------
    void
    wait_until_needed() {
        std::unique_lock<std::mutex> lock(mutables_);

        //block execution until woken up
        //(while loop protects against spurious wake-ups)
        while(status_.load() == status::idle || (first_ == last_)) {
            wakeUp_.wait(lock);
        }
    }

    //---------------------------------------------------------------
    /// @brief execute next unfinished task (or signal finished)
    void
    try_run_next_task() {
        iterator current;
        {
            std::lock_guard<std::mutex> lock(mutables_);

            if(first_ == last_) return;

            current = first_;
            ++first_;
        }

        //run task: current->operator()(args_...)
        apply(*current, args_);

        --numUnfinished_;
        //finished? => unblock everything
        if(numUnfinished_.load() < 1) {
            status_.store(status::idle);
            isDone_.notify_one();
        }
    }


    //---------------------------------------------------------------
    mutable std::mutex mutables_;
    std::condition_variable isDone_;
    std::condition_variable wakeUp_;
    iterator last_;
    iterator first_;
    std::atomic<int> numUnfinished_;
    std::atomic<status> status_;
    args_type args_;
    std::vector<std::thread> workers_;
};


}  // namespace am


#endif
