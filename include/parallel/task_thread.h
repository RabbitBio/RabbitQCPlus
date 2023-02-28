/*****************************************************************************
 *
 * AM utilities
 *
 * released under MIT license
 *
 * 2008-2018 André Müller
 *
 *****************************************************************************/

#ifndef AM_PARALLEL_TASK_THREAD_H_
#define AM_PARALLEL_TASK_THREAD_H_

#include <mutex>
#include <condition_variable>
#include <thread>
#include <atomic>


namespace am {


/*****************************************************************************
 *
 * @brief  thread that executes a single task and waits after completion
 *         until execution is demanded again
 *
 *****************************************************************************/
template<class Task>
class task_thread
{
    enum class status {busy = 0,  idle = 1, terminate = 2 };

public:
    //---------------------------------------------------------------
    using task_type = Task;


    //---------------------------------------------------------------
    explicit
    task_thread():
        mutables_{}, wakeup_{}, isIdle_{},
        status_{status::idle},
        task_{},
        thread_{ [&] { work(); }}
    {}

    //---------------------------------------------------------------
    template<
        class Arg, class... Args, class =
        typename std::enable_if<!std::is_same<typename std::decay<Arg>::type,task_thread>::value,Arg>::type
    >
    explicit
    task_thread(Arg&& arg, Args&&... args):
        mutables_{}, wakeup_{}, isIdle_{},
        status_{status::idle},
        task_{std::forward<Arg>(arg), std::forward<Args>(args)...},
        thread_{ [&] { work(); }}
    {}


    //---------------------------------------------------------------
    ~task_thread() {
        //signal thread to terminate
        status_.store(status::terminate);
        wakeup_.notify_all();
        isIdle_.notify_all();

        wait();

        thread_.join();
    }


    //---------------------------------------------------------------
    task_thread(const task_thread& source):
        mutables_{}, wakeup_{}, isIdle_{},
        status_{status::idle},
        task_{source.task_},
        thread_{ [&] { work(); }}
    {}

    //---------------------------------------------------------------
    task_thread(task_thread&& source)
        noexcept(noexcept(task_type(std::declval<task_type>())))
    :
        mutables_{}, wakeup_{}, isIdle_{},
        status_{status::idle},
        task_{std::move(source.task_)},
        thread_{ [&] { work(); }}
    {}


    //---------------------------------------------------------------
    task_thread& operator = (const task_thread& source)
    {
        wait();

        std::lock_guard<std::mutex> lock(mutables_);
        task_ = source.task_;

        return *this;
    }

    //---------------------------------------------------------------
    task_thread& operator = (task_thread&& source)
        noexcept(noexcept(std::declval<task_type>() = source.task_))
    {
        wait();

        std::lock_guard<std::mutex> lock(mutables_);
        task_ = std::move(source.task_);

        return *this;
    }


    //---------------------------------------------------------------
    bool available() const noexcept {
        return bool(status_.load());
    }


    //---------------------------------------------------------------
    /**
     * @brief block execution of calling thread until task is finished
     */
    void wait()
    {
        std::unique_lock<std::mutex> lock(mutables_);

        //(while loop protects against spurious wake-ups)
        while(status_.load() == status::busy) {
            isIdle_.wait(lock);
        }
    }


    //---------------------------------------------------------------
    /**
     * @brief non-thread-safe access to current task object
     */
    const task_type&
    unsafe_task() const & noexcept {
        return task_;
    }
    /**
     * @brief non-thread-safe access to current task object
     */
    task_type&
    unsafe_task() & noexcept {
        return task_;
    }
    /**
     * @brief extract current task object
     */
    task_type&&
    extract_task() && noexcept {
        return std::move(task_);
    }


    //---------------------------------------------------------------
    bool operator () ()
    {
        if(status_.load() != status::idle) return false;
        std::lock_guard<std::mutex> lock{mutables_};
        status_.store(status::busy);

        wakeup_.notify_all();
        return true;
    }

    //---------------------------------------------------------------
    bool operator () (const task_type& task)
    {
        if(status_.load() != status::idle) return false;
        std::lock_guard<std::mutex> lock{mutables_};
        status_.store(status::busy);

        task_ = task;

        wakeup_.notify_all();
        return true;
    }

    //---------------------------------------------------------------
    bool operator () (task_type&& task)
    {
        if(status_.load() != status::idle) return false;
        std::lock_guard<std::mutex> lock{mutables_};
        status_.store(status::busy);

        task_ = std::move(task);

        wakeup_.notify_all();
        return true;
    }

    //---------------------------------------------------------------
    template<class Arg, class... Args>
    bool operator () (Arg&& arg, Args&&... args)
    {
        if(status_.load() != status::idle) return false;
        std::lock_guard<std::mutex> lock{mutables_};
        status_.store(status::busy);

        task_ = task_type{std::forward<Arg>(arg),
                          std::forward<Args>(args)...};

        wakeup_.notify_all();
        return true;
    }


private:

    //---------------------------------------------------------------
    void work()
    {
         while(!abort_requested()) {
            wait_until_needed();

            if(!abort_requested()) {
                task_();
                status_.store(status::idle);
                isIdle_.notify_all();
            }
         }
         isIdle_.notify_all();
    }


    //---------------------------------------------------------------
    bool abort_requested() const noexcept {
        return (status_.load() == status::terminate);
    }


    //---------------------------------------------------------------
    void wait_until_needed() {
        std::unique_lock<std::mutex> lock{mutables_};

        //block execution of calling thread until woken up
        //(while loop protects against spurious wake-ups)
        while(status_.load() == status::idle) {
            wakeup_.wait(lock);
        }
    }


    //---------------------------------------------------------------
    mutable std::mutex mutables_;
    std::condition_variable wakeup_;
    std::condition_variable isIdle_;
    std::atomic<status> status_;
    task_type task_;
    std::thread thread_;
};


}  // namespace am


#endif
