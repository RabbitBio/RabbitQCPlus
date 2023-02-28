/*****************************************************************************
 *
 * AM utilities
 *
 * released under MIT license
 *
 * 2008-2018 André Müller
 *
 *****************************************************************************/

#ifndef AMLIB_PARALLEL_TASK_QUEUE_H_
#define AMLIB_PARALLEL_TASK_QUEUE_H_

#include <deque>
#include <vector>
#include <algorithm>
#include <iostream>
#include <functional>
#include <memory>
#include <atomic>

#include "task_thread.h"


namespace am {


/*****************************************************************************
 *
 * @brief    simple queue that holds tasks and executes them in parallel
 *
 * @tparam   Task : assumed to be cheap to copy or move
 *                  (e.g. std::function, std::reference_wrapper, ...)
 *                  Task::operator()() must be defined
 *
 *****************************************************************************/
template<class Task>
class parallel_task_queue
{
public:
    //---------------------------------------------------------------
    using task_type = Task;


private:
    //---------------------------------------------------------------
    class task_executor {
    public:
        explicit
        task_executor(parallel_task_queue* callback = nullptr,
                      task_type&& task = task_type{}) noexcept
        :
            queue_{callback},
            task_{std::move(task)}
        {}

        void operator () () {
            task_();

            std::unique_lock<std::mutex> lock{queue_->busyMtx_};
            --queue_->running_;
            //queue_->notify_task_complete();
            (*queue_).notify_task_complete();
        }

    private:
        parallel_task_queue* queue_;
        task_type task_;
    };

    friend class task_executor;


public:
    //---------------------------------------------------------------
    explicit
    parallel_task_queue(unsigned int concurrency =
        std::thread::hardware_concurrency())
    :
        waitMtx_{}, enqueueMtx_{},
        active_{true}, hasWaiting_{false},
        running_{0},
        waiting_{},
        workers_(std::min(concurrency, std::thread::hardware_concurrency())),
        isDone_{},
        busyMtx_{}, 
        isBusy_{},
        isWaitingForTask_{},
        scheduler_{ [&] { schedule(); }}        
    {}

    //-----------------------------------------------------
    parallel_task_queue(const parallel_task_queue&) = delete;
    parallel_task_queue(parallel_task_queue&&) = delete;


    //-----------------------------------------------------
    ~parallel_task_queue()
    {
        clear();
        //make sure that scheduler wakes up and terminates
        std::unique_lock<std::recursive_mutex> lock(enqueueMtx_);
        active_.store(false);
        isWaitingForTask_.notify_one();
        lock.unlock();
        //wait for scheduler to terminate
        scheduler_.join();
    }


    //---------------------------------------------------------------
    parallel_task_queue& operator = (const parallel_task_queue&) = delete;
    parallel_task_queue& operator = (parallel_task_queue&&) = delete;


    //---------------------------------------------------------------
    void
    enqueue(const task_type& t)
    {
        std::lock_guard<std::recursive_mutex> lock(enqueueMtx_);
        waiting_.push_back(t);
        hasWaiting_.store(true);
        isWaitingForTask_.notify_one();
    }
    //-----------------------------------------------------
    void
    enqueue(task_type&& t)
    {
        std::lock_guard<std::recursive_mutex> lock(enqueueMtx_);
        waiting_.push_back(std::move(t));
        hasWaiting_.store(true);
        isWaitingForTask_.notify_one();
    }
    //-----------------------------------------------------
    template<class InputIterator>
    void
    enqueue(InputIterator first, InputIterator last)
    {
        std::lock_guard<std::recursive_mutex> lock(enqueueMtx_);
        waiting_.insert(begin(waiting_), first, last);
        hasWaiting_.store(true);
        isWaitingForTask_.notify_one();
    }
    //-----------------------------------------------------
    void
    enqueue(std::initializer_list<task_type> il)
    {
        std::lock_guard<std::recursive_mutex> lock(enqueueMtx_);
        waiting_.insert(waiting_.end(), il);
        hasWaiting_.store(true);
        isWaitingForTask_.notify_one();
    }
    //-----------------------------------------------------
    template<class... Args>
    void
    enqueue_emplace(Args&&... args)
    {
        std::lock_guard<std::recursive_mutex> lock(enqueueMtx_);
        waiting_.emplace_back(std::forward<Args>(args)...);
        hasWaiting_.store(true);
        isWaitingForTask_.notify_one();
    }


    /****************************************************************
     * @brief removes first waiting task that compares equal to 't'
     */
    bool
    try_remove(const task_type& t)
    {
        std::lock_guard<std::recursive_mutex> lock(enqueueMtx_);
        auto it = std::find(waiting_.begin(), waiting_.end(), t);
        if(it != waiting_.end()) {
            waiting_.erase(it);
            if(waiting_.empty()) hasWaiting_.store(false);
            return true;
        }
        return false;
    }


    //---------------------------------------------------------------
    void
    clear() {
        std::lock_guard<std::recursive_mutex> lock(enqueueMtx_);
        waiting_.clear();
        hasWaiting_.store(false);
    }


    //---------------------------------------------------------------
    unsigned int
    concurrency() const noexcept {
        return workers_.size();
    }

    //---------------------------------------------------------------
    bool
    empty() const noexcept {
        return !hasWaiting_.load();
    }


    //-----------------------------------------------------
    /// @return number of waiting tasks
    std::size_t
    waiting() const noexcept {
        std::lock_guard<std::recursive_mutex> lock{enqueueMtx_};
        return waiting_.size();
    }
    std::size_t
    unsafe_waiting() const noexcept {
        return waiting_.size();
    }

    //-----------------------------------------------------
    /// @return number of running tasks
    std::size_t
    running() const noexcept {
        return running_.load();
    }

    //-----------------------------------------------------
    /// @return true, if all threads are working on a task
    bool
    busy() const noexcept {
        return running_.load() >= int(workers_.size());
    }

    //-----------------------------------------------------
    bool
    complete() const noexcept {
        std::lock_guard<std::recursive_mutex> lock{enqueueMtx_};
        return empty() && (running() < 1);
    }


    //---------------------------------------------------------------
    /**
     * @brief block execution of calling thread until all tasks are completed which are currently pending or running 
     */
    void wait()
    {
        std::unique_lock<std::recursive_mutex> enqueuelock{enqueueMtx_};

        auto barrierCount = std::make_shared<std::atomic<int>>(int(concurrency()));
        
        std::condition_variable cv;

        auto barrierFunc = [&, barrierCount](){
            std::unique_lock<std::mutex> lock{waitMtx_};
            --(*barrierCount);

            if((*barrierCount) == 0){
                cv.notify_all();
            }else{
                cv.wait(lock, [&](){return (*barrierCount) == 0;});
            }
        };

        for(int i = 0; i < int(concurrency()); i++){
            enqueue(barrierFunc);
        }

        enqueuelock.unlock(); 

        std::unique_lock<std::mutex> lock{waitMtx_};
        if((*barrierCount) != 0){
            cv.wait(lock, [&](){return (*barrierCount) == 0;});
        }
    }


private:
    //-----------------------------------------------------
    void notify_task_complete() {
        isBusy_.notify_one();
    }


    //---------------------------------------------------------------
    void try_assign_tasks()
    {
        for(auto& worker : workers_) {
            if(worker.available()) {
                std::lock_guard<std::recursive_mutex> lock{enqueueMtx_};
                if(waiting_.empty()) {
                    hasWaiting_.store(false);
                    return;
                }
                if(worker(this, std::move(waiting_.front()))) {
                    ++running_;
                    waiting_.pop_front();
                    if(waiting_.empty()) {
                        hasWaiting_.store(false);
                        return;
                    }
                    else if(busy()) {
                        return;
                    }
                }
            }
        }
    }

    //---------------------------------------------------------------
    /// @brief this will run in a separate, dedicated thread
    void schedule()
    {
        while(active_.load()) {
            if(busy()) {
                std::unique_lock<std::mutex> lock{busyMtx_};
                if(busy()){
                    isBusy_.wait(lock, [this]{ return !busy(); });
                }
            }
            else if(!empty()) {
                try_assign_tasks();
            }
            else{
                std::unique_lock<std::recursive_mutex> lock{enqueueMtx_};
                isWaitingForTask_.wait(lock, [this](){return !active_.load() || !empty();});
            }
        }
    }


    //---------------------------------------------------------------
    mutable std::mutex waitMtx_;
    mutable std::recursive_mutex enqueueMtx_;
    std::atomic_bool active_;
    std::atomic_bool hasWaiting_;
    std::atomic_int running_;
    std::deque<task_type> waiting_;
    std::vector<task_thread<task_executor>> workers_;
    std::condition_variable isDone_;
    std::mutex busyMtx_;
    std::condition_variable isBusy_;
    std::condition_variable_any isWaitingForTask_;
    std::thread scheduler_;
};



//-------------------------------------------------------------------
/// @brief convenience alias
using parallel_queue = parallel_task_queue<std::function<void()>>;


} // namespace am


#endif
