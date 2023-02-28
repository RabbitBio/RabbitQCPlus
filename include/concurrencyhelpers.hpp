#ifndef CARE_CONCURRENCY_HELPERS_HPP
#define CARE_CONCURRENCY_HELPERS_HPP

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <cassert>

#include <moodycamel/readerwriterqueue/readerwriterqueue.h>
#include <moodycamel/concurrentqueue/blockingconcurrentqueue.h>

namespace care{

struct SyncFlag{
    std::atomic<bool> busy{false};
    std::mutex m;
    std::condition_variable cv;

    void setBusy(){
        assert(busy == false);
        busy = true;
    }

    bool isBusy() const{
        return busy;
    }

    void wait(){
        if(isBusy()){
            std::unique_lock<std::mutex> l(m);
            while(isBusy()){
                cv.wait(l);
            }
        }
    }

    void signal(){
        std::unique_lock<std::mutex> l(m);
        busy = false;
        cv.notify_all();
    }        
};

template<class T>
struct WaitableData{
    T data;
    SyncFlag syncFlag;

    void setBusy(){
        syncFlag.setBusy();
    }

    bool isBusy() const{
        return syncFlag.isBusy();
    }

    void wait(){
        syncFlag.wait();
    }

    void signal(){
        syncFlag.signal();
    } 
};


template<class T>
struct SimpleSingleProducerSingleConsumerQueue{
    bool useDefaultElement = false;

    std::queue<T> queue;
    std::mutex mutex;
    std::condition_variable cv;

    T defaultElement{};

    void push(T item){
        std::lock_guard<std::mutex> lg(mutex);
        queue.emplace(std::move(item));
        cv.notify_one();
    }

    //wait until queue is not empty, then remove first element from queue and return it
    T pop(){
        std::unique_lock<std::mutex> ul(mutex);

        while(queue.empty() && !useDefaultElement){
            cv.wait(ul);
        }

        if(!queue.empty()){
            T item = queue.front();
            queue.pop();
            return item;
        }else{
            assert(useDefaultElement);
            return defaultElement;
        }
    }

    void enableDefaultElement(T def){
        std::lock_guard<std::mutex> lg(mutex);
        if(!useDefaultElement){
            defaultElement = def;
            useDefaultElement = true;
            cv.notify_all();
        }
    }

    //Wait until func() == false or queue is not empty.
    //if func() == false, return defaultElement, else the first element in queue
    template<class Func>
    T popOrDefault(Func func, T defaultElement){
        std::unique_lock<std::mutex> ul(mutex);

        if(func() || !queue.empty()){
            while(func() && queue.empty()){
                cv.wait(ul);
            }
            if(!queue.empty()){
                T item = queue.front();
                queue.pop();
                return item;
            }else{
                return defaultElement;
            }
        }else{
            return defaultElement;
        }
    }
};


template<class T>
struct SimpleMultiProducerMultiConsumerQueue{

    int numActiveProducers = 0;
    std::queue<T> queue;
    std::mutex mutex;
    std::condition_variable cv;

    SimpleMultiProducerMultiConsumerQueue() : SimpleMultiProducerMultiConsumerQueue(1){}

    SimpleMultiProducerMultiConsumerQueue(int prod) : numActiveProducers(prod){
        
    }

    void push(T item){
        std::lock_guard<std::mutex> lg(mutex);
        queue.emplace(std::move(item));
        cv.notify_one();
    }

    //wait until queue is not empty, then remove first element from queue and return it
    T pop(){
        std::unique_lock<std::mutex> ul(mutex);

        while(queue.empty()){
            cv.wait(ul);
        }

        T item = queue.front();
        queue.pop();
        return item;
    }

    void decreaseActiveProducers(){
        std::lock_guard<std::mutex> lg(mutex);
        if(numActiveProducers > 0){
            numActiveProducers--;
            if(numActiveProducers == 0){
                cv.notify_all();
            }
        }
    }

    T pop(T defaultElement){
        std::unique_lock<std::mutex> ul(mutex);

        while(queue.empty() && numActiveProducers > 0){
            cv.wait(ul);
        }

        if(!queue.empty()){
            T item = queue.front();
            queue.pop();
            return item;
        }else{
            return defaultElement;
        }
    }
};








template<class T>
struct SingleProducerSingleConsumerQueue{
    moodycamel::BlockingReaderWriterQueue<T> queue;

    void push(T item){
        queue.enqueue(item); 
    }

    //wait until queue is not empty, then remove first element from queue and return it
    T pop(){
        T item{};
        queue.wait_dequeue(item);
        return item;
    }

    //Wait until func() == false or queue is not empty.
    //if func() == false, return defaultElement, else the first element in queue
    template<class Func>
    T popOrDefault(Func func, T defaultElement){
        T item;
        bool success = false;

        while(!success && (func() || queue.size_approx() > 0)){
            success = queue.try_dequeue(item);
            
            if(!success){
                std::this_thread::sleep_for(std::chrono::milliseconds{1});
            }
        }

        if(!success){
            return defaultElement;
        }else{
            return item;
        }
    }
};


template<class T>
struct MultiProducerMultiConsumerQueue{
    moodycamel::BlockingConcurrentQueue<T> queue;

    void push(T item){
        queue.enqueue(item); 
    }

    //wait until queue is not empty, then remove first element from queue and return it
    T pop(){
        T item;
        queue.wait_dequeue(item);
        return item;
    }
};

template<class T>
using SimpleConcurrentQueue = SimpleSingleProducerSingleConsumerQueue<T>;
//using SimpleConcurrentQueue = SingleProducerSingleConsumerQueue<T>;


}

#endif