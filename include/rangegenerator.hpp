#ifndef CARE_RANGE_GENERATOR_HPP
#define CARE_RANGE_GENERATOR_HPP

#include <config.hpp>

#include <atomic>
#include <mutex>
#include <stdexcept>
#include <vector>
#include <numeric>

namespace care{

    template<class Count_t>
    struct RangeGenerator{
    private:
        Count_t begin;
        Count_t end;
        Count_t current;
        bool isEmpty;
        std::mutex mutex;

    public:
        RangeGenerator(Count_t end) : RangeGenerator(Count_t{0}, end, Count_t{0}) {}

        RangeGenerator(Count_t begin, Count_t end, Count_t initial_value)
            : begin(begin), end(end), current(initial_value), isEmpty(initial_value >= end) {}

        bool empty(){
            std::lock_guard<std::mutex> lm(mutex);
            return isEmpty;
        }

        void skip(Count_t n){
            std::lock_guard<std::mutex> lm(mutex);
            Count_t remaining = end - current;
            Count_t resultsize = std::min(remaining, n);
            current += resultsize;

            if(current == end){
                isEmpty = true;
            }
        }

        std::vector<Count_t> next_n(Count_t n){
            std::lock_guard<std::mutex> lm(mutex);
            if(isEmpty)
                return {};

            Count_t remaining = end - current;

            Count_t resultsize = std::min(remaining, n);

            std::vector<Count_t> result(resultsize);
            std::iota(result.begin(), result.end(), current);

            current += resultsize;

            if(current == end)
                isEmpty = true;

            return result;
        }

        //buffer must point to memory location of at least n elements
        //returns past the end iterator
        template<class Iter>
        Iter next_n_into_buffer(Count_t n, Iter buffer){
            std::lock_guard<std::mutex> lm(mutex);
            if(isEmpty)
                return buffer;

            Count_t remaining = end - current;

            Count_t resultsize = std::min(remaining, n);

            std::iota(buffer, buffer + resultsize, current);

            current += resultsize;

            if(current == end)
                isEmpty = true;

            return buffer + resultsize;
        }

        Count_t getBegin() const{
            return begin;
        }

        Count_t getEnd() const{
            return end;
        }

        Count_t getCurrentUnsafe() const{
            return current;
        }
    };

    template<class Iterator>
    struct IteratorRangeTraversal{
    private:

        bool isEmpty;
        Iterator begin;
        Iterator end;
        Iterator current;
        std::mutex mutex;

    public:
        IteratorRangeTraversal(Iterator begin_, Iterator end_) : isEmpty(std::distance(begin_, end_) <= 0), begin(begin_), end(end_), current(begin_){}

        void reset(Iterator begin_, Iterator end_){
            std::lock_guard<std::mutex> lm(mutex);
            isEmpty = std::distance(begin_, end_) <= 0;
            begin = begin_;
            end = end_;
            current = begin_;
        }

        bool empty(){
            std::lock_guard<std::mutex> lm(mutex);
            return isEmpty;
        }

        void skip(std::size_t n){
            std::lock_guard<std::mutex> lm(mutex);
            const std::size_t remaining = std::distance(current, end);
            const std::size_t resultsize = std::min(remaining, n);
            std::advance(current, resultsize);

            if(current == end){
                isEmpty = true;
            }
        }

        template<class Func>
        void process_next_n(std::size_t n, Func callback){
            std::unique_lock<std::mutex> lock(mutex);
            if(!isEmpty){
                const std::size_t remaining = std::distance(current, end);
                const std::size_t resultsize = std::min(remaining, n);
                auto callbackbegin = current;
                std::advance(current, resultsize);
                auto callbackend = current;

                if(current == end){
                    isEmpty = true;
                }
                lock.unlock();
                callback(callbackbegin, callbackend);
            }else{
                lock.unlock();
                callback(end, end);
            }
        }
    };

    template<class Iterator>
    IteratorRangeTraversal<Iterator> makeIteratorRangeTraversal(Iterator begin, Iterator end){
        return IteratorRangeTraversal<Iterator>(begin, end);
    }

}



#endif
