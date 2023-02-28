#ifndef CARE_UTIL_HPP
#define CARE_UTIL_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <iterator>
#include <functional>
#include <cmath>
#include <numeric>
#include <vector>
#include <cassert>
#include <sstream>
#include <fstream>
#include <thread>
#include <atomic>
#include <chrono>
#include <type_traits>
#include <iostream>
#include <queue>


#include <cstdlib>


template<class T>
struct ProgressThread{

    template<class ProgressFunc>
    ProgressThread(T maxProgress_, ProgressFunc&& pfunc)
            : ProgressThread(maxProgress_, std::move(pfunc), [](auto seconds){return seconds;}){

    }

    template<class ProgressFunc, class SleepUpdateFunc>
    ProgressThread(T maxProgress_, ProgressFunc&& pfunc, SleepUpdateFunc&& sfunc)
            : starttime{std::chrono::system_clock::now()},
            sleepPeriod{std::chrono::seconds{1}},
            currentProgress{0},
            maxProgress{maxProgress_},
            showProgress{std::move(pfunc)},
            updateSleepPeriod{std::move(sfunc)},
            thread{[&](){threadFunc();}}{

        showProgress(0,0);
    }

    ~ProgressThread(){
        doCancel = true;
        thread.join();
    }

    ProgressThread(const ProgressThread&) = delete;
    ProgressThread(ProgressThread&&) = delete;
    ProgressThread operator=(const ProgressThread&) = delete;
    ProgressThread operator=(ProgressThread&&) = delete;

    void threadFunc(){
        std::this_thread::sleep_for(sleepPeriod);
        
        while(!doCancel && currentProgress < maxProgress){
            auto now = std::chrono::system_clock::now();
            const std::chrono::duration<double> duration = now - starttime;
            showProgress(currentProgress, duration.count());

            int sleepiters = sleepPeriod.count();
            while(sleepiters > 0 && !doCancel){
                std::this_thread::sleep_for(std::chrono::seconds{1});
                sleepiters--;
            }
            
            sleepPeriod = updateSleepPeriod(sleepPeriod);
        }
    }

    void cancel(){
        doCancel = true;
    }

    void finished(){
        doCancel = true;
        auto now = std::chrono::system_clock::now();
        const std::chrono::duration<double> duration = now - starttime;
        showProgress(currentProgress, duration.count());
    }

    void setProgress(T newProgress){
        assert(newProgress >= currentProgress);
        currentProgress = newProgress;
    }

    void addProgress(T add){
        // if(std::is_signed<T>::value){
        //     assert(add >= 0);
        // }
        
        currentProgress += add;
    }

    //std::atomic<bool> doCancel = false;
    bool doCancel = false;
    std::chrono::time_point<std::chrono::system_clock> starttime;
    std::chrono::seconds sleepPeriod{1};
    std::atomic<T> currentProgress;
    std::atomic<T> maxProgress;
    std::function<void(T, double)> showProgress;
    std::function<std::chrono::seconds(std::chrono::seconds)> updateSleepPeriod;
    std::thread thread;
};





template<class T>
class View{
private:
    const T* ptr;
    int nElements;
public:
    View() : View(nullptr, 0){}
    View(const T* p, int n) : ptr(p), nElements(n){}

    const T& operator[](int i) const{
        if(i >= nElements){
            throw std::runtime_error("Out-of-bounds view access!!!");
        }
        return ptr[i];
    }

    const T* begin() const{
        return ptr;
    }

    const T* end() const{
        return ptr + nElements;
    }

    int size() const{
        return int(std::distance(begin(), end()));
    }
};

inline
std::vector<std::string> split(const std::string& str, char c){
	std::vector<std::string> result;

	std::stringstream ss(str);
	std::string s;

	while (std::getline(ss, s, c)) {
		result.emplace_back(s);
	}

	return result;
}

// https://devblogs.microsoft.com/oldnewthing/20170102-00/?p=95095
// permutes data according to indices. indices is left in unspecified state
template<class T, class Index>
void permute(T* data, Index* indices, std::size_t N){
    using std::swap;

    for (size_t i = 0; i < N; i++) {
        auto current = i;
        while (i != indices[current]) {
            auto next = indices[current];
            swap(data[current], data[next]);
            indices[current] = current;
            current = next;
        }
        indices[current] = current;
    }
}

template<class T, std::size_t N>
struct FixedCapacityVector{
    std::size_t size_{};
    std::array<T,N> array{};

    void push_back(T elem){
        assert(size() < N);
        array[size_++] = std::move(elem);
    }

    std::size_t size() const noexcept{
        return size_;
    }

    auto begin() noexcept{
        return array.begin();
    }

    auto end() noexcept{
        return array.begin() + size();
    }

    T& operator[](std::size_t i) noexcept{
        assert(i < size());
        return array[i];
    }

    const T& operator[](std::size_t i) const noexcept{
        assert(i < size());
        return array[i];
    }
};


/*
    Performs a set union of multiple ranges into a single output range
*/

struct SetUnionHandle{
    std::vector<char> buffer{};
};

template<class OutputIt, class Iter>
OutputIt k_way_set_union(
        SetUnionHandle& handle,
        OutputIt outputbegin, 
        std::pair<Iter, Iter>* ranges,
        int numRanges){

    using OutputType = typename std::iterator_traits<OutputIt>::value_type;
    using InputType = typename std::iterator_traits<Iter>::value_type;

    static_assert(std::is_same<OutputType, InputType>::value, "");

    using T = InputType;

    //handle simple cases

    if(numRanges == 0){
        return outputbegin;
    }

    if(numRanges == 1){
        return std::copy(ranges[0].first, ranges[0].second, outputbegin);
    }

    if(numRanges == 2){
        return std::set_union(ranges[0].first,
                              ranges[0].second,
                              ranges[1].first,
                              ranges[1].second,
                              outputbegin);
    }

    //handle generic case

    //sort ranges by size
    std::sort(ranges, ranges + numRanges, [](const auto& l, const auto& r){
        return std::distance(l.first, l.second) < std::distance(r.first, r.second);
    });

    int totalElements = 0;
    for(int i = 0; i < numRanges; i++){
        const auto& range = ranges[i];
        totalElements += std::distance(range.first, range.second);
    }

    auto& temp = handle.buffer;
    temp.resize(sizeof(T) * totalElements);

    T* tempbegin = reinterpret_cast<T*>(temp.data());
    T* tempend = tempbegin;
    auto outputend = outputbegin;

    //to avoid a final copy from temp to outputrange, both ranges are swapped in the beginning if number of ranges is odd.
    if(numRanges % 2 == 1){
        std::swap(tempbegin, outputbegin);
        std::swap(tempend, outputend);
    }

    for(int k = 0; k < numRanges; k++){
        tempend = std::set_union(ranges[k].first,
                                  ranges[k].second,
                                  outputbegin,
                                  outputend,
                                  tempbegin);

        std::swap(tempbegin, outputbegin);
        std::swap(tempend, outputend);
    }

    return outputend;
}



template<class OutputIt, class Iter>
OutputIt k_way_set_union_with_priorityqueue(OutputIt outputbegin, std::vector<std::pair<Iter,Iter>>& ranges){
    using OutputType = typename std::iterator_traits<OutputIt>::value_type;
    using InputType = typename std::iterator_traits<Iter>::value_type;

    static_assert(std::is_same<OutputType, InputType>::value, "");

    //handle simple cases

    if(ranges.empty()){
        return outputbegin;
    }

    if(ranges.size() == 1){
        return std::copy(ranges[0].first, ranges[0].second, outputbegin);
    }

    if(ranges.size() == 2){
        return std::set_union(ranges[0].first,
                              ranges[0].second,
                              ranges[1].first,
                              ranges[1].second,
                              outputbegin);
    }

    //handle generic case

    struct PQval{
        int rangeIndex;
        Iter dataIter;

        bool operator<(const PQval& rhs) const{
            return *dataIter > *(rhs.dataIter); //order such that smallest element comes first in priority queue
        }
    };


    std::priority_queue<PQval> pq;

    for(int i = 0; i < int(ranges.size()); i++){
        const auto& range = ranges[i];
        if(std::distance(range.first, range.second) > 0){
            pq.emplace(PQval{i, range.first});
        }
    }

    //position of the previously added output element
    auto currentOutputIter = outputbegin;
    //points behind the last element in output range
    auto outputEnd = outputbegin;

    while(!pq.empty()){
        auto cur = pq.top();
        pq.pop();

        if(currentOutputIter != outputEnd){
            if(*currentOutputIter < *(cur.dataIter)){
                ++currentOutputIter;
                *currentOutputIter = *(cur.dataIter);
                ++outputEnd;
            }
        }else{
            *currentOutputIter = *(cur.dataIter); //the first output element can always be inserted
            ++outputEnd;
        }

         //if there is another element from the same range, add it to the priority queue
        ++cur.dataIter;
        if(cur.dataIter != ranges[cur.rangeIndex].second){
            pq.emplace(cur);
        }        
    }

    return outputEnd;
}



template<class OutputIt, class Iter>
OutputIt k_way_set_union_complicatedsort(OutputIt outputbegin, std::vector<std::pair<Iter,Iter>>& ranges){
    using OutputType = typename std::iterator_traits<OutputIt>::value_type;
    using InputType = typename std::iterator_traits<Iter>::value_type;

    static_assert(std::is_same<OutputType, InputType>::value, "");

    using T = InputType;

    //handle simple cases

    if(ranges.empty()){
        return outputbegin;
    }

    if(ranges.size() == 1){
        return std::copy(ranges[0].first, ranges[0].second, outputbegin);
    }

    if(ranges.size() == 2){
        return std::set_union(ranges[0].first,
                              ranges[0].second,
                              ranges[1].first,
                              ranges[1].second,
                              outputbegin);
    }

    //handle generic case
    auto sortAscending = [](const auto& l, const auto& r){
        return std::distance(l.first, l.second) > std::distance(r.first, r.second);
    };

    //sort ranges by size
    std::sort(ranges.begin(), ranges.end(), sortAscending);

    int totalElements = 0;
    for(const auto& range : ranges){
        totalElements += std::distance(range.first, range.second);
    }

    std::vector<std::vector<T>> tempunions;
    tempunions.resize(ranges.size() - 2);

    int iteration = 0;

    while(ranges.size() > 2){
        auto& range2 = ranges[ranges.size()-1];
        auto& range1 = ranges[ranges.size()-2];

        auto& result = tempunions[iteration];

        result.resize(std::distance(range1.first, range1.second) + std::distance(range2.first, range2.second));

        auto resEnd = std::set_union(range1.first,
                                  range1.second,
                                  range2.first,
                                  range2.second,
                                  result.data());

        auto newRange = std::make_pair(result.data(), resEnd);
        ranges.pop_back();
        ranges.pop_back();
        auto insertpos = std::lower_bound(ranges.begin(), ranges.end(), newRange, sortAscending);
        ranges.insert(insertpos, newRange);

        iteration++;
    }

    auto& range1 = ranges[0];
    auto& range2 = ranges[1];

    auto outputend = std::set_union(range1.first,
                              range1.second,
                              range2.first,
                              range2.second,
                              outputbegin);

    return outputend;
}






/*
    Merge ranges [first1, last1) and [first2, last2) into range beginning at d_first.
    If more than max_num_elements unique elements in the result range
    would occur >= threshold times, an empty range is returned.
*/
template<class OutputIt, class Iter1, class Iter2>
OutputIt merge_with_count_theshold(Iter1 first1, Iter1 last1,
                        Iter2 first2, Iter2 last2,
                        std::size_t threshold,
                        std::size_t max_num_elements,
                        OutputIt d_first){
    //static_assert(std::is_same<typename Iter1::value_type, typename Iter2::value_type>::value, "");
    static_assert(std::is_same<typename std::iterator_traits<Iter1>::value_type, typename std::iterator_traits<Iter2>::value_type>::value, "");

    using T = typename std::iterator_traits<Iter1>::value_type;

    OutputIt d_first_orig = d_first;

    T previous{};
    std::size_t count = 0;
    bool foundone = false;

    auto update = [&](){
        if(*d_first == previous){
            ++count;
        }else{
            if(count >= threshold){
                if(foundone){
                    --max_num_elements;
                }
                foundone = true;
            }
            previous = *d_first;
            count = 1;
        }
    };

    for (; first1 != last1 && max_num_elements > 0; ++d_first) {
        if (first2 == last2) {
            while(first1 != last1 && max_num_elements > 0){
                *d_first = *first1;
                update();
                ++d_first;
                ++first1;
            }
            break;
        }
        if (*first2 < *first1) {
            *d_first = *first2;
            ++first2;
        } else {
            *d_first = *first1;
            ++first1;
        }

        update();
    }

    while(first2 != last2 && max_num_elements > 0){
        *d_first = *first2;
        update();
        ++d_first;
        ++first2;
    }

    if(max_num_elements == 0 || (max_num_elements == 1 && count >= threshold))
        return d_first_orig;
    else
        return d_first;
}

template<class OutputIt, class Iter>
OutputIt k_way_set_intersection_naive(OutputIt destinationbegin, const std::vector<Iter>& iters){
    static_assert(std::is_same<typename OutputIt::value_type, typename Iter::value_type>::value, "");

    using T = typename Iter::value_type;

    // at least one range is invalid
    if(iters.size() % 2 == 1)
        return destinationbegin;

    if(iters.size() == 0)
        return destinationbegin;

    if(iters.size() == 4)
        return std::set_intersection(iters[0], iters[1], iters[2], iters[3], destinationbegin);

    std::size_t nranges = iters.size()/2;

    std::size_t num_elements = 0;
    for(std::size_t i = 0; i < iters.size() / 2; i++){
        num_elements += std::distance(iters[2*i + 0], iters[2*i+1]);
    }

    std::vector<T> tmpbuffer(num_elements);
    std::size_t merged_elements = 0;

    for(std::size_t i = 0; i < nranges; i++){
        auto src_begin = i % 2 == 0 ? tmpbuffer.begin() : destinationbegin;
        auto src_end = src_begin + merged_elements;
        auto dest_begin = i % 2 == 0 ? destinationbegin : tmpbuffer.begin();

        auto dest_end = std::set_intersection(src_begin, src_end,
                                    iters[2*i + 0], iters[2*i+1],
                                    dest_begin);

        merged_elements = std::distance(dest_begin, dest_end);
    }

    auto destinationend = destinationbegin + merged_elements;

    if(nranges % 2 == 0){
        destinationend = std::copy(tmpbuffer.begin(), tmpbuffer.begin() + merged_elements, destinationbegin);
    }

    return destinationend;
}


template<class OutputIt, class Iter>
OutputIt k_way_merge_naive(OutputIt destinationbegin, const std::vector<Iter>& iters){
    static_assert(std::is_same<typename OutputIt::value_type, typename Iter::value_type>::value, "");

    using T = typename Iter::value_type;

    // at least one range is invalid
    if(iters.size() % 2 == 1)
        return destinationbegin;

    if(iters.size() == 0)
        return destinationbegin;

    if(iters.size() == 2)
        return std::copy(iters[0], iters[1], destinationbegin);

    if(iters.size() == 4)
        return std::merge(iters[0], iters[1], iters[2], iters[3], destinationbegin);

    std::size_t nranges = iters.size()/2;

    std::size_t num_elements = 0;
    for(std::size_t i = 0; i < iters.size() / 2; i++){
        num_elements += std::distance(iters[2*i + 0], iters[2*i+1]);
    }

    std::vector<T> tmpbuffer(num_elements);
    std::size_t merged_elements = 0;

    for(std::size_t i = 0; i < nranges; i++){
        auto src_begin = i % 2 == 0 ? tmpbuffer.begin() : destinationbegin;
        auto src_end = src_begin + merged_elements;
        auto dest_begin = i % 2 == 0 ? destinationbegin : tmpbuffer.begin();

        auto dest_end = std::merge(src_begin, src_end,
                                    iters[2*i + 0], iters[2*i+1],
                                    dest_begin);

        merged_elements = std::distance(dest_begin, dest_end);
    }

    auto destinationend = destinationbegin + merged_elements;

    if(nranges % 2 == 0){
        destinationend = std::copy(tmpbuffer.begin(), tmpbuffer.begin() + merged_elements, destinationbegin);
    }

    return destinationend;
}

template<class OutputIt, class Iter>
OutputIt k_way_merge_naive_sortonce(OutputIt destinationbegin, const std::vector<Iter>& iters){
    static_assert(std::is_same<typename std::iterator_traits<OutputIt>::value_type, typename std::iterator_traits<Iter>::value_type>::value, "");

    using T = typename std::iterator_traits<Iter>::value_type;

    // at least one range is invalid
    if(iters.size() % 2 == 1)
        return destinationbegin;

    if(iters.size() == 0)
        return destinationbegin;

    if(iters.size() == 2)
        return std::copy(iters[0], iters[1], destinationbegin);

    if(iters.size() == 4)
        return std::merge(iters[0], iters[1], iters[2], iters[3], destinationbegin);

    std::size_t nranges = iters.size()/2;

    std::size_t num_elements = 0;
    for(std::size_t i = 0; i < iters.size() / 2; i++){
        num_elements += std::distance(iters[2*i + 0], iters[2*i+1]);
    }

    std::vector<int> indices(nranges);
    std::iota(indices.begin(), indices.end(), int(0));

    std::sort(indices.begin(), indices.end(), [&](auto l, auto r){
        auto ldist = std::distance(iters[2*l + 0], iters[2*l+1]);
        auto rdist = std::distance(iters[2*r + 0], iters[2*r+1]);
        return ldist < rdist;
    });

    std::vector<T> tmpbuffer(num_elements);
    std::size_t merged_elements = 0;

    for(std::size_t i = 0; i < nranges; i++){
        const int rangeid = indices[i];

        auto src_begin = i % 2 == 0 ? tmpbuffer.begin() : destinationbegin;
        auto src_end = src_begin + merged_elements;
        auto dest_begin = i % 2 == 0 ? destinationbegin : tmpbuffer.begin();

        auto dest_end = std::merge(src_begin, src_end,
                                    iters[2*rangeid + 0], iters[2*rangeid+1],
                                    dest_begin);

        merged_elements = std::distance(dest_begin, dest_end);
    }

    auto destinationend = destinationbegin + merged_elements;

    if(nranges % 2 == 0){
        destinationend = std::copy(tmpbuffer.begin(), tmpbuffer.begin() + merged_elements, destinationbegin);
    }

    return destinationend;
}

template<class OutputIt, class Iter>
OutputIt k_way_merge_sorted(OutputIt destinationbegin, std::vector<Iter> iters){
    static_assert(std::is_same<typename OutputIt::value_type, typename Iter::value_type>::value, "");

    using T = typename Iter::value_type;

    // at least one range is invalid
    if(iters.size() % 2 == 1)
        return destinationbegin;

    if(iters.size() == 0)
        return destinationbegin;

    if(iters.size() == 2)
        return std::copy(iters[0], iters[1], destinationbegin);

    if(iters.size() == 4)
        return std::merge(iters[0], iters[1], iters[2], iters[3], destinationbegin);

    std::size_t nranges = iters.size()/2;

    std::size_t num_elements = 0;
    for(std::size_t i = 0; i < iters.size() / 2; i++){
        num_elements += std::distance(iters[2*i + 0], iters[2*i+1]);
    }

    std::vector<int> indices(nranges);

    std::vector<std::vector<T>> buffers(nranges-1);
    //for(auto& buffer : buffers)
	//buffer.reserve(num_elements);

    auto destinationend = destinationbegin;

    int pending_merges = nranges-1;

    while(pending_merges > 0){
	//for(int i = 0; i < pending_merges+1; i++)
	//	std::cout << "range " << i << ", " << std::distance(iters[2*i + 0], iters[2*i+1]) << " elements" << std::endl;

    	if(pending_merges > 1){
    		indices.resize(pending_merges+1);
    	    	std::iota(indices.begin(), indices.end(), int(0));

    		std::sort(indices.begin(), indices.end(), [&](auto l, auto r){
    			auto ldist = std::distance(iters[2*l + 0], iters[2*l+1]);
    			auto rdist = std::distance(iters[2*r + 0], iters[2*r+1]);
    			return ldist < rdist;
    		});

    		int lindex = indices[0];
    		int rindex = indices[1];

    		std::size_t ldist = std::distance(iters[2*lindex + 0], iters[2*lindex+1]);
    		std::size_t rdist = std::distance(iters[2*rindex + 0], iters[2*rindex+1]);

    		int bufferindex = nranges-1 - pending_merges;
    		auto& buffer = buffers[bufferindex];
    		buffer.resize(ldist+rdist);

    		std::merge(iters[2*lindex + 0], iters[2*lindex+1],
    		           iters[2*rindex + 0], iters[2*rindex+1],
    		           buffer.begin());

    		iters[2*lindex+0] = buffer.begin();
    		iters[2*lindex+1] = buffer.end();
    		iters.erase(iters.begin() + (2*rindex+0), iters.begin() + (2*rindex+1) + 1);

    	}else{
    		int lindex = 0;
    		int rindex = 1;
    		destinationend = std::merge(iters[2*lindex + 0], iters[2*lindex+1],
    					   iters[2*rindex + 0], iters[2*rindex+1],
    					   destinationbegin);
    	}

    	--pending_merges;
    }

    return destinationend;
}


template<class OutputIter, class Iter, class FlagIter>
OutputIter select_if(Iter begin, Iter end, FlagIter flags, OutputIter output){

    while(begin != end){
        if(*flags){
            *output = *begin;
            ++output;
        }

        ++begin;
        ++flags;
    }

    return output;
}



/*
    Removes elements from sorted range which occure less than k times.
    Returns end of the new range
*/

template<class Iter, class Equal>
Iter remove_by_count(Iter first,
                    Iter last,
                    std::size_t k,
                    Equal equal){

    using T = typename Iter::value_type;

    if(std::distance(first, last) == 0) return first;

    Iter copytobegin = first;
    Iter copyfrombegin = first;
    Iter copyfromend = first;
    ++copyfromend;

    const T* curElem = &(*copyfrombegin);
    std::size_t count = 1;

    while(copyfromend != last){

        if(equal(*curElem, *copyfromend)){
            ++count;
        }else{
            if(count < k){
                copyfrombegin = copyfromend;
            }else{
                copytobegin = std::copy(copyfrombegin, copyfromend, copytobegin);
                copyfrombegin = copyfromend;
            }

            curElem = &(*copyfromend);
            count = 1;
        }

        ++copyfromend;
    }

    //handle last element
    if(count >= k)
        copytobegin = std::copy(copyfrombegin, copyfromend, copytobegin);

    return copytobegin;
}

template<class Iter>
Iter remove_by_count(Iter first,
                        Iter last,
                        std::size_t k){
    using T = typename Iter::value_type;
    return remove_by_count(first, last, k, std::equal_to<T>{});
}

/*
    Removes elements from sorted range which occure less than k times.
    If a range of equals elements is greater than or equal to k, only the first element of this range is kept.
    Returns end of the new range.
    The new range is empty if there are more than max_num_elements unique elements
*/

template<class Iter, class Equal>
Iter remove_by_count_unique_with_limit(Iter first,
                    Iter last,
                    std::size_t k,
                    std::size_t max_num_elements,
                    Equal equal){
    using T = typename Iter::value_type;

    constexpr std::size_t elements_to_copy = 1;

    if(std::distance(first, last) == 0) return first;
    if(elements_to_copy > k) return first;

    Iter copytobegin = first;
    Iter copyfrombegin = first;
    Iter copyfromend = first;
    ++copyfromend;

    std::size_t num_copied_elements = 0;

    const T* curElem = &(*copyfrombegin);
    std::size_t count = 1;

    while(copyfromend != last && num_copied_elements <= max_num_elements){

        if(equal(*curElem, *copyfromend)){
            ++count;
        }else{
            if(count < k){
                copyfrombegin = copyfromend;
            }else{
                copytobegin = std::copy_n(copyfrombegin, elements_to_copy, copytobegin);
                copyfrombegin = copyfromend;
                num_copied_elements += elements_to_copy;
            }

            curElem = &(*copyfromend);
            count = 1;
        }

        ++copyfromend;
    }

    //handle last element
    if(count >= k)
        copytobegin = std::copy_n(copyfrombegin, elements_to_copy, copytobegin);

    if(num_copied_elements > max_num_elements)
        return first;

    return copytobegin;
}

template<class Iter>
Iter remove_by_count_unique_with_limit(Iter first,
                    Iter last,
                    std::size_t k,
                    std::size_t max_num_elements){

    using T = typename Iter::value_type;
    return remove_by_count_unique_with_limit(first, last, k, max_num_elements, std::equal_to<T>{});
}


/*
    Essentially performs std::set_union(first1, last1, first2, last2, d_first)
    but limits the allowed result size to n.
    If the result would contain more than n elements, d_first is returned, i.e. the result range is empty
*/
template<class InputIt1, class InputIt2, class OutputIt>
OutputIt set_union_n_or_empty(InputIt1 first1, InputIt1 last1,
                   InputIt2 first2, InputIt2 last2,
                   std::size_t n,
                   OutputIt d_first){

    const OutputIt d_first_old = d_first;

    for (; first1 != last1 && n > 0; ++d_first, --n) {
        if (first2 == last2){
            const std::size_t remaining = std::distance(first1, last1);
            if(remaining > n)
                return d_first_old;
            else
                return std::copy_n(first1, std::min(remaining, n), d_first);
        }

        if (*first2 < *first1) {
            *d_first = *first2++;
        } else {
            *d_first = *first1;
            if (!(*first1 < *first2))
                ++first2;
            ++first1;
        }
    }

    const std::size_t remaining = std::distance(first2, last2);
    if(remaining > n)
        return d_first_old;
    else
        return std::copy_n(first2, std::min(remaining, n), d_first);
}

/*
    Essentially performs std::set_union(first1, last1, first2, last2, d_first)
    but limits the allowed result size to n.
*/
template<class InputIt1, class InputIt2, class OutputIt>
OutputIt set_union_n(InputIt1 first1, InputIt1 last1,
                   InputIt2 first2, InputIt2 last2,
                   std::size_t n,
                   OutputIt d_first){

    for (; first1 != last1 && n > 0; ++d_first, --n) {
        if (first2 == last2){
            const std::size_t remaining = std::distance(first1, last1);
            return std::copy_n(first1, std::min(remaining, n), d_first);
        }

        if (*first2 < *first1) {
            *d_first = *first2++;
        } else {
            *d_first = *first1;
            if (!(*first1 < *first2))
                ++first2;
            ++first1;
        }
    }

    const std::size_t remaining = std::distance(first2, last2);
    return std::copy_n(first2, std::min(remaining, n), d_first);
}


/*
    Essentially performs std::set_intersection(first1, last1, first2, last2, d_first)
    but limits the allowed result size to n.
    If the result would contain more than n elements, d_first is returned, i.e. the result range is empty
*/
template<class InputIt1, class InputIt2, class OutputIt>
OutputIt set_intersection_n_or_empty(InputIt1 first1, InputIt1 last1,
                   InputIt2 first2, InputIt2 last2,
                   std::size_t n,
                   OutputIt d_first){

        const OutputIt d_first_old = d_first;
        ++n;
        while (first1 != last1 && first2 != last2 && n > 0) {
            if (*first1 < *first2) {
                ++first1;
            } else {
                if (!(*first2 < *first1)) {
                    *d_first++ = *first1++;
                    --n;
                }
                ++first2;
            }
        }
        if(n == 0){
           //intersection contains at least n+1 elements, return empty range
           return d_first_old;
        }

        return d_first;
}


template<class InputIt1, class InputIt2,
         class OutputIt1, class OutputIt2,
         class CompareOther, class CompareSame1, class CompareSame2,
         class LessThan12, class LessThan21>
std::pair<OutputIt1, OutputIt2> flagPairedCandidates(
    InputIt1 first1, 
    InputIt1 last1,
    InputIt2 first2, 
    InputIt2 last2,
    OutputIt1 d_first1, 
    OutputIt2 d_first2, 
    CompareOther isSameReadPairInOtherRange,
    CompareSame1 isSameReadPairInSameRange1,
    CompareSame2 isSameReadPairInSameRange2,
    LessThan12 lessThan12,
    LessThan21 lessThan21
){
    auto outputbegin1 = d_first1;
    auto outputbegin2 = d_first2;

    while (first1 != last1 && first2 != last2) {
        const auto nextfirst1 = std::next(first1);
        const auto nextfirst2 = std::next(first2);

        int elems1 = 1;
        if(nextfirst1 != last1){
            if(isSameReadPairInSameRange1(*first1, *nextfirst1)){
                elems1 = 2;
            }
        }

        int elems2 = 1;
        if(nextfirst2 != last2){
            if(isSameReadPairInSameRange2(*first2, *nextfirst2)){
                elems2 = 2;
            }
        }

        if(elems1 == 1 && elems2 == 1){
            if(isSameReadPairInOtherRange(*first1, *first2)){
                //if *first1 != *first2
                if(lessThan12(*first1,*first2) || lessThan21(*first2, *first1)){
                    *d_first1++ = *first1++;
                    *d_first2++ = *first2++;
                }else{
                    ++first1;
                    ++first2;
                }
            }else{
                if(lessThan12(*first1,*first2)){
                    ++first1;
                }else{
                    ++first2;
                }
            }
        }else if (elems1 == 2 && elems2 == 2){
            if(isSameReadPairInOtherRange(*first1, *first2)){
                *d_first1++ = *first1++;
                *d_first2++ = *first2++;
                *d_first1++ = *first1++;
                *d_first2++ = *first2++;
            }else{
                if(lessThan12(*first1,*first2)){
                    ++first1;
                    ++first1;
                }else{
                    ++first2;
                    ++first2;
                }
            }

        }else if (elems1 == 2 && elems2 == 1){
            if(isSameReadPairInOtherRange(*first1, *first2)){
                //if *first1 == *first2 , e.g (4,5) , (4)
                if(!lessThan12(*first1,*first2) && !lessThan21(*first2, *first1)){
                    //discard first entry of first range, keep rest
                    ++first1;
                    *d_first1++ = *first1++;
                    *d_first2++ = *first2++;
                }else{ // e.g (4,5) , (5)
                    //keep first entry of first range, discard second entry
                    *d_first1++ = *first1++;
                    *d_first2++ = *first2++;
                    ++first1;
                }
            }else{
                if(lessThan12(*first1,*first2)){
                    ++first1;
                    ++first1;
                }else{
                    ++first2;
                }
            }
            
        }else {
            //(elems1 == 1 && elems2 == 2)

            if(isSameReadPairInOtherRange(*first1, *first2)){
                //if *first1 == *first2 , e.g (4) , (4,5)
                if(!lessThan12(*first1,*first2) && !lessThan21(*first2, *first1)){
                    //discard first entry of second range, keep rest
                    ++first2;
                    *d_first1++ = *first1++;
                    *d_first2++ = *first2++;
                }else{
                    //keep first entry of second range, discard second entry of second range
                    *d_first1++ = *first1++;
                    *d_first2++ = *first2++;
                    ++first2;
                }
            }else{
                if(lessThan12(*first1,*first2)){
                    ++first1;
                }else{
                    ++first2;
                    ++first2;
                }
            }            
        }
    }
    
    if(std::distance(outputbegin1, d_first1) != std::distance(outputbegin2, d_first2))
        throw std::runtime_error("Error flagging paired candidates");

    return std::make_pair(d_first1, d_first2);
}


/*
    Input: Two sorted ranges of read ids.
    Two read ids x,y form a pair if x / 2 == y / 2

    Output ranges will contain positions of those read ids for which its pair id exists in the respective other range
*/
template<class InputIt1, class InputIt2,
         class OutputIt1, class OutputIt2>
std::pair<OutputIt1, OutputIt2> findPositionsOfPairedReadIds(
    InputIt1 first1, 
    InputIt1 last1,
    InputIt2 first2, 
    InputIt2 last2,
    OutputIt1 d_first1, 
    OutputIt2 d_first2
){
    std::vector<int> indexlist1(std::distance(first1, last1));
    std::iota(indexlist1.begin(), indexlist1.end(), 0);

    std::vector<int> indexlist2(std::distance(first2, last2));
    std::iota(indexlist2.begin(), indexlist2.end(), 0);

    auto isSameReadPairInOtherRange = [&](int l, int r){
        return *(first1 + l) / 2 == *(first2 + r) / 2;
    };

    auto isSameReadPairInSameRange1 = [&](int l, int r){
        return*(first1 + l)  / 2 == *(first1 + r) / 2;
    };

    auto isSameReadPairInSameRange2 = [&](int l, int r){
        return *(first2 + l) / 2 == *(first2 + r) / 2;
    };

    auto lessThan12 = [&](int l, int r){
        return *(first1 + l) < *(first2 + r);
    };

    auto lessThan21 = [&](int l, int r){
        return *(first2 + l) < *(first1 + r);
    };

    return flagPairedCandidates(
        indexlist1.begin(), 
        indexlist1.end(),
        indexlist2.begin(), 
        indexlist2.end(),
        d_first1, 
        d_first2, 
        isSameReadPairInOtherRange,
        isSameReadPairInSameRange1,
        isSameReadPairInSameRange2,
        lessThan12,
        lessThan21
    );
}



/*
    Input: Two sorted ranges of read ids.
    Two read ids x,y form a pair if x / 2 == y / 2

    Output ranges will contain read ids for which its pair id exists in the respective other range
*/
template<class InputIt1, class InputIt2,
         class OutputIt1, class OutputIt2>
std::pair<OutputIt1, OutputIt2> findPairedReadIds(
    InputIt1 first1, 
    InputIt1 last1,
    InputIt2 first2, 
    InputIt2 last2,
    OutputIt1 d_first1, 
    OutputIt2 d_first2
){

    auto isSameReadPairInOtherRange = [](const auto& l, const auto& r){
        return l / 2 == r / 2;
    };

    auto isSameReadPairInSameRange1 = [](const auto& l, const auto& r){
        return l / 2 == r / 2;
    };

    auto isSameReadPairInSameRange2 = [](const auto& l, const auto& r){
        return l / 2 == r / 2;
    };

    auto lessThan12 = [](const auto& l, const auto& r){
        return l < r;
    };

    auto lessThan21 = [](const auto& l, const auto& r){
        return l < r;
    };

    return flagPairedCandidates(
        first1, 
        last1,
        first2, 
        last2,
        d_first1, 
        d_first2, 
        isSameReadPairInOtherRange,
        isSameReadPairInSameRange1,
        isSameReadPairInSameRange2,
        lessThan12,
        lessThan21
    );
}



#endif
