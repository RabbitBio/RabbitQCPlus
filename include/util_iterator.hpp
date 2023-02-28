#ifndef CARE_UTIL_ITERATOR_HPP
#define CARE_UTIL_ITERATOR_HPP

#include <hpc_helpers.cuh>

#include <iterator>


// iterator multiplier

template<class Iter>
struct IteratorMultiplier{
    using value_type = typename std::iterator_traits<Iter>::value_type;

    int factor{};
    Iter data{};

    HOSTDEVICEQUALIFIER
    IteratorMultiplier(Iter data_, int factor_)
        : factor(factor_), data(data_){

    }

    HOSTDEVICEQUALIFIER
    value_type operator()(int i) const{
        return *(data + (i / factor));
    }
};

template<class Iter>
HOSTDEVICEQUALIFIER
IteratorMultiplier<Iter> make_iterator_multiplier(Iter data, int factor){
    return IteratorMultiplier<Iter>{data, factor};
}


// strided iterator

template<class Iter>
struct StridedIterator{
    using value_type = typename std::iterator_traits<Iter>::value_type;
    using reference = typename std::iterator_traits<Iter>::reference;

    Iter data;
    int stride;

    HOSTDEVICEQUALIFIER
    StridedIterator(Iter data_, int stride_) 
        : data(data_), stride(stride_){
    }

    HOSTDEVICEQUALIFIER
    reference operator*(){
        return *data;
    }
    
    HOSTDEVICEQUALIFIER
    const reference operator*() const{
        return *data;
    }

    HOSTDEVICEQUALIFIER
    StridedIterator operator+(int i) const{
        return StridedIterator{data + i * stride, stride};
    }
};

#endif