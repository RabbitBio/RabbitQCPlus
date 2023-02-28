#ifndef CARE_MINHASHERLIMIT_HPP
#define CARE_MINHASHERLIMIT_HPP

#include <config.hpp>

#include <limits>
#include <cmath>

inline
int calculateResultsPerMapThreshold(int coverage){
    int result = int(coverage * 2.5f);
    result = std::min(result, int(std::numeric_limits<care::BucketSize>::max()));
    result = std::max(10, result);
    return result;
}


#endif