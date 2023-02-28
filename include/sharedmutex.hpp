#ifndef CARE_SHARED_MUTEX_HPP
#define CARE_SHARED_MUTEX_HPP


#include <shared_mutex>

#if __cplusplus >= 201703L

using SharedMutex = std::shared_mutex;

#elif __cplusplus == 201402L

using SharedMutex = std::shared_timed_mutex;

#endif



#endif