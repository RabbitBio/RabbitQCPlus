#ifndef CARE_MEMORY_MANAGEMENT_HPP
#define CARE_MEMORY_MANAGEMENT_HPP

#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <string>
#include <cstdint>
#include <cassert>
#include <fstream>
#include <limits>
#include <map>
#include <iostream>

    struct MemoryUsage{
        std::size_t host = 0;
        std::map<int, std::size_t> device{};

        MemoryUsage operator+(const MemoryUsage rhs) const{
            MemoryUsage res{};
            res.host = host + rhs.host;
            res.device = device;
            for(auto pair : rhs.device){
                res.device[pair.first] += pair.second;
            }

            return res;
        }

        MemoryUsage& operator+=(const MemoryUsage rhs){
            MemoryUsage m = operator+(rhs);
            *this = m;
            return *this;
        }
    };





    __inline__
    std::size_t getAvailableMemoryInKB_linux(){
        //https://stackoverflow.com/questions/349889/how-do-you-determine-the-amount-of-linux-system-ram-in-c
        std::string token;
        std::ifstream file("/proc/meminfo");
        assert(bool(file));
        while(file >> token) {
            if(token == "MemAvailable:") {
                std::size_t mem;
                if(file >> mem) {
                    return mem;
                } else {
                    return 0;       
                }
            }
            file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        return 0;
    }

    __inline__ 
    std::size_t getCurrentRSS_linux(){
        std::ifstream in("/proc/self/statm");
        std::size_t tmp, rss;
        in >> tmp >> rss;
        
        return rss * sysconf(_SC_PAGESIZE);
    }

    __inline__
    std::size_t getRSSLimit_linux(){
        rlimit rlim;
        int ret = getrlimit(RLIMIT_RSS, &rlim);
        if(ret != 0){
            std::perror("Could not get RSS limit!");
            return 0;
        }
        return rlim.rlim_cur;    
    }


    __inline__
    std::size_t getAvailableMemoryInKB(){
        //return getAvailableMemoryInKB_linux();

        return std::min(getAvailableMemoryInKB_linux(), (getRSSLimit_linux() - getCurrentRSS_linux()) / 1024);
    }

    __inline__ 
    std::size_t getMaxRSSUsageInKB(){
        struct rusage usage;

        int ret = getrusage(RUSAGE_SELF, &usage);
        if(ret != 0){
            std::perror("Could not get max RSS size!");
            return 0;
        }

        return usage.ru_maxrss;
    }

    __inline__
    void compareMaxRssToLimit(std::size_t memoryLimitBytes, const std::string& errorMessage){
        std::size_t maxRssKB = getMaxRSSUsageInKB();
        std::size_t maxRssBytes = maxRssKB * 1024;

        if(maxRssBytes > memoryLimitBytes){
            std::cerr << errorMessage << ": maxRssBytes = " << maxRssBytes << ", memoryLimitBytes = " << memoryLimitBytes << "\n";
        }
        else{
            std::cerr << "memorylimit ok. maxRssBytes = " << maxRssBytes << ", memoryLimitBytes = " << memoryLimitBytes << "\n";
        }
    }



#endif