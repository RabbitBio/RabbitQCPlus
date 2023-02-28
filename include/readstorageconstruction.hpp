#ifndef CARE_READSTORAGECONSTRUCTION_HPP 
#define CARE_READSTORAGECONSTRUCTION_HPP

#include <util.hpp>
#include <config.hpp>
#include <threadpool.hpp>
#include <readlibraryio.hpp>

#include <vector>
#include <array>
#include <mutex>
#include <condition_variable>
#include <algorithm>
#include <string>
#include <iostream>

namespace care{


template<class MakeInserterFunc, class MakeReadContainsNFunc>
void constructReadStorageFromFiles(
    const std::vector<std::string>& inputfiles,
    bool /*useQualityScores*/,
    read_number expectedNumberOfReads,
    int /*expectedMinimumReadLength*/,
    int expectedMaximumReadLength,
    int threads,
    bool showProgress,
    MakeInserterFunc makeInserterFunc,
    MakeReadContainsNFunc makeReadContainsNFunc
){
            
    auto setReadContainsN = makeReadContainsNFunc();

    auto checkRead = [&](read_number readIndex, Read& read, int& Ncount){
        const int readLength = int(read.sequence.size());

        if(readIndex >= expectedNumberOfReads){
            throw std::runtime_error("Error! Expected " + std::to_string(expectedNumberOfReads)
                                    + " reads, but file contains at least "
                                    + std::to_string(readIndex+1) + " reads.");
        }

        if(readLength > expectedMaximumReadLength){
            throw std::runtime_error("Error! Expected maximum read length = "
                                    + std::to_string(expectedMaximumReadLength)
                                    + ", but read " + std::to_string(readIndex)
                                    + "has length " + std::to_string(readLength));
        }

        auto isValidBase = [](char c){
            constexpr std::array<char, 10> validBases{'A','C','G','T','a','c','g','t'};
            return validBases.end() != std::find(validBases.begin(), validBases.end(), c);
        };

        const int undeterminedBasesInRead = std::count_if(read.sequence.begin(), read.sequence.end(), [&](char c){
            return !isValidBase(c);
        });

        //nmap[undeterminedBasesInRead]++;

        if(undeterminedBasesInRead > 0){
            setReadContainsN(readIndex, true);
        }

        constexpr std::array<char, 4> bases = {'A', 'C', 'G', 'T'};

        for(auto& c : read.sequence){
            if(c == 'a') c = 'A';
            else if(c == 'c') c = 'C';
            else if(c == 'g') c = 'G';
            else if(c == 't') c = 'T';
            else if(!isValidBase(c)){
                c = bases[Ncount];
                Ncount = (Ncount + 1) % 4;
            }
        }
    };

    constexpr int maxbuffersize = 65536;
    constexpr int numBuffers = 3;

    std::array<std::vector<read_number>, numBuffers> indicesBuffers;
    std::array<std::vector<Read>, numBuffers> readsBuffers;
    std::array<int, numBuffers> bufferSizes;
    std::array<bool, numBuffers> canBeUsed;
    std::array<std::mutex, numBuffers> mutex;
    std::array<std::condition_variable, numBuffers> cv;

    using Inserter_t = decltype(makeInserterFunc());

    std::array<Inserter_t, 3> inserterFuncs{makeInserterFunc(), makeInserterFunc(), makeInserterFunc()};

    ThreadPool threadPool(numBuffers);
    ThreadPool pforPool(std::max(1, threads - numBuffers));

    for(int i = 0; i < numBuffers; i++){
        indicesBuffers[i].resize(maxbuffersize);
        readsBuffers[i].resize(maxbuffersize);
        canBeUsed[i] = true;
        bufferSizes[i] = 0;

    }

    int bufferindex = 0;
    read_number globalReadId = 0;

    auto showProgressFunc = [show = showProgress](auto totalCount, auto seconds){
        if(show){
            std::cout << "Processed " << totalCount << " reads in file. Elapsed time: " 
                            << seconds << " seconds." << std::endl;
        }
    };

    auto updateShowProgressInterval = [](auto duration){
        return duration * 2;
    };

    ProgressThread<read_number> progressThread(
        expectedNumberOfReads, 
        showProgressFunc, 
        updateShowProgressInterval
    );

    for(const auto& inputfile : inputfiles){
        std::cout << "Converting reads of file " << inputfile << ", storing them in memory\n";

        forEachReadInFile(inputfile,
                        [&](auto /*readnum*/, auto& read){

                if(!canBeUsed[bufferindex]){
                    std::unique_lock<std::mutex> ul(mutex[bufferindex]);
                    if(!canBeUsed[bufferindex]){
                        //std::cerr << "waiting for other buffer\n";
                        //nvtx::push_range("wait for doublebuffer", 0);
                        cv[bufferindex].wait(ul, [&](){ return canBeUsed[bufferindex]; });
                        //nvtx::pop_range();
                    }
                }

                auto indicesBufferPtr = &indicesBuffers[bufferindex];
                auto readsBufferPtr = &readsBuffers[bufferindex];
                auto& buffersize = bufferSizes[bufferindex];

                (*indicesBufferPtr)[buffersize] = globalReadId;
                std::swap((*readsBufferPtr)[buffersize], read);

                ++buffersize;

                ++globalReadId;

                progressThread.addProgress(1);

                if(buffersize >= maxbuffersize){
                    canBeUsed[bufferindex] = false;

                    //std::cerr << "launch other thread\n";
                    //nvtx::push_range("enqeue", 0);
                    threadPool.enqueue([&, indicesBufferPtr, readsBufferPtr, bufferindex](){
                        //std::cerr << "buffer " << bufferindex << " running\n";
                        const int buffersize = bufferSizes[bufferindex];

                        //nvtx::push_range("process read batch", 0);
                        int nmodcounter = 0;

                        //nvtx::push_range("check read batch", 0);
                        for(int i = 0; i < buffersize; i++){
                            read_number readId = (*indicesBufferPtr)[i];
                            auto& read = (*readsBufferPtr)[i];
                            checkRead(readId, read, nmodcounter);
                        }
                        //nvtx::pop_range();

                        //nvtx::push_range("insert read batch", 0);
                        inserterFuncs[bufferindex](
                            &pforPool, 
                            indicesBufferPtr->data(), 
                            readsBufferPtr->data(), 
                            buffersize
                        );
                        //nvtx::pop_range();

                        //TIMERSTARTCPU(clear);
                        bufferSizes[bufferindex] = 0;
                        //TIMERSTOPCPU(clear);
                        
                        std::lock_guard<std::mutex> l(mutex[bufferindex]);
                        canBeUsed[bufferindex] = true;
                        cv[bufferindex].notify_one();

                        //nvtx::pop_range();

                        //std::cerr << "buffer " << bufferindex << " finished\n";
                    });

                    bufferindex = (bufferindex + 1) % numBuffers; //swap buffers

                    //nvtx::pop_range();
                }

        });
    }

    auto indicesBufferPtr = &indicesBuffers[bufferindex];
    auto readsBufferPtr = &readsBuffers[bufferindex];
    auto& buffersize = bufferSizes[bufferindex];

    if(buffersize > 0){
        if(!canBeUsed[bufferindex]){
            std::unique_lock<std::mutex> ul(mutex[bufferindex]);
            if(!canBeUsed[bufferindex]){
                //std::cerr << "waiting for other buffer\n";
                cv[bufferindex].wait(ul, [&](){ return canBeUsed[bufferindex]; });
            }
        }

        int nmodcounter = 0;

        for(int i = 0; i < buffersize; i++){
            read_number readId = (*indicesBufferPtr)[i];
            auto& read = (*readsBufferPtr)[i];
            checkRead(readId, read, nmodcounter);
        }

        inserterFuncs[bufferindex](&pforPool, indicesBufferPtr->data(), readsBufferPtr->data(), buffersize);

        buffersize = 0;
    }

    for(int i = 0; i < numBuffers; i++){
        std::unique_lock<std::mutex> ul(mutex[i]);
        if(!canBeUsed[i]){
            //std::cerr << "Reading file completed. Waiting for buffer " << i << "\n";
            cv[i].wait(ul, [&](){ return canBeUsed[i]; });
        }
    }

    progressThread.finished();

}







template<class MakeInserterFunc, class MakeReadContainsNFunc>
void constructReadStorageFromPairedEndFiles(
    const std::vector<std::string>& inputfiles,
    bool /*useQualityScores*/,
    read_number expectedNumberOfReads,
    int /*expectedMinimumReadLength*/,
    int expectedMaximumReadLength,
    int threads,
    bool showProgress,
    MakeInserterFunc makeInserterFunc,
    MakeReadContainsNFunc makeReadContainsNFunc
){
            
    auto setReadContainsN = makeReadContainsNFunc();

    auto checkRead = [&](read_number readIndex, Read& read, int& Ncount){
        const int readLength = int(read.sequence.size());

        if(readIndex >= expectedNumberOfReads){
            throw std::runtime_error("Error! Expected " + std::to_string(expectedNumberOfReads)
                                    + " reads, but file contains at least "
                                    + std::to_string(readIndex+1) + " reads.");
        }

        if(readLength > expectedMaximumReadLength){
            throw std::runtime_error("Error! Expected maximum read length = "
                                    + std::to_string(expectedMaximumReadLength)
                                    + ", but read " + std::to_string(readIndex)
                                    + "has length " + std::to_string(readLength));
        }

        auto isValidBase = [](char c){
            constexpr std::array<char, 10> validBases{'A','C','G','T','a','c','g','t'};
            return validBases.end() != std::find(validBases.begin(), validBases.end(), c);
        };

        const int undeterminedBasesInRead = std::count_if(read.sequence.begin(), read.sequence.end(), [&](char c){
            return !isValidBase(c);
        });

        //nmap[undeterminedBasesInRead]++;

        if(undeterminedBasesInRead > 0){
            setReadContainsN(readIndex, true);
        }

        constexpr std::array<char, 4> bases = {'A', 'C', 'G', 'T'};

        for(auto& c : read.sequence){
            if(c == 'a') c = 'A';
            else if(c == 'c') c = 'C';
            else if(c == 'g') c = 'G';
            else if(c == 't') c = 'T';
            else if(!isValidBase(c)){
                c = bases[Ncount];
                Ncount = (Ncount + 1) % 4;
            }
        }
    };

    constexpr int maxbuffersize = 65536;
    constexpr int numBuffers = 3;



    std::array<std::vector<read_number>, numBuffers> indicesBuffers;
    std::array<std::vector<Read>, numBuffers> readsBuffers;
    std::array<int, numBuffers> bufferSizes;
    std::array<bool, numBuffers> canBeUsed;
    std::array<std::mutex, numBuffers> mutex;
    std::array<std::condition_variable, numBuffers> cv;

    using Inserter_t = decltype(makeInserterFunc());

    std::array<Inserter_t, 3> inserterFuncs{makeInserterFunc(), makeInserterFunc(), makeInserterFunc()};

    ThreadPool threadPool(numBuffers);
    ThreadPool pforPool(std::max(1, threads - numBuffers));

    for(int i = 0; i < numBuffers; i++){
        indicesBuffers[i].resize(maxbuffersize);
        readsBuffers[i].resize(maxbuffersize);
        canBeUsed[i] = true;
        bufferSizes[i] = 0;

    }

    int bufferindex = 0;
    read_number globalReadId = 0;

    auto showProgressFunc = [show = showProgress](auto totalCount, auto seconds){
        if(show){
            std::cout << "Processed " << totalCount << " reads in file. Elapsed time: " 
                            << seconds << " seconds." << std::endl;
        }
    };

    auto updateShowProgressInterval = [](auto duration){
        return duration * 2;
    };

    ProgressThread<read_number> progressThread(
        expectedNumberOfReads, 
        showProgressFunc, 
        updateShowProgressInterval
    );

    auto work = [&](auto /*readnum*/, auto& read){

        if(!canBeUsed[bufferindex]){
            std::unique_lock<std::mutex> ul(mutex[bufferindex]);
            if(!canBeUsed[bufferindex]){
                //std::cerr << "waiting for other buffer\n";
                //nvtx::push_range("wait for doublebuffer", 0);
                cv[bufferindex].wait(ul, [&](){ return canBeUsed[bufferindex]; });
                //nvtx::pop_range();
            }
        }

        auto indicesBufferPtr = &indicesBuffers[bufferindex];
        auto readsBufferPtr = &readsBuffers[bufferindex];
        auto& buffersize = bufferSizes[bufferindex];

        (*indicesBufferPtr)[buffersize] = globalReadId;
        std::swap((*readsBufferPtr)[buffersize], read);

        ++buffersize;

        ++globalReadId;

        progressThread.addProgress(1);

        if(buffersize >= maxbuffersize){
            canBeUsed[bufferindex] = false;

            //std::cerr << "launch other thread\n";
            //nvtx::push_range("enqeue", 0);
            threadPool.enqueue([&, indicesBufferPtr, readsBufferPtr, bufferindex](){
                //std::cerr << "buffer " << bufferindex << " running\n";
                const int buffersize = bufferSizes[bufferindex];

                //nvtx::push_range("process read batch", 0);
                int nmodcounter = 0;

                //nvtx::push_range("check read batch", 0);
                for(int i = 0; i < buffersize; i++){
                    read_number readId = (*indicesBufferPtr)[i];
                    auto& read = (*readsBufferPtr)[i];
                    checkRead(readId, read, nmodcounter);
                }
                //nvtx::pop_range();

                //nvtx::push_range("insert read batch", 0);
                inserterFuncs[bufferindex](
                    &pforPool, 
                    indicesBufferPtr->data(), 
                    readsBufferPtr->data(), 
                    buffersize
                );
                //nvtx::pop_range();

                //TIMERSTARTCPU(clear);
                bufferSizes[bufferindex] = 0;
                //TIMERSTOPCPU(clear);
                
                std::lock_guard<std::mutex> l(mutex[bufferindex]);
                canBeUsed[bufferindex] = true;
                cv[bufferindex].notify_one();

                //nvtx::pop_range();

                //std::cerr << "buffer " << bufferindex << " finished\n";
            });

            bufferindex = (bufferindex + 1) % numBuffers; //swap buffers

            //nvtx::pop_range();
        }

    };

    /*
        accept either a single input file or two input files
    
        If two input files are present, it is assumed that the mate of i-th read in file 1
        is the i-th read in file 2. The first read in file 1 has id 0. The first read in file 2 has id 1.
        Read ids within a file increase by 2.


        If one input file is present, it is assumed that the i-th read pair consists of reads 2*i and 2*i+1
    */
    assert(inputfiles.size() > 0);
    assert(inputfiles.size() <= 2);

    if(inputfiles.size() == 1){

        const auto& filename1 = inputfiles[0];

        std::cout << "Converting paired reads of files " 
            << filename1  << ", storing them in memory\n";

        forEachReadInFile(filename1, work);
    
    }else{
        assert(inputfiles.size() == 2);

        const auto& filename1 = inputfiles[0];
        const auto& filename2 = inputfiles[1];

        std::cout << "Converting paired reads of files " 
            << filename1 << " and " << filename2 << ", storing them in memory\n";

        forEachReadInPairedFiles(filename1, filename2, work);
    
    }

    auto indicesBufferPtr = &indicesBuffers[bufferindex];
    auto readsBufferPtr = &readsBuffers[bufferindex];
    auto& buffersize = bufferSizes[bufferindex];

    if(buffersize > 0){
        if(!canBeUsed[bufferindex]){
            std::unique_lock<std::mutex> ul(mutex[bufferindex]);
            if(!canBeUsed[bufferindex]){
                //std::cerr << "waiting for other buffer\n";
                cv[bufferindex].wait(ul, [&](){ return canBeUsed[bufferindex]; });
            }
        }

        int nmodcounter = 0;

        for(int i = 0; i < buffersize; i++){
            read_number readId = (*indicesBufferPtr)[i];
            auto& read = (*readsBufferPtr)[i];
            checkRead(readId, read, nmodcounter);
        }

        inserterFuncs[bufferindex](&pforPool, indicesBufferPtr->data(), readsBufferPtr->data(), buffersize);

        buffersize = 0;
    }

    for(int i = 0; i < numBuffers; i++){
        std::unique_lock<std::mutex> ul(mutex[i]);
        if(!canBeUsed[i]){
            //std::cerr << "Reading file completed. Waiting for buffer " << i << "\n";
            cv[i].wait(ul, [&](){ return canBeUsed[i]; });
        }
    }

    progressThread.finished();

}



} //namespace care


#endif
