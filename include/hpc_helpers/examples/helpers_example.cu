#include <iostream>
#include <vector>
#include <chrono>
#include <thread>

#include "../include/cuda_helpers.cuh"
#include "../include/hpc_helpers.h"
#include "../include/io_helpers.h"
#include "../include/timers.cuh"

#include "../include/cuda_raiiwrappers.cuh"
#include "../include/coop_group_helpers.cuh"

#include <cooperative_groups.h>

namespace cg = cooperative_groups;

#define N ((1L)<<(28))

GLOBALQUALIFIER
void reverse_kernel(int * array, size_t n) {

    size_t thid = helpers::global_thread_id();

    if (thid < n/2) {
        const int lower = array[thid];
        const int upper = array[n-thid-1];
        array[thid] = upper;
        array[n-thid-1] = lower;
    }
}

template<class Group>
HOSTDEVICEQUALIFIER
void reduce(Group& group, int* data, int n){
    int val = 0;
    for(int i = group.thread_rank(); i < n; i += group.size()){
        val += data[i];
    }

    for (int i = group.size() / 2; i > 0; i /= 2) {
        val += group.shfl_down(val, i);
    }

    group.sync();


    if(group.thread_rank() == 0){
        printf("reduced value: %d\n", val);
    }
}

void useReusableTimer(){

    {
        using namespace std::chrono_literals;

        helpers::CpuTimer t("ttt");

        t.start();

        std::this_thread::sleep_for(1s);

        t.stop();

        t.start();

        std::this_thread::sleep_for(1s);

        t.stop();

        t.print();
    }

    helpers::CpuTimer timer("reusable");

    std::vector<int> vec(1024*1024);

    for(size_t i = 0; i < vec.size(); i++){
        vec[i] = i;
    }

    timer.stop();

    timer.print();

    timer.start();

    helpers::CpuTimer timer2("oneshot");

    size_t sum = 0;

    for(size_t i = 0; i < vec.size(); i++){
        sum += vec[i];
    }

    std::cout << "sum = " << sum << "\n";

    timer.stop();
    timer.print();

    timer2.stop();
    timer2.print();

}


void useReusableGpuTimer(){

    {
        using namespace std::chrono_literals;

        helpers::GpuTimer t("reuseable gputimer");

        t.start();

        helpers::lambda_kernel<<<65535, 128>>>(
            [=] DEVICEQUALIFIER {

            });

        t.stop();

        t.print();

        t.start();

        helpers::lambda_kernel<<<65535, 128>>>(
            [=] DEVICEQUALIFIER {

            });

        t.stop();

        t.print();
    }

    helpers::CpuTimer timer("reusable");

    std::vector<int> vec(1024*1024);

    for(size_t i = 0; i < vec.size(); i++){
        vec[i] = i;
    }

    timer.stop();

    timer.print();

    timer.start();

    helpers::CpuTimer timer2("oneshot");

    size_t sum = 0;

    for(size_t i = 0; i < vec.size(); i++){
        sum += vec[i];
    }

    std::cout << "sum = " << sum << "\n";

    timer.stop();
    timer.print();

    timer2.stop();
    timer2.print();

}

void benchmarkAtomicAggAdd(std::uint64_t num_threads = 1000000000)
{
    int * ctr_d = nullptr;
    cudaMalloc(&ctr_d, sizeof(int)); CUERR

    auto reset_ctr = [ctr_d] () { cudaMemset(ctr_d, 0, sizeof(int)); CUERR };

    helpers::GpuTimer naive("naive atomicAdd");
    helpers::GpuTimer aggregated("warp-aggregated atomicAdd");

    for(std::uint8_t i = 1; i < WARPSIZE; i <<= 1)
    {
        const float selectivity = WARPSIZE / (float(i) * WARPSIZE);
        std::cout << selectivity * 100.0 << "% selectivity:" << std::endl;

        reset_ctr();
        naive.start();
        helpers::lambda_kernel
        <<<SDIV(num_threads, MAXBLOCKSIZE), MAXBLOCKSIZE>>>(
            [=] DEVICEQUALIFIER {
                const auto tid = helpers::global_thread_id();

                if(tid >= num_threads) return;

                if(tid % i == 0)
                {
                    atomicAdd(ctr_d, int((tid % 32) + 1));
                }
            });
        naive.stop();
        naive.print();
        naive.reset(); CUERR

        reset_ctr();
        aggregated.start();
        helpers::lambda_kernel
        <<<SDIV(num_threads, MAXBLOCKSIZE), MAXBLOCKSIZE>>>(
            [=] DEVICEQUALIFIER {
                const auto tid = helpers::global_thread_id();

                if(tid >= num_threads) return;

                if(tid % i == 0)
                {
                    helpers::atomicAggAdd(ctr_d, int((tid % 32) + 1));
                }
            });
        aggregated.stop();
        aggregated.print();
        aggregated.reset(); CUERR
    }

    cudaFree(ctr_d); CUERR
}

int main () {
    using namespace helpers;

    useReusableTimer();

    useReusableGpuTimer();

    CpuTimer allovertimer("allover");

    init_cuda_context();                                                  CUERR

    debug_printf("this message will only be shown if NDEBUG is undefined");

    std::vector<int> hostVector(N);
    for (size_t i = 0; i < N; i++)
        hostVector[i] = i;

    std::cout << "Host reduction\n";

    SingleThreadGroup hostgroup;
    reduce(hostgroup, hostVector.data(), 10);


    int * deviceArray = NULL;
    cudaMalloc(&deviceArray, sizeof(int)*N);                              CUERR
    cudaMemcpy(deviceArray, hostVector.data(), sizeof(int)*N, H2D);       CUERR

    std::cout << "Device reduction\n";
    lambda_kernel<<<1, 32>>>(
        [=] DEVICEQUALIFIER {
            auto tile = cg::tiled_partition<32>(cg::this_thread_block());

            reduce(tile, deviceArray, 10);
        });

    cudaDeviceSynchronize(); CUERR


    {
        GpuTimer gtimer("kernel");

        reverse_kernel<<<SDIV(N, MAXBLOCKSIZE), MAXBLOCKSIZE>>>(deviceArray, N);
                                                                          CUERR

        gtimer.stop();
        gtimer.print();
    }

    CudaStream stream;

    {
        CpuTimer timer2("scopedcputimer");
        GpuTimer gtimer2(stream, "lambda_kernel", std::cerr);

        lambda_kernel<<<SDIV(N, MAXBLOCKSIZE), MAXBLOCKSIZE, 0, stream>>>(
            [=] DEVICEQUALIFIER {
                size_t thid = global_thread_id();

                if (thid < N/2) {
                    const int lower = deviceArray[thid];
                    const int upper = deviceArray[N-thid-1];
                    deviceArray[thid] = upper;
                    deviceArray[N-thid-1] = lower;
                }
            });                                                               CUERR

        gtimer2.print();

        stream.synchronize();

        GpuTimer gthroughput("copy");

        cudaMemcpy(hostVector.data(), deviceArray, sizeof(int)*N, D2H);       CUERR

        gthroughput.print_throughput(sizeof(int), N);

        cudaFree(deviceArray);                                                CUERR

        allovertimer.stop();

    }

    allovertimer.print();

    // store data
    const std::vector<std::uint32_t> vectorToStore{1, 2, 3, 4, 42, 6};
    const std::string fname{"test_file.dump"};
    dump_binary(vectorToStore, fname);
    // load stored data
    const auto load = load_binary<std::uint32_t>(fname);
    for(const auto& x : load)
        std::cout << x << " ";
    std::cout << std::endl;

    std::cout
        << "available GPU memory: " << B2GB(available_gpu_memory()) << " GB\n";

    std::cout << "naive vs. warp-aggregated atomic add:" << std::endl;
    benchmarkAtomicAggAdd();

    std::cout << "causing memory error by allocating 2^60 bytes" << std::endl;
    cudaMalloc(&deviceArray, (1L<<60));                                   CUERR
    cudaDeviceSynchronize();

}
