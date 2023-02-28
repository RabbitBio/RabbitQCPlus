#ifndef CUDAGRAPHHELPERS_CUH
#define CUDAGRAPHHELPERS_CUH

#include <gpu/cudaerrorcheck.cuh>

#include <utility>
#include <memory>
#include <algorithm>

struct CudaGraph{
    bool valid = false;
    cudaGraphExec_t execgraph = nullptr;

    CudaGraph() = default;
    
    CudaGraph(const CudaGraph&) = delete;
    CudaGraph& operator=(const CudaGraph&) = delete;

    CudaGraph(CudaGraph&& rhs){
        *this = std::move(rhs);
    }

    CudaGraph& operator=(CudaGraph&& rhs){
        destroy();

        valid = std::exchange(rhs.valid, false);
        execgraph = std::exchange(rhs.execgraph, nullptr);

        return *this;
    }

    ~CudaGraph(){
        destroy();
    }

    void destroy(){
        if(execgraph != nullptr){
            CUDACHECK(cudaGraphExecDestroy(execgraph));
        }
        execgraph = nullptr;
        valid = false;
    }

    template<class Func>
    void capture(Func&& func){
        if(execgraph != nullptr){
            CUDACHECK(cudaGraphExecDestroy(execgraph));
        }

        cudaStream_t stream;
        CUDACHECK(cudaStreamCreate(&stream));
        
        CUDACHECK(cudaStreamBeginCapture(stream, cudaStreamCaptureModeRelaxed));

        func(stream);

        cudaGraph_t graph;
        CUDACHECK(cudaStreamEndCapture(stream, &graph));
        
        cudaGraphExec_t execGraph;
        cudaGraphNode_t errorNode;
        auto logBuffer = std::make_unique<char[]>(1025);
        std::fill_n(logBuffer.get(), 1025, 0);
        cudaError_t status = cudaGraphInstantiate(&execGraph, graph, &errorNode, logBuffer.get(), 1025);
        if(status != cudaSuccess){
            if(logBuffer[1024] != '\0'){
                std::cerr << "cudaGraphInstantiate: truncated error message: ";
                std::copy_n(logBuffer.get(), 1025, std::ostream_iterator<char>(std::cerr, ""));
                std::cerr << "\n";
            }else{
                std::cerr << "cudaGraphInstantiate: error message: ";
                std::cerr << logBuffer.get();
                std::cerr << "\n";
            }            
        }

        CUDACHECK(status);

        CUDACHECK(cudaGraphDestroy(graph));

        execgraph = execGraph;
        valid = true;

        CUDACHECK(cudaStreamDestroy(stream));
    }

    void execute(cudaStream_t stream){
        assert(valid);
        
        CUDACHECK(cudaGraphLaunch(execgraph, stream));
    }
};

#endif