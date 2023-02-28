#ifndef CARE_FOREST_GPU_CUH
#define CARE_FOREST_GPU_CUH


#include <hpc_helpers.cuh>
#include <classification.hpp>
#include <memorymanagement.hpp>
#include <gpu/cudaerrorcheck.cuh>

#include <utility>

#include <cub/cub.cuh>

//#include <forest.hpp>

#if __CUDACC_VER_MAJOR__ >= 11

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>

namespace cg = cooperative_groups;

#endif

namespace care{
namespace gpu{

struct GpuForest {

    using Features = anchor_extractor::features_t;

    union IndexOrProb{
        uint32_t idx;
        float prob; 
    };

    struct Node {
        uint8_t att;
        uint8_t flag;
        double split;
        union {
            uint32_t idx;
            float prob; 
        } lhs, rhs;
    };

    struct Clf{
        using NodeType = Node;

        int numTrees;
        const Node** data;

        template<class Iter>
        HOSTDEVICEQUALIFIER
        float decideTree(Iter features, const Node* tree, size_t i = 0) const {
            if (*(features + tree[i].att) <= tree[i].split) {
                if (tree[i].flag / 2)
                    return tree[i].lhs.prob;
                return decideTree(features, tree, tree[i].lhs.idx);
            } else {
                if (tree[i].flag % 2)
                    return tree[i].rhs.prob;
                return decideTree(features, tree, tree[i].rhs.idx);
            }
        }

        template<class Iter>
        HOSTDEVICEQUALIFIER
        bool decide(Iter features, float thresh) const {
            float prob = 0.f;
            for(int t = 0; t < numTrees; t++){
                prob += decideTree(features, data[t]);
                if(prob / numTrees >= thresh){
                    return true;
                }
            }
            return false;
        }

        //use group to parallelize over trees. one thread per tree
        //correct result is returned by the first thread in the group
        template<class Group, class Iter, class GroupReduceFloatSum>
        DEVICEQUALIFIER
        float decide(Group& g, Iter features, float thresh, GroupReduceFloatSum reduce) const{
            #if 1
            const int limit = SDIV(numTrees, g.size()) * g.size();
            float accumulatedProb = 0.f;

            for(int t = g.thread_rank(); t < limit; t += g.size()){
                float currentProb = 0.0;
                if(t < numTrees){
                    currentProb += decideTree(features, data[t]);
                }
                accumulatedProb += reduce(currentProb);
                if(accumulatedProb / numTrees >= thresh){
                    return true;
                }
            }
            return false;
            #else

            float prob = 0.f;
            for(int t = g.thread_rank(); t < numTrees; t += g.size()){
                prob += decideTree(features, data[t]);
            }
            prob = reduce(prob);
            return prob / numTrees >= thresh;
            #endif
        }
        
    };

public:
    int deviceId{};
    int numTrees{};
    Node** d_data{};
    std::vector<std::size_t> numNodesPerTree{};

    MemoryUsage getMemoryInfo() const noexcept{
        MemoryUsage result{};

        for(auto num : numNodesPerTree){
            result.device[deviceId] += sizeof(Node) * num;
        }

        return result;
    }


    GpuForest() = default;

    template<class CpuForest>
    GpuForest(const CpuForest& clf, int gpuId)
        : deviceId(gpuId) {

        cub::SwitchDevice sd{deviceId};

        numTrees = clf.forest_.size();

        std::vector<Node*> devicePointers(numTrees);

        for(int t = 0; t < numTrees; t++){
            const int numNodes = clf.forest_[t].size();

            numNodesPerTree.emplace_back(numNodes);

            std::vector<Node> nodes(numNodes);

            for(int n = 0; n < numNodes; n++){
                const auto forestNode = clf.forest_[t][n];

                Node node;
                node.att = forestNode.att;
                node.flag = forestNode.flag;
                node.split = forestNode.split;
                node.lhs.idx = forestNode.lhs.idx;
                node.rhs.idx = forestNode.rhs.idx;

                nodes[n] = node;
            }

            CUDACHECK(cudaMalloc(&devicePointers[t], sizeof(Node) * numNodes));
            CUDACHECK(cudaMemcpy(
                devicePointers[t],
                nodes.data(),
                sizeof(Node) * numNodes,
                H2D
            ));
        }

        CUDACHECK(cudaMalloc(&d_data, sizeof(Node*) * numTrees));

        CUDACHECK(cudaMemcpy(
            d_data,
            devicePointers.data(),
            sizeof(Node*) * numTrees,
            H2D
        ));
    }

    GpuForest(GpuForest&& rhs){
        deviceId = std::exchange(rhs.deviceId, 0);
        numTrees = std::exchange(rhs.numTrees, 0);
        d_data = std::exchange(rhs.d_data, nullptr);
        numNodesPerTree = std::exchange(rhs.numNodesPerTree, std::vector<uint64_t>{});
    }

    GpuForest(const GpuForest& rhs) = delete;

    GpuForest& operator=(GpuForest&& rhs){
        deviceId = std::exchange(rhs.deviceId, 0);
        numTrees = std::exchange(rhs.numTrees, 0);
        d_data = std::exchange(rhs.d_data, nullptr);
        numNodesPerTree = std::exchange(rhs.numNodesPerTree, std::vector<uint64_t>{});
        return *this;
    }

    GpuForest& operator=(const GpuForest& rhs) = delete;

    ~GpuForest(){
        cub::SwitchDevice sd{deviceId};

        std::vector<Node*> devicePointers(numTrees);
        CUDACHECK(cudaMemcpy(
            devicePointers.data(),
            d_data,
            sizeof(Node*) * numTrees,
            D2H
        ));

        for(int t = 0; t < numTrees; t++){
            CUDACHECK(cudaFree(devicePointers[t]));
        }
        CUDACHECK(cudaFree(d_data));
    }

    operator Clf() const{
        return getClf();
    }

    Clf getClf() const{
        return Clf{numTrees, (const Node**)d_data};
    }

};



}
}


#endif