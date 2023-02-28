#ifndef CARE_FOREST_HPP
#define CARE_FOREST_HPP


#include <fstream>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <array>
#include <deserialize.hpp>

#include <queue>

namespace care {

namespace gpu{
    class GpuForest;
}

template<typename extractor_t>
class ForestClf {
    friend class care::gpu::GpuForest;

    struct Node {
        uint8_t att;
        uint8_t flag;
        double split;
        union {
            uint32_t idx;
            float prob; 
        } lhs, rhs;
    };

    using Tree = std::vector<Node>;
    using features_t = typename extractor_t::features_t;

    void populate(std::ifstream& is, Tree& tree) {
        Node& node = *tree.emplace(tree.end());
        read_one(is, node.att);
        read_one(is, node.split);
        read_one(is, node.flag);
        if (node.flag / 2)
            new(&node.lhs.prob) float(read_one<float>(is));
        else {
            new(&node.lhs.idx) uint32_t(tree.size());
            populate(is, tree);
        }
        if (node.flag % 2)
            new(&node.rhs.prob) float(read_one<float>(is));
        else {
            new(&node.rhs.idx) uint32_t(tree.size());
            populate(is, tree);
        }
    }

    float prob(const features_t& features, const Tree& tree, size_t i = 0) const {
        if (features[tree[i].att] <= tree[i].split) {
            if (tree[i].flag / 2)
                return tree[i].lhs.prob;
            return prob(features, tree, tree[i].lhs.idx);
        } else {
            if (tree[i].flag % 2)
                return tree[i].rhs.prob;
            return prob(features, tree, tree[i].rhs.idx);
        }
    }

    using Forest = std::vector<Tree>;
    
    Forest forest_;
    float thresh_;

    void sortNodesViaBFS(Tree& tree){
        if(tree.size() > 0){
            std::vector<bool> visited(tree.size(), false);
            std::vector<int> indices;
            std::queue<int> queue;

            indices.reserve(tree.size());

            auto pushNode = [&](int i){
                if(!visited[i]){
                    queue.push(i);                
                    visited[i] = true;
                }
            };

            pushNode(0);

            while(!queue.empty()){
                int nodeIndex = queue.front();
                queue.pop();
                indices.push_back(nodeIndex);

                const auto& node = tree[nodeIndex];

                if(node.flag / 2 == 0){
                    pushNode(node.lhs.idx);
                }

                if(node.flag % 2 == 0){
                    pushNode(node.rhs.idx);
                }
            }

            // assert(indices.size() == tree.size());
            // auto indicescopy = indices;
            // std::sort(indicescopy.begin(), indicescopy.end());
            // auto iter = std::unique(indicescopy.begin(), indicescopy.end());
            // assert(iter == indicescopy.end());

            std::map<int, int> mapOldToNew;
            for(size_t i = 0; i < indices.size(); i++){
                mapOldToNew[indices[i]] = i;
            }

            Tree newtree(tree.size());
            for(size_t i = 0; i < indices.size(); i++){
                newtree[i] = tree[indices[i]];

                if(newtree[i].flag / 2 == 0){
                    newtree[i].lhs.idx = mapOldToNew[newtree[i].lhs.idx];
                }

                if(newtree[i].flag % 2 == 0){
                    newtree[i].rhs.idx = mapOldToNew[newtree[i].rhs.idx];
                }
            }

            //std::swap(newtree, tree);

            // for(int i = 0; i < 32; i++){
            //     std::cerr << indices[i] << "\n";
            // }
        }
    }

public:

    ForestClf (const std::string& path, uint32_t max_trees, float t = 0.5f) : 
        thresh_(t) 
    {
        std::ifstream is(path, std::ios::binary);
        if (!is)
            throw std::runtime_error("Loading classifier file failed! " + path);
        
        auto desc = read_str(is);
        auto expected = std::string(extractor_t());
        if (desc != expected)
            throw std::runtime_error("Classifier and extractor descriptors do not match! Expected: " + expected + " Received: " + desc);

        const auto n_trees = read_one<uint32_t>(is);
        max_trees = std::min(max_trees, n_trees);
        // std::cerr << "Using " << max_trees << " of " << n_trees << "trees.\n";
        forest_ = Forest(max_trees);
        for (Tree& tree: forest_) {
            tree.reserve(read_one<uint32_t>(is)); // reserve space for nodes
            populate(is, tree);

            //sortNodesViaBFS(tree);
        }

        //sortNodesViaBFS(forest_[0]);

        // std::ofstream os("tempforestout.txt");
        // const auto& tree = forest_[0];
        // size_t n = tree.size();
        // for(size_t i = 0; i < n; i++){
        //     os << i;
        //     if(tree[i].flag / 2 == 0){
        //         os << " " << tree[i].lhs.idx;
        //     }else{
        //         os << " -1";
        //     }

        //     if(tree[i].flag % 2 == 0){
        //         os << " " << tree[i].rhs.idx;
        //     }else{
        //         os << " -1";
        //     }
        //     os << "\n";
        // }
    }

    void threshold(float t) {
        thresh_ = t;
    }

    float threshold() const {
        return thresh_;
    }

    float prob_full(const features_t& features) const {
        float sum = 0.f;
        for (const Tree& tree: forest_)
            sum += prob(features, tree);
        return sum/forest_.size();
    }

    float prob_full_debug(const features_t& features) const {
        float sum = 0.f;
        std::cerr << "# ";
        for (const Tree& tree: forest_) {
            auto tmp = prob(features, tree);
            sum += tmp;
            std::cerr << tmp << " ";
        }
        std::cerr << "# ";
        return sum/forest_.size();
    }

    bool decide_full(const features_t& features) const {
        return prob_full(features) >= thresh_;
    }

    bool decide(const features_t& features) const {
        float sum = 0.f;
        for (const Tree& tree: forest_){
            sum += prob(features, tree);
            if(sum / forest_.size() >= thresh_){
                return true;
            }
        }
        return false;
    }
};




} // namespace care


#endif
