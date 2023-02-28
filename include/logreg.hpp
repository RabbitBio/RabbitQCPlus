#ifndef CARE_LOGREG_HPP
#define CARE_LOGREG_HPP

#include <deserialize.hpp>
#include <array>

namespace care {

template<size_t num_features>
class LogRegClf {

    constexpr float logit(float p) const {
        return std::log(p/(1-p));
    }

    using weights_t = std::array<float, num_features>;

    weights_t w_;
    float b_;
    float thresh_;

public:

    LogRegClf (const std::string& path, float thresh = .5f) : 
        thresh_(logit(thresh)) 
    {
        std::ifstream is(path, std::ios::binary);
        if (read_one<uint32_t>(is) != num_features)
            throw std::runtime_error("Could not load LogRegClf! Feature shape does not match!");

        read_one<weights_t>(is, w_);
        read_one<float>(is, b_);
        is.close();
    }

    void threshold(float t) {
        thresh_ = logit(t);
    }

    float threshold() const {
        return thresh_;
    }

    template<typename features_t>
    bool decide(const features_t& features) const {
        float z = b_;
        for (size_t i = 0; i < num_features; ++i)
            z += w_[i] * features[i];
        return z >= thresh_;
    }
};

} // namespace care

#endif
