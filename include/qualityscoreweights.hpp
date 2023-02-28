#ifndef CARE_CPU_QUALITY_SCORE_WEIGHTS_HPP
#define CARE_CPU_QUALITY_SCORE_WEIGHTS_HPP

#include <config.hpp>

#include <array>
#include <cmath>

namespace care{
namespace cpu{

    struct QualityScoreConversion{
    public:
        QualityScoreConversion(){
            constexpr int size = 256;
            constexpr int ASCII_BASE = 33;
            constexpr float MIN_WEIGHT = 0.001f;

            for(int i = 0; i < size; i++){
                if(i < ASCII_BASE)
                    qscore_to_error_prob[i] = 1.0;
                else
                    qscore_to_error_prob[i] = std::pow(10.0f, -(i-ASCII_BASE)/10.0f);
            }

            for(int i = 0; i < size; i++){
                qscore_to_weight[i] = std::max(MIN_WEIGHT, 1.0f - qscore_to_error_prob[i]);
            }
        }

        float getErrorProbability(char c) const{
            return qscore_to_error_prob[static_cast<unsigned char>(c)];
        }

        float getWeight(char c) const{
            return qscore_to_weight[static_cast<unsigned char>(c)];
        }

    private:
        using Array_t = std::array<float, 256>;

        Array_t qscore_to_error_prob;
        Array_t qscore_to_weight;
    };

}
}

#endif
