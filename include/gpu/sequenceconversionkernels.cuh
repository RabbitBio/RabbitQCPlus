#ifndef CARE_SEQUENCE_CONVERSION_KERNELS_CUH
#define CARE_SEQUENCE_CONVERSION_KERNELS_CUH

#include <sequencehelpers.hpp>

namespace care{
namespace gpu{

template<class Group, class InputIndexTrafo, class OutputIndexTrafo>
__device__
void convert2BitTo2BitHiLo(
    Group& group,
    unsigned int* out,
    const unsigned int* in,
    int length,
    InputIndexTrafo inindextrafo,
    OutputIndexTrafo outindextrafo
){
    if(length <= 0) return;

    const int inInts = SequenceHelpers::getEncodedNumInts2Bit(length);
    const int outInts = SequenceHelpers::getEncodedNumInts2BitHiLo(length);

    for(int i = group.thread_rank(); i < outInts; i += group.size()){
        const int outIndex = outindextrafo(i);
        const int inindex1 = inindextrafo(((i % (outInts / 2))*2));
        const unsigned int data1 = in[inindex1];
        const int indexInHalf = i % (outInts / 2);
        const int whichHalf = (i < outInts / 2) ? 0 : 1; 

        const unsigned int bits16_1 = SequenceHelpers::extractEvenBits(data1 >> (1 - whichHalf));
        unsigned int result = bits16_1 << 16;
        if((indexInHalf < outInts / 2 - 1) || ((length-1) % 32) >= 16){
            const int inindex2 = inindex1 + inindextrafo(1);
            const unsigned int data2 = in[inindex2];
            const unsigned int bits16_2 = SequenceHelpers::extractEvenBits(data2 >> (1 - whichHalf));
            result = result | bits16_2;
        }
        out[outIndex] = result;
    }
};



void callConversionKernel2BitTo2BitHiLoNN(
            const unsigned int* d_inputdata,
            size_t inputpitchInInts,
            unsigned int* d_outputdata,
            size_t outputpitchInInts,
            const int* d_sequenceLengths,
            int numSequences,
            cudaStream_t stream);

void callConversionKernel2BitTo2BitHiLoNT(
            const unsigned int* d_inputdata,
            size_t inputpitchInInts,
            unsigned int* d_outputdata,
            size_t outputpitchInInts,
            const int* d_sequenceLengths,
            int numSequences,
            cudaStream_t stream);

void callConversionKernel2BitTo2BitHiLoTT(
            const unsigned int* d_inputdata,
            size_t inputpitchInInts,
            unsigned int* d_outputdata,
            size_t outputpitchInInts,
            const int* d_sequenceLengths,
            int numSequences,
            cudaStream_t stream);            

void callEncodeSequencesTo2BitKernel(
    unsigned int* d_encodedSequences,
    const char* d_decodedSequences,
    const int* d_sequenceLengths,
    size_t decodedSequencePitchInBytes,
    size_t encodedSequencePitchInInts,
    int numSequences,
    int groupsize,
    cudaStream_t stream
);

void callDecodeSequencesFrom2BitKernel(
    char* d_decodedSequences,
    const unsigned int* d_encodedSequences,
    const int* d_sequenceLengths,
    size_t decodedSequencePitchInBytes,
    size_t encodedSequencePitchInInts,
    int numSequences,
    int groupsize,
    cudaStream_t stream
);

}
}


#endif