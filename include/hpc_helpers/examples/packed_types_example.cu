#include <iostream>
#include "../include/packed_types.cuh"

int main()
{
    using namespace packed_types;

    // packed type with three fields: 13 bits, 8 bits, 34 bits, (+ 9 bits padding)
    using pack_t = PackedTriple<13, 8, 34>;
    // or using pack_t = PackedPair<_, _>;
    // or using pack_t = PackedQuadruple<_, _, _, _>;

    // payload of packed data type
    using base_t = pack_t::base_type;

    std::cout
        << "size of pack_t="
        << sizeof(base_t)
        << " is equal to the size of base_t="
        << sizeof(base_t)
        << std::endl;

    // print data layout
    std::cout
        << "bit partition: MSB->{padding="
        << +pack_t::padding_bits() << "bits}[third="
        << +pack_t::third_bits() << "bits][second="
        << +pack_t::second_bits() << "][first="
        << +pack_t::first_bits() << "bits]<-LSB"
        << std::endl;

    // build packed triple
    pack_t triple{1234, 12, 123};

    // get fields
    std::cout
        << "triple = ("
        << triple.first()
        << ", "
        << triple.second()
        << ", "
        << triple.third()
        << ")"
        << std::endl;

    // also supports std::get<>
    std::cout
        << "std::get<> triple = ("
        << get<0>(triple)
        << ", "
        << get<1>(triple)
        << ", "
        << get<2>(triple)
        << ")"
        << std::endl;

    // and member functions with similar syntax i.e. get<>()/set<>()
    // NOTE you may need to provide additional info to the compiler
    // e.g. "triple.template get<0>()"
    std::cout
        << "triple.get<> triple = ("
        << triple.get<0>()
        << ", "
        << triple.get<1>()
        << ", "
        << triple.get<2>()
        << ")"
        << std::endl;

    // update third field
    triple.third(42);
    // or triple.set<2>(42) or triple.template set<2>(42)
    std::cout << "third = " << triple.third() << std::endl;

    // update third field with using a float as input
    // implicitly reinterpretates bits as base type
    float pi = 3.14;
    triple.set<2>(pi);
    // get float from pack
    float pi_too = triple.third_as<float>(); // or triple.get<2, float>()
    std::cout << pi << " == " << pi_too << std::endl;

    // following line should trigger an assertion error since 12345 needs more than 8 bit
    // triple.second(12345);
    std::cout
        << std::boolalpha
        << "should be false: "
        << pack_t::is_valid_second(12345) // or is_valid<1>(12345)
        << std::endl;

    // support for atomic updates:
    // comes with specializations for CUDA's atomicCAS() and atomicExch()
    // also valid: std::atomic<pack_t>
}
