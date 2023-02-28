#include "catch.hpp"
#include "packed_types.cuh"
#include "cuda_helpers.cuh"

using helpers::lambda_kernel;
using namespace packed_types;

TEMPLATE_TEST_CASE_SIG(
    "PackedPair with variable split",
    "[pack][pair][packedpair][variablesplit][template]",
    ((std::uint8_t FirstBits, std::uint8_t SecondBits),
        FirstBits, SecondBits),
        (16, 16),
        (15, 17),
        (18, 14),
        (7, 10),
        (32, 32),
        (31, 33),
        (34, 30),
        (10, 7))
{
    using pack_t = PackedPair<FirstBits, SecondBits>;
    using base_t = typename pack_t::base_type;

    REQUIRE(FirstBits + SecondBits <= sizeof(base_t) * CHAR_BIT);

    const base_t first_max = (base_t{1} << pack_t::first_bits()) - base_t{1};
    const base_t second_max = (base_t{1} << pack_t::second_bits()) - base_t{1};
    const base_t first  = GENERATE(as<base_t>{}, 0, 1, 2, 42);
    const base_t second  = GENERATE(as<base_t>{}, 0, 1, 2, 42);
    const base_t update = 60;

    CAPTURE(first, second, update, first_max, second_max);

    REQUIRE(first  <= first_max);
    REQUIRE(second <= second_max);
    REQUIRE(update <= first_max);
    REQUIRE(update <= second_max);

    CHECK(pack_t::is_valid_first(first));
    CHECK(pack_t::is_valid_second(second));
    CHECK(pack_t::is_valid_first(update));
    CHECK(pack_t::is_valid_second(update));

    SECTION("pack size")
    {
        CHECK(sizeof(pack_t) == sizeof(base_t));
    }

    SECTION("empty pack")
    {
        pack_t empty     = pack_t();
        pack_t empty_too = pack_t::empty();

        CHECK(empty == empty_too);
        CHECK(empty.first() == 0);
        CHECK(empty.second() == 0);
    }

    SECTION("set and get pack")
    {
        pack_t pack(first, second);

        CHECK(pack.first() == first);
        CHECK(pack.second() == second);

        CHECK(pack.template get<0>() == first);
        CHECK(pack.template get<1>() == second);

        CHECK(get<0>(pack) == first);
        CHECK(get<1>(pack) == second);

        SECTION("equality operator")
        {
            pack_t pack_too = pack;

            CHECK(pack == pack_too);
        }

        SECTION("update first")
        {
            pack.first(update);

            CHECK(pack.first() == update);
            CHECK(pack.second() == second);
        }

        SECTION("update second")
        {
            pack.template set<1>(update);

            CHECK(pack.first() == first);
            CHECK(pack.second() == update);
        }

        SECTION("maximum first value")
        {
            pack.first(first_max);

            CHECK(first_max != 0);
            CHECK(pack.first() == first_max);
            CHECK(pack.second() == second);
        }

        SECTION("maximum second value")
        {
            pack.second(second_max);

            CHECK(second_max != 0);
            CHECK(pack.first() == first);
            CHECK(pack.second() == second_max);
        }

        SECTION("atomic operations")
        {
            pack_t val = pack_t(update, update);

            pack_t * pack_d = nullptr;
            cudaMalloc(&pack_d, sizeof(pack_t));
            REQUIRE(cudaGetLastError() == cudaSuccess);

            SECTION("atomic CAS")
            {
                pack_t compare = pack_t(pack);

                cudaMemcpy(pack_d, &pack, sizeof(pack_t), H2D);
                REQUIRE(cudaGetLastError() == cudaSuccess);

                lambda_kernel
                <<<1, 1>>>([=] DEVICEQUALIFIER
                {
                    atomicCAS(pack_d, compare, val);
                });

                cudaMemcpy(&pack, pack_d, sizeof(pack_t), D2H);
                REQUIRE(cudaGetLastError() == cudaSuccess);

                CHECK(pack == val);
            }

            SECTION("atomic exchange")
            {
                cudaMemcpy(pack_d, &pack, sizeof(pack_t), H2D);
                REQUIRE(cudaGetLastError() == cudaSuccess);

                lambda_kernel
                <<<1, 1>>>([=] DEVICEQUALIFIER
                {
                    atomicExch(pack_d, val);
                });

                cudaMemcpy(&pack, pack_d, sizeof(pack_t), D2H);
                REQUIRE(cudaGetLastError() == cudaSuccess);

                CHECK(pack == val);
            }

            cudaFree(pack_d);
            CHECK(cudaGetLastError() == cudaSuccess);
        }

    }
}

TEMPLATE_TEST_CASE_SIG(
    "PackedTriple with variable split",
    "[pack][triple][packedtriple][variablesplit][template]",
    ((std::uint8_t FirstBits, std::uint8_t SecondBits, std::uint8_t ThirdBits),
        FirstBits, SecondBits, ThirdBits),
        (10, 10, 12),
        (8, 9, 15),
        (13, 8, 11),
        (7, 7, 7),
        (20, 20, 24),
        (18, 19, 27),
        (23, 19, 22))
{
    using pack_t = PackedTriple<FirstBits, SecondBits, ThirdBits>;
    using base_t = typename pack_t::base_type;

    REQUIRE(FirstBits + SecondBits + ThirdBits <= sizeof(base_t) * 8);

    const base_t first_max = (base_t{1} << pack_t::first_bits()) - base_t{1};
    const base_t second_max = (base_t{1} << pack_t::second_bits()) - base_t{1};
    const base_t third_max = (base_t{1} << pack_t::third_bits()) - base_t{1};
    const base_t first  = GENERATE(as<base_t>{}, 0, 1, 2, 42);
    const base_t second  = GENERATE(as<base_t>{}, 0, 1, 2, 42);
    const base_t third  = GENERATE(as<base_t>{}, 0, 1, 2, 42);
    const base_t update = 60;

    CAPTURE(first, second, third, update, first_max, second_max, third_max);

    REQUIRE(first  <= first_max);
    REQUIRE(second <= second_max);
    REQUIRE(third  <= third_max);
    REQUIRE(update <= first_max);
    REQUIRE(update <= second_max);
    REQUIRE(update <= third_max);

    CHECK(pack_t::is_valid_first(first));
    CHECK(pack_t::is_valid_second(second));
    CHECK(pack_t::is_valid_third(third));
    CHECK(pack_t::is_valid_first(update));
    CHECK(pack_t::is_valid_second(update));
    CHECK(pack_t::is_valid_third(update));

    SECTION("pack size")
    {
        CHECK(sizeof(pack_t) == sizeof(base_t));
    }

    SECTION("empty pack")
    {
        pack_t empty     = pack_t();
        pack_t empty_too = pack_t::empty();

        CHECK(empty == empty_too);
        CHECK(empty.first() == 0);
        CHECK(empty.second() == 0);
        CHECK(empty.third() == 0);
    }

    SECTION("set and get pack")
    {
        pack_t pack(first, second, third);

        CHECK(pack.first() == first);
        CHECK(pack.second() == second);
        CHECK(pack.third() == third);

        CHECK(pack.template get<0>() == first);
        CHECK(pack.template get<1>() == second);
        CHECK(pack.template get<2>() == third);

        CHECK(get<0>(pack) == first);
        CHECK(get<1>(pack) == second);
        CHECK(get<2>(pack) == third);

        SECTION("equality operator")
        {
            pack_t pack_too = pack;

            CHECK(pack == pack_too);
        }

        SECTION("update first")
        {
            pack.first(update);

            CHECK(pack.first() == update);
            CHECK(pack.second() == second);
            CHECK(pack.third() == third);
        }

        SECTION("update second")
        {
            pack.template set<1>(update);

            CHECK(pack.first() == first);
            CHECK(pack.second() == update);
            CHECK(pack.third() == third);
        }

        SECTION("update third")
        {
            pack.third(update);

            CHECK(pack.first() == first);
            CHECK(pack.second() == second);
            CHECK(pack.third() == update);
        }

        SECTION("maximum first value")
        {
            pack.first(first_max);

            CHECK(first_max != 0);
            CHECK(pack.first() == first_max);
            CHECK(pack.second() == second);
            CHECK(pack.third() == third);
        }

        SECTION("maximum second value")
        {
            pack.second(second_max);

            CHECK(second_max != 0);
            CHECK(pack.first() == first);
            CHECK(pack.second() == second_max);
            CHECK(pack.third() == third);
        }

        SECTION("maximum third value")
        {
            pack.third(third_max);

            CHECK(third_max != 0);
            CHECK(pack.first() == first);
            CHECK(pack.second() == second);
            CHECK(pack.third() == third_max);
        }

        SECTION("atomic operations")
        {
            pack_t val = pack_t(update, update, update);

            pack_t * pack_d = nullptr;
            cudaMalloc(&pack_d, sizeof(pack_t));
            REQUIRE(cudaGetLastError() == cudaSuccess);

            SECTION("atomic CAS")
            {
                pack_t compare = pack_t(pack);

                cudaMemcpy(pack_d, &pack, sizeof(pack_t), H2D);
                REQUIRE(cudaGetLastError() == cudaSuccess);

                lambda_kernel
                <<<1, 1>>>([=] DEVICEQUALIFIER
                {
                    atomicCAS(pack_d, compare, val);
                });

                cudaMemcpy(&pack, pack_d, sizeof(pack_t), D2H);
                REQUIRE(cudaGetLastError() == cudaSuccess);

                CHECK(pack == val);
            }

            SECTION("atomic exchange")
            {
                cudaMemcpy(pack_d, &pack, sizeof(pack_t), H2D);
                REQUIRE(cudaGetLastError() == cudaSuccess);

                lambda_kernel
                <<<1, 1>>>([=] DEVICEQUALIFIER
                {
                    atomicExch(pack_d, val);
                });

                cudaMemcpy(&pack, pack_d, sizeof(pack_t), D2H);
                REQUIRE(cudaGetLastError() == cudaSuccess);

                CHECK(pack == val);
            }

            cudaFree(pack_d);
            CHECK(cudaGetLastError() == cudaSuccess);
        }
    }
}

TEMPLATE_TEST_CASE_SIG(
    "PackedQuadruple with variable split",
    "[pack][quadruple][packedquadruple][variablesplit][template]",
    ((std::uint8_t FirstBits, std::uint8_t SecondBits, std::uint8_t ThirdBits, std::uint8_t FourthBits),
        FirstBits, SecondBits, ThirdBits, FourthBits),
        (8, 8, 8, 8),
        (7, 9, 9, 7),
        (9, 8, 7, 8),
        (7, 7, 7, 7),
        (16, 16, 16, 16),
        (15, 17, 13, 19),
        (8, 8, 32, 16))
{
    using pack_t = PackedQuadruple<FirstBits, SecondBits, ThirdBits, FourthBits>;
    using base_t = typename pack_t::base_type;

    REQUIRE(FirstBits + SecondBits + ThirdBits + FourthBits <= sizeof(base_t) * 8);

    const base_t first_max = (base_t{1} << pack_t::first_bits()) - base_t{1};
    const base_t second_max = (base_t{1} << pack_t::second_bits()) - base_t{1};
    const base_t third_max = (base_t{1} << pack_t::third_bits()) - base_t{1};
    const base_t fourth_max = (base_t{1} << pack_t::fourth_bits()) - base_t{1};
    const base_t first  = GENERATE(as<base_t>{}, 0, 1, 2, 42);
    const base_t second  = GENERATE(as<base_t>{}, 0, 1, 2, 42);
    const base_t third  = GENERATE(as<base_t>{}, 0, 1, 2, 42);
    const base_t fourth  = GENERATE(as<base_t>{}, 0, 1, 2, 42);
    const base_t update = 60;

    CAPTURE(first, second, third, fourth, update, first_max, second_max, third_max, fourth_max);

    REQUIRE(first  <= first_max);
    REQUIRE(second <= second_max);
    REQUIRE(third  <= third_max);
    REQUIRE(fourth <= fourth_max);
    REQUIRE(update <= first_max);
    REQUIRE(update <= second_max);
    REQUIRE(update <= third_max);
    REQUIRE(update <= fourth_max);

    CHECK(pack_t::is_valid_first(first));
    CHECK(pack_t::is_valid_second(second));
    CHECK(pack_t::is_valid_third(third));
    CHECK(pack_t::is_valid_fourth(fourth));
    CHECK(pack_t::is_valid_first(update));
    CHECK(pack_t::is_valid_second(update));
    CHECK(pack_t::is_valid_third(update));
    CHECK(pack_t::is_valid_fourth(update));

    SECTION("pack size")
    {
        CHECK(sizeof(pack_t) == sizeof(base_t));
    }

    SECTION("empty pack")
    {
        pack_t empty     = pack_t();
        pack_t empty_too = pack_t::empty();

        CHECK(empty == empty_too);
        CHECK(empty.first() == 0);
        CHECK(empty.second() == 0);
        CHECK(empty.third() == 0);
        CHECK(empty.fourth() == 0);
    }

    SECTION("set and get pack")
    {
        pack_t pack(first, second, third, fourth);

        CHECK(pack.first() == first);
        CHECK(pack.second() == second);
        CHECK(pack.third() == third);
        CHECK(pack.fourth() == fourth);

        CHECK(pack.template get<0>() == first);
        CHECK(pack.template get<1>() == second);
        CHECK(pack.template get<2>() == third);
        CHECK(pack.template get<3>() == fourth);

        CHECK(get<0>(pack) == first);
        CHECK(get<1>(pack) == second);
        CHECK(get<2>(pack) == third);
        CHECK(get<3>(pack) == fourth);

        SECTION("equality operator")
        {
            pack_t pack_too = pack;

            CHECK(pack == pack_too);
        }

        SECTION("update first")
        {
            pack.first(update);

            CHECK(pack.first() == update);
            CHECK(pack.second() == second);
            CHECK(pack.third() == third);
            CHECK(pack.fourth() == fourth);
        }

        SECTION("update second")
        {
            pack.template set<1>(update);

            CHECK(pack.first() == first);
            CHECK(pack.second() == update);
            CHECK(pack.third() == third);
            CHECK(pack.fourth() == fourth);
        }

        SECTION("update third")
        {
            pack.third(update);

            CHECK(pack.first() == first);
            CHECK(pack.second() == second);
            CHECK(pack.third() == update);
            CHECK(pack.fourth() == fourth);
        }

        SECTION("update fourth")
        {
            pack.fourth(update);

            CHECK(pack.first() == first);
            CHECK(pack.second() == second);
            CHECK(pack.third() == third);
            CHECK(pack.fourth() == update);
        }

        SECTION("maximum first value")
        {
            pack.first(first_max);

            CHECK(first_max != 0);
            CHECK(pack.first() == first_max);
            CHECK(pack.second() == second);
            CHECK(pack.third() == third);
            CHECK(pack.fourth() == fourth);
        }

        SECTION("maximum second value")
        {
            pack.second(second_max);

            CHECK(second_max != 0);
            CHECK(pack.first() == first);
            CHECK(pack.second() == second_max);
            CHECK(pack.third() == third);
            CHECK(pack.fourth() == fourth);
        }

        SECTION("maximum third value")
        {
            pack.third(third_max);

            CHECK(third_max != 0);
            CHECK(pack.first() == first);
            CHECK(pack.second() == second);
            CHECK(pack.third() == third_max);
            CHECK(pack.fourth() == fourth);
        }

        SECTION("maximum third value")
        {
            pack.fourth(fourth_max);

            CHECK(third_max != 0);
            CHECK(pack.first() == first);
            CHECK(pack.second() == second);
            CHECK(pack.third() == third);
            CHECK(pack.fourth() == fourth_max);
        }

        SECTION("atomic operations")
        {
            pack_t val = pack_t(update, update, update, update);

            pack_t * pack_d = nullptr;
            cudaMalloc(&pack_d, sizeof(pack_t));
            REQUIRE(cudaGetLastError() == cudaSuccess);

            SECTION("atomic CAS")
            {
                pack_t compare = pack_t(pack);

                cudaMemcpy(pack_d, &pack, sizeof(pack_t), H2D);
                REQUIRE(cudaGetLastError() == cudaSuccess);

                lambda_kernel
                <<<1, 1>>>([=] DEVICEQUALIFIER
                {
                    atomicCAS(pack_d, compare, val);
                });

                cudaMemcpy(&pack, pack_d, sizeof(pack_t), D2H);
                REQUIRE(cudaGetLastError() == cudaSuccess);

                CHECK(pack == val);
            }

            SECTION("atomic exchange")
            {
                cudaMemcpy(pack_d, &pack, sizeof(pack_t), H2D);
                REQUIRE(cudaGetLastError() == cudaSuccess);

                lambda_kernel
                <<<1, 1>>>([=] DEVICEQUALIFIER
                {
                    atomicExch(pack_d, val);
                });

                cudaMemcpy(&pack, pack_d, sizeof(pack_t), D2H);
                REQUIRE(cudaGetLastError() == cudaSuccess);

                CHECK(pack == val);
            }

            cudaFree(pack_d);
            CHECK(cudaGetLastError() == cudaSuccess);
        }
    }
}

/* TODO
TEST_CASE(
    "Data type reinterpretation",
    "[pack][reinterpret][template]")
{
    SECTION("pack of two std::uint8_t and one float")
    {
        using pack_t = PackedTriple<std::uint64_t, 8, 8, 32>;

        const float a = GENERATE(as<std::uint8_t>{}, 0, 1, 3, 255);
        const float b = GENERATE(as<std::uint8_t>{}, 0, 1, 3, 255);
        const float c = GENERATE(as<float>{}, 1.0, 1.22, 3.14, 9999999.9999);

        CAPTURE(a, b, c);

        pack_t pack{a, b, c};

        CHECK(pack.first<std::uint8_t>() == a);
        CHECK(pack.get<1, std::uint8_t>() == b);
        //CHECK(pack.third<float>() == c);
        //CHECK(pack.get<2, float>() == c);
    }
}
*/