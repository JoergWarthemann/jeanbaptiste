#include <algorithm>
#include <array>
#include <complex>
#include <ranges>
#include <span>

#include <gtest/gtest.h>

#include "core/Radix2.hpp"
#include "AlgorithmFixture.hpp"

namespace jeanbaptiste::testing {
class Radix2Test
    : public ::testing::Test
    , AlgorithmFixture {
public:
    Radix2Test()
    {
        constexpr std::string_view testName{"square pulse (n=64)"};
    }

    ~Radix2Test() override = default;

protected:
    bool mInitialized{};
};

TEST_F(Radix2Test, SuccessfullyCalculatesRadix2Dit)
{
    GTEST_SKIP() << "TODO: Add Radix2 test for DIT FFT";
}

TEST_F(Radix2Test, SuccessfullyCalculatesRadix2Dif)
{
    GTEST_SKIP() << "TODO: Add Radix2 test for DIF FFT";
}

} // namespace jeanbaptiste::testing