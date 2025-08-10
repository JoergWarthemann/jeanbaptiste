#include <algorithm>
#include <array>
#include <complex>
#include <ranges>
#include <span>

#include <gtest/gtest.h>

#include "core/Radix2.hpp"

namespace JB::Testing {
class Radix2Test
    : public ::testing::Test {
public:
    Radix2Test(void) = default;
    ~Radix2Test(void) = default;
};

TEST_F(Radix2Test, SuccessfullyCalculatesRadix2Dit)
{
    GTEST_SKIP() << "TODO: Add Radix2 test for DIT FFT";
}

TEST_F(Radix2Test, SuccessfullyCalculatesRadix2Dif)
{
    GTEST_SKIP() << "TODO: Add Radix2 test for DIF FFT";
}

} // namespace JB::Testing