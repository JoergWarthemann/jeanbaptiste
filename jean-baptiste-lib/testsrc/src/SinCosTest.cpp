#include <algorithm>
#include <cmath>
#include <complex>

#include <gtest/gtest.h>

#include "tools/SineCosine.hpp"

#include <iostream>

namespace JB::Testing {
class SinCosTest
    : public ::testing::Test {
public:
    SinCosTest(void) = default;
    ~SinCosTest(void) = default;
};

TEST_F(SinCosTest, SuccessfullyCalculatesSineFromMinus2PiToPlus2Pi)
{
    constexpr auto floatPrecision{0.00001};
    constexpr auto doublePrecision{0.000000000001};

    // Running sine computation in float range [-2pi ... 2pi].
    for (float x = -constants::two_pi<float>(); x <= constants::two_pi<float>(); x += 0.1) {
        EXPECT_NEAR(std::sin(x), jeanbaptiste::tools::sine<float>(x), floatPrecision);
    }

    // Running sine computation in double range [-2pi ... 2pi].
    for (double x = -constants::two_pi<double>(); x <= constants::two_pi<double>(); x += 0.1) {
        EXPECT_NEAR(std::sin(x), jeanbaptiste::tools::sine<double>(x), doublePrecision);
    }
}

TEST_F(SinCosTest, SuccessfullyCalculatesCosineFromMinus2PiToPlus2Pi)
{
    constexpr auto floatPrecision{0.0001};
    constexpr auto doublePrecision{0.00000000001};

    // Running cosine computation in float range [-2pi ... 2pi].
    for (float x = -constants::two_pi<float>(); x <= constants::two_pi<float>(); x += 0.1) {
        EXPECT_NEAR(std::cos(x), jeanbaptiste::tools::cosine<float>(x), floatPrecision);
    }

    // Running cosine computation in double range [-2pi ... 2pi].
    for (double x = -constants::two_pi<double>(); x <= constants::two_pi<double>(); x += 0.1) {
        EXPECT_NEAR(std::cos(x), jeanbaptiste::tools::cosine<double>(x), doublePrecision);
    }
}

TEST_F(SinCosTest, SuccessfullyReducesAnArgumentsRange)
{
    constexpr auto precision{0.000000001};

    // Check some specific values.
    ASSERT_NEAR(3.716814693, jeanbaptiste::tools::internal::reduceRange<double>(10), precision);
    ASSERT_NEAR(6.01770285, jeanbaptiste::tools::internal::reduceRange<double>(50), precision);
    ASSERT_NEAR(5.752220392, jeanbaptiste::tools::internal::reduceRange<double>(100), precision);
    ASSERT_NEAR(3.628360733, jeanbaptiste::tools::internal::reduceRange<double>(500), precision);
    ASSERT_NEAR(0.973536158, jeanbaptiste::tools::internal::reduceRange<double>(1000), precision);
}

TEST_F(SinCosTest, FailsWhenTryingToReduceAnOutOfScopeArgumentsRange)
{
    constexpr auto precision{0.000000001};
    constexpr auto expected{1.967870092};
    constexpr auto difference = std::abs(expected - jeanbaptiste::tools::internal::reduceRange<double>(1001));

    ASSERT_GT(jeanbaptiste::tools::internal::reduceRange<double>(difference), precision);
}

} // namespace JB::Testing