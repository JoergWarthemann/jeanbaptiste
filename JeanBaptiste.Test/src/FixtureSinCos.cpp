#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <complex>
#include <iostream>
#include "../../JeanBaptiste/include/basic/SineCosine.h"
#include <memory>
#include <string>

namespace constants = boost::math::constants;
namespace jb = jeanbaptiste;
namespace jbb = jeanbaptiste::basic;
namespace ut = boost::unit_test;

class SinCosFixture
{
public:
    SinCosFixture()
    {}

    ~SinCosFixture()
    {}
};


BOOST_FIXTURE_TEST_SUITE(SinCosTestSuite, SinCosFixture)

    BOOST_AUTO_TEST_CASE(float_sine_from_minus_2pi_to_plus_2_pi)
    {
        BOOST_TEST_MESSAGE("Running sine computation in float range [-2pi ... 2pi].");

        for (float x = -constants::two_pi<float>(); x <=constants::two_pi<float>(); x += 0.1)
            BOOST_TEST(std::fabs(std::sin(x) - jbb::sine<float>(x)) < 0.00001);
    }

    BOOST_AUTO_TEST_CASE(double_sine_from_minus_2pi_to_plus_2_pi)
    {
        BOOST_TEST_MESSAGE("Running sine computation in double range [-2pi ... 2pi].");

        for (double x = -constants::two_pi<double>(); x <=constants::two_pi<double>(); x += 0.1)
            BOOST_TEST(std::fabs(std::sin(x) - jbb::sine<double>(x)) < 0.000000000001);
    }

    BOOST_AUTO_TEST_CASE(float_cosine_from_minus_2pi_to_plus_2_pi)
    {
        BOOST_TEST_MESSAGE("Running cosine computation in float range [-2pi ... 2pi].");

        for (float x = -constants::two_pi<float>(); x <=constants::two_pi<float>(); x += 0.1)
            BOOST_TEST(std::fabs(std::cos(x) - jbb::cosine<float>(x)) < 0.0001);
    }

    BOOST_AUTO_TEST_CASE(double_cosine_from_minus_2pi_to_plus_2_pi)
    {
        BOOST_TEST_MESSAGE("Running cosine computation in double range [-2pi ... 2pi].");

        for (double x = -constants::two_pi<double>(); x <=constants::two_pi<double>(); x += 0.1)
            BOOST_TEST(std::fabs(std::cos(x) - jbb::cosine<double>(x)) < 0.00000000001);
    }

    BOOST_AUTO_TEST_CASE(wait)
    {
        std::cin.get();
    }

BOOST_AUTO_TEST_SUITE_END()