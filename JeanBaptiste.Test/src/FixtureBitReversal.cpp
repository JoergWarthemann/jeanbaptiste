#include <boost/test/unit_test.hpp>
#include <complex>
#include "../../JeanBaptiste/include/basic/BitReversalIndexSwapping.h"
#include <string>

namespace ut = boost::unit_test;
namespace jb = jeanbaptiste;
namespace jbb = jeanbaptiste::basic;

class BitReversalFixture
{
public:
    BitReversalFixture()
    {}

    ~BitReversalFixture()
    {}
};


BOOST_FIXTURE_TEST_SUITE(BitReversalTestSuite, BitReversalFixture)

    BOOST_AUTO_TEST_CASE(reverse_on_4_bit_base)
    {
        BOOST_TEST_MESSAGE("Running bit reversal based on 4 bit indices.");
        jbb::BitReversalIndexSwapping<typename decltype(std::integral_constant<int, 256>{})::type, std::complex<double>> reversal_;

    }

    BOOST_AUTO_TEST_CASE(wait)
    {
        std::cin.get();
    }

BOOST_AUTO_TEST_SUITE_END()