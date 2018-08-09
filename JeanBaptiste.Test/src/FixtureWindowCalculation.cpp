#include <boost/test/unit_test.hpp>
#include <complex>
#include "../include/AlgorithmResultAnalysis.h"
#include "../../JeanBaptiste/include/windowing/BartlettWindow.h"
#include <string>

namespace ut = boost::unit_test;
namespace jb = jeanbaptiste;
namespace jbo = jeanbaptiste::options;

class WindowCalculationFixture
{
    std::vector<std::complex<double>> workingSet_;
	std::vector<std::complex<double>> expectedOut_;
    Utilities::AlgorithmResultAnalysis<double> algorithmResult_;

public:
    WindowCalculationFixture()
    {
    }

    ~WindowCalculationFixture()
    {}
};


BOOST_FIXTURE_TEST_SUITE(WindowCalculationTestSuite, WindowCalculationFixture)

    BOOST_AUTO_TEST_CASE(bartlett)
    {
        BOOST_TEST_MESSAGE("Checking bartlett window samples.");
        algorithmResult_.initialize("././test cases/WinBartlettTest.xml", "win.in", workingSet_, expectedOutIFF_, "win.out", expectedOutFFT_);
    }

    BOOST_AUTO_TEST_CASE(wait)
    {
        std::cin.get();
    }

BOOST_AUTO_TEST_SUITE_END()