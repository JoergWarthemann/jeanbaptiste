#include <boost/test/unit_test.hpp>
#include <complex>
#include "../include/WindowingAnalysis.h"
#include "../../JeanBaptiste/include/windowing/BartlettWindow.h"
#include <string>

namespace ut = boost::unit_test;
namespace jb = jeanbaptiste;
namespace jw = jeanbaptiste::windowing;

class WindowCalculationFixture
{
protected:
    std::vector<double> workingSet_;
    std::vector<double> expectedOut_;
    static const unsigned kSampleCnt_ = 128;
    Utilities::WindowingAnalysis<double> _analysis;

public:
    WindowCalculationFixture()
    {
        workingSet_.clear();
        expectedOut_.clear();
    }

    ~WindowCalculationFixture()
    {}
};


BOOST_FIXTURE_TEST_SUITE(WindowCalculationTestSuite, WindowCalculationFixture)

    BOOST_AUTO_TEST_CASE(bartlett)
    {
        BOOST_TEST_MESSAGE("Checking bartlett window samples.");
        _analysis.initialize("././test cases/WinBartlettTest.xml", "win.in", workingSet_, "win.out", expectedOut_);

        jw::BartlettWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> bartlett;

        int i = 0;
    }

    BOOST_AUTO_TEST_CASE(wait)
    {
        std::cin.get();
    }

BOOST_AUTO_TEST_SUITE_END()