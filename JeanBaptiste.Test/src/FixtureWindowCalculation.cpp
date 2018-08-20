#include <boost/test/unit_test.hpp>
#include <complex>
#include "../include/WindowingAnalysis.h"
#include "../../JeanBaptiste/include/tools/RealComplexConversion.h"
#include "../../JeanBaptiste/include/windowing/BartlettWindow.h"
#include "../../JeanBaptiste/include/windowing/BlackmanHarrisWindow.h"
#include <string>

namespace ut = boost::unit_test;
namespace jb = jeanbaptiste;
namespace jt = jeanbaptiste::tools;
namespace jw = jeanbaptiste::windowing;

class WindowCalculationFixture
{
protected:
    std::vector<double> workingSet_;
    std::vector<double> expectedOut_;
    static const unsigned kSampleCnt_ = 128;
    Utilities::WindowingAnalysis<double> _analysis;
    jt::Real2Complex<std::complex<double>> _real2ComplexConverter;
    jt::Complex2Real<std::complex<double>> _complex2RealConverter;

public:
    WindowCalculationFixture()
    {
        workingSet_.clear();
        expectedOut_.clear();
    }

    ~WindowCalculationFixture()
    {}

    //void check()
    //{
    //    std::vector<std::complex<double>> 
    //    _converter()(expectedOut_, );
    //}
};


BOOST_FIXTURE_TEST_SUITE(WindowCalculationTestSuite, WindowCalculationFixture)

    BOOST_AUTO_TEST_CASE(bartlett)
    {
        BOOST_TEST_MESSAGE("Checking bartlett window samples.");
        _analysis.initialize("././test cases/WinBartlettTest.xml", "win.in", workingSet_, "win.out", expectedOut_);

        jw::BartlettWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> bartlett;

        auto complexData = _real2ComplexConverter(workingSet_);
        bartlett(&complexData[0]);
        auto realData = _complex2RealConverter(complexData);

        _analysis.checkOutput(realData, expectedOut_);
    }

    BOOST_AUTO_TEST_CASE(blackman_harris)
    {
        BOOST_TEST_MESSAGE("Checking blackman harris window samples.");
        _analysis.initialize("././test cases/WinBlackmanHarrisTest.xml", "win.in", workingSet_, "win.out", expectedOut_);

        jw::BlackmanHarrisWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> bartlett;

        //auto complexData = _real2ComplexConverter(workingSet_);
        //bartlett(&complexData[0]);
        //auto realData = _complex2RealConverter(complexData);
//
        //_analysis.checkOutput(realData, expectedOut_);
    }

    BOOST_AUTO_TEST_CASE(wait)
    {
        std::cin.get();
    }

BOOST_AUTO_TEST_SUITE_END()