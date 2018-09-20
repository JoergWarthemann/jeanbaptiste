#include <boost/test/unit_test.hpp>
#include <complex>
#include "../include/WindowingAnalysis.h"
#include "../../JeanBaptiste/include/tools/RealComplexConversion.h"
#include "../../JeanBaptiste/include/windowing/BartlettWindow.h"
#include "../../JeanBaptiste/include/windowing/BlackmanWindow.h"
#include "../../JeanBaptiste/include/windowing/BlackmanHarrisWindow.h"
#include "../../JeanBaptiste/include/windowing/CosineWindow.h"
#include "../../JeanBaptiste/include/windowing/FlatTopWindow.h"
#include "../../JeanBaptiste/include/windowing/HammingWindow.h"
#include "../../JeanBaptiste/include/windowing/VonHannWindow.h"
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
        BOOST_TEST_MESSAGE("Checking Bartlett window samples.");
        _analysis.initialize("././test cases/WinBartlettTest.xml", "win.in", workingSet_, "win.out", expectedOut_);

        jw::BartlettWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> bartlett;

        auto complexData = _real2ComplexConverter(workingSet_);
        bartlett(&complexData[0]);
        auto realData = _complex2RealConverter(complexData);

        _analysis.checkOutput(realData, expectedOut_);
    }

    BOOST_AUTO_TEST_CASE(blackman_harris)
    {
        BOOST_TEST_MESSAGE("Checking Blackman Harris window samples.");
        _analysis.initialize("././test cases/WinBlackmanHarrisTest.xml", "win.in", workingSet_, "win.out", expectedOut_);

        jw::BlackmanHarrisWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> blackmanHarris;

        auto complexData = _real2ComplexConverter(workingSet_);
        blackmanHarris(&complexData[0]);
        auto realData = _complex2RealConverter(complexData);

        _analysis.checkOutput(realData, expectedOut_);
    }

    BOOST_AUTO_TEST_CASE(blackman)
    {
        BOOST_TEST_MESSAGE("Checking Blackman window samples.");
        _analysis.initialize("././test cases/WinBlackmanTest.xml", "win.in", workingSet_, "win.out", expectedOut_);

        jw::BlackmanWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> blackman;

        auto complexData = _real2ComplexConverter(workingSet_);
        blackman(&complexData[0]);
        auto realData = _complex2RealConverter(complexData);

        _analysis.checkOutput(realData, expectedOut_);
    }

    BOOST_AUTO_TEST_CASE(cosine)
    {
        BOOST_TEST_MESSAGE("Checking Cosine window samples.");
        _analysis.initialize("././test cases/WinCosineTest.xml", "win.in", workingSet_, "win.out", expectedOut_);

        jw::CosineWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> cosineWin;

        auto complexData = _real2ComplexConverter(workingSet_);
        cosineWin(&complexData[0]);
        auto realData = _complex2RealConverter(complexData);

        _analysis.checkOutput(realData, expectedOut_);
    }

    BOOST_AUTO_TEST_CASE(flat_top)
    {
        BOOST_TEST_MESSAGE("Checking Flat Top window samples.");
        _analysis.initialize("././test cases/WinFlatTopTest.xml", "win.in", workingSet_, "win.out", expectedOut_);

        jw::FlatTopWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> flatTopWin;

        auto complexData = _real2ComplexConverter(workingSet_);
        flatTopWin(&complexData[0]);
        auto realData = _complex2RealConverter(complexData);

        _analysis.checkOutput(realData, expectedOut_);
    }

    BOOST_AUTO_TEST_CASE(hamming)
    {
        BOOST_TEST_MESSAGE("Checking Hamming window samples.");
        _analysis.initialize("././test cases/WinHammingTest.xml", "win.in", workingSet_, "win.out", expectedOut_);

        jw::HammingWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> hammingWin;

        auto complexData = _real2ComplexConverter(workingSet_);
        hammingWin(&complexData[0]);
        auto realData = _complex2RealConverter(complexData);

        _analysis.checkOutput(realData, expectedOut_);
    }

    BOOST_AUTO_TEST_CASE(von_hann)
    {
        BOOST_TEST_MESSAGE("Checking von Hann window samples.");
        _analysis.initialize("././test cases/WinvonHannTest.xml", "win.in", workingSet_, "win.out", expectedOut_);

        jw::VonHannWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> vonHannWin;

        auto complexData = _real2ComplexConverter(workingSet_);
        vonHannWin(&complexData[0]);
        auto realData = _complex2RealConverter(complexData);

        _analysis.checkOutput(realData, expectedOut_);
    }

    BOOST_AUTO_TEST_CASE(wait)
    {
        std::cin.get();
    }

BOOST_AUTO_TEST_SUITE_END()