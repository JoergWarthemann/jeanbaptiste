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
#include "../../JeanBaptiste/include/windowing/WelchWindow.h"
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
    Utilities::WindowingAnalysis<double> analysis_;
    jt::Real2Complex<std::complex<double>> real2ComplexConverter_;
    jt::Complex2Real<std::complex<double>> complex2RealConverter_;
    bool initialized_;

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
        BOOST_TEST((initialized_ = analysis_.initialize("../../test cases/WinBartlettTest.xml", "win.in", workingSet_, "win.out", expectedOut_)), "Loading test data failed.");
        if (!initialized_)
            return;

        BOOST_TEST_MESSAGE("Checking Bartlett window samples.");

        jw::BartlettWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> bartlett;

        auto complexData = real2ComplexConverter_(workingSet_);
        bartlett(&complexData[0]);
        auto realData = complex2RealConverter_(complexData);

        analysis_.checkOutput(realData, expectedOut_);
    }

    BOOST_AUTO_TEST_CASE(blackman_harris)
    {
        BOOST_TEST((initialized_ = analysis_.initialize("../../test cases/WinBlackmanHarrisTest.xml", "win.in", workingSet_, "win.out", expectedOut_)), "Loading test data failed.");
        if (!initialized_)
            return;

        BOOST_TEST_MESSAGE("Checking Blackman Harris window samples.");

        jw::BlackmanHarrisWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> blackmanHarris;

        auto complexData = real2ComplexConverter_(workingSet_);
        blackmanHarris(&complexData[0]);
        auto realData = complex2RealConverter_(complexData);

        analysis_.checkOutput(realData, expectedOut_);
    }

    BOOST_AUTO_TEST_CASE(blackman)
    {
        BOOST_TEST((initialized_ = analysis_.initialize("../../test cases/WinBlackmanTest.xml", "win.in", workingSet_, "win.out", expectedOut_)), "Loading test data failed.");
        if (!initialized_)
            return;

        BOOST_TEST_MESSAGE("Checking Blackman window samples.");

        jw::BlackmanWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> blackman;

        auto complexData = real2ComplexConverter_(workingSet_);
        blackman(&complexData[0]);
        auto realData = complex2RealConverter_(complexData);

        analysis_.checkOutput(realData, expectedOut_);
    }

    BOOST_AUTO_TEST_CASE(cosine)
    {
        BOOST_TEST((initialized_ = analysis_.initialize("../../test cases/WinCosineTest.xml", "win.in", workingSet_, "win.out", expectedOut_)), "Loading test data failed.");
        if (!initialized_)
            return;

        BOOST_TEST_MESSAGE("Checking Cosine window samples.");

        jw::CosineWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> cosineWin;

        auto complexData = real2ComplexConverter_(workingSet_);
        cosineWin(&complexData[0]);
        auto realData = complex2RealConverter_(complexData);

        analysis_.checkOutput(realData, expectedOut_);
    }

    BOOST_AUTO_TEST_CASE(flat_top)
    {
        BOOST_TEST((initialized_ = analysis_.initialize("../../test cases/WinFlatTopTest.xml", "win.in", workingSet_, "win.out", expectedOut_)), "Loading test data failed.");
        if (!initialized_)
            return;

        BOOST_TEST_MESSAGE("Checking Flat Top window samples.");

        jw::FlatTopWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> flatTopWin;

        auto complexData = real2ComplexConverter_(workingSet_);
        flatTopWin(&complexData[0]);
        auto realData = complex2RealConverter_(complexData);

        analysis_.checkOutput(realData, expectedOut_);
    }

    BOOST_AUTO_TEST_CASE(hamming)
    {
        BOOST_TEST((initialized_ = analysis_.initialize("../../test cases/WinHammingTest.xml", "win.in", workingSet_, "win.out", expectedOut_)), "Loading test data failed.");
        if (!initialized_)
            return;

        BOOST_TEST_MESSAGE("Checking Hamming window samples.");

        jw::HammingWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> hammingWin;

        auto complexData = real2ComplexConverter_(workingSet_);
        hammingWin(&complexData[0]);
        auto realData = complex2RealConverter_(complexData);

        analysis_.checkOutput(realData, expectedOut_);
    }

    BOOST_AUTO_TEST_CASE(von_hann)
    {
        BOOST_TEST((initialized_ = analysis_.initialize("../../test cases/WinvonHannTest.xml", "win.in", workingSet_, "win.out", expectedOut_)), "Loading test data failed.");
        if (!initialized_)
            return;

        BOOST_TEST_MESSAGE("Checking von Hann window samples.");

        jw::VonHannWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> vonHannWin;

        auto complexData = real2ComplexConverter_(workingSet_);
        vonHannWin(&complexData[0]);
        auto realData = complex2RealConverter_(complexData);

        analysis_.checkOutput(realData, expectedOut_);
    }

    BOOST_AUTO_TEST_CASE(welch)
    {
        BOOST_TEST((initialized_ = analysis_.initialize("../../test cases/WinWelchTest.xml", "win.in", workingSet_, "win.out", expectedOut_)), "Loading test data failed.");
        if (!initialized_)
            return;

        BOOST_TEST_MESSAGE("Checking Welch window samples.");

        jw::WelchWindow<std::integral_constant<int, kSampleCnt_>, std::complex<double>> welchWin;

        auto complexData = real2ComplexConverter_(workingSet_);
        welchWin(&complexData[0]);
        auto realData = complex2RealConverter_(complexData);

        analysis_.checkOutput(realData, expectedOut_);
    }

BOOST_AUTO_TEST_SUITE_END()