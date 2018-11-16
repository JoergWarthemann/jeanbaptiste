#include <boost/test/unit_test.hpp>
#include <complex>
#include "../../JeanBaptiste/include/AlgorithmFactory.h"
#include "../include/AlgorithmFixture.h"
#include <string>

namespace ut = boost::unit_test;
namespace jb = jeanbaptiste;
namespace jbo = jeanbaptiste::options;

class Radix2Fixture
    : public AlgorithmFixture
{
protected:
    bool initialized_;

public:
    Radix2Fixture()
        : AlgorithmFixture(),
          initialized_(false)
    {
        BOOST_TEST_MESSAGE("Setup fixture: square pulse of 64 samples.");
        BOOST_TEST((initialized_ = algorithmResult_.initialize("./test cases/square pulse (n=64).xml", "fft.in", workingSet_, expectedOutIFFT_, "fft.out", expectedOutFFT_)), "Error loading test data.");
    }

    ~Radix2Fixture()
    {}
};


BOOST_FIXTURE_TEST_SUITE(Radix2TestSuite, Radix2Fixture)

    BOOST_AUTO_TEST_CASE(fft_radix_2_dif)
    {
        if (!initialized_)
            return;

        BOOST_TEST_MESSAGE("Running radix 2 DIF FFT and IFFT.");

        // Create Radix-2 DIF FFT algorithms for sample counts 2 ... 256.
        jb::AlgorithmFactory<1, 8, jbo::Radix_2, jbo::Decimation_In_Frequency, jbo::Direction_Forward, jbo::Window_None,
            jbo::Normalization_Square_Root, std::complex<double>> fftFactory;

        // Create Radix-2 DIF IFFT algorithms for sample counts 2 ... 256.
        jb::AlgorithmFactory<1, 8, jbo::Radix_2, jbo::Decimation_In_Frequency, jbo::Direction_Backward, jbo::Window_None,
            jbo::Normalization_Square_Root, std::complex<double>> ifftFactory;

        runAlgorithms(fftFactory.getAlgorithm(6), ifftFactory.getAlgorithm(6));
    }

    BOOST_AUTO_TEST_CASE(fft_radix_2_dit)
    {
        if (!initialized_)
            return;

        BOOST_TEST_MESSAGE("Running radix 2 DIT FFT and IFFT.");

        // Create Radix-2 DIT FFT algorithms for sample counts 2 ... 256.
        jb::AlgorithmFactory<1, 8, jbo::Radix_2, jbo::Decimation_In_Time, jbo::Direction_Forward, jbo::Window_None,
            jbo::Normalization_Square_Root, std::complex<double>> fftFactory;

        // Create Radix-2 DIT IFFT algorithms for sample counts 2 ... 256.
        jb::AlgorithmFactory<1, 8, jbo::Radix_2, jbo::Decimation_In_Time, jbo::Direction_Backward, jbo::Window_None,
            jbo::Normalization_Square_Root, std::complex<double>> ifftFactory;

        runAlgorithms(fftFactory.getAlgorithm(6), ifftFactory.getAlgorithm(6));
    }

    BOOST_AUTO_TEST_CASE(wait)
    {
        std::cin.get();
    }

BOOST_AUTO_TEST_SUITE_END()