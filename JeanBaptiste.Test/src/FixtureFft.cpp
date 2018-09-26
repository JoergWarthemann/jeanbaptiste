#include <boost/test/unit_test.hpp>
#include <complex>
#include "../../JeanBaptiste/include/AlgorithmFactory.h"
#include "../include/AlgorithmFixture.h"
#include <string>

namespace ut = boost::unit_test;
namespace jb = jeanbaptiste;
namespace jbo = jeanbaptiste::options;

class FFTFixture
    : public AlgorithmFixture
{
public:
    FFTFixture()
        : AlgorithmFixture()
    {
        BOOST_TEST_MESSAGE("Setup fixture: cos(2*pi*n/10) of 128 samples.");

		//algorithmResult_.initialize("././test cases/cosine (n=128) + Bartlett.xml", "fft.in", workingSet_, expectedOutIFFT_, "fft.out", expectedOutFFT_);
    }

    ~FFTFixture()
    {}
};


BOOST_FIXTURE_TEST_SUITE(FftTestSuite, FFTFixture)

    BOOST_AUTO_TEST_CASE(fft_radix_2_dif_bartlett)
    {
        BOOST_TEST_MESSAGE("Running radix 2 DIF FFT and IFFT.");

        algorithmResult_.initialize("././test cases/cosine (n=128) + Bartlett.xml", "fft.in", workingSet_, expectedOutIFFT_, "fft.out", expectedOutFFT_);

        // Create Radix-2 DIF FFT algorithms for sample counts 2 ... 256.
        jb::AlgorithmFactory<1, 8, jbo::Radix_2, jbo::Decimation_In_Frequency, jbo::Direction_Forward, jbo::Window_Bartlett,
            jbo::Normalization_Square_Root, std::complex<double>> fftFactory;

        runAlgorithm(fftFactory.getAlgorithm(7));

        // Create Radix-2 DIF IFFT algorithms for sample counts 2 ... 256.
        jb::AlgorithmFactory<1, 8, jbo::Radix_2, jbo::Decimation_In_Frequency, jbo::Direction_Backward, jbo::Window_None,
            jbo::Normalization_No, std::complex<double>> ifftFactory;

        //runAlgorithms(fftFactory.getAlgorithm(7), ifftFactory.getAlgorithm(7));
    }

    BOOST_AUTO_TEST_CASE(wait)
    {
        std::cin.get();
    }

BOOST_AUTO_TEST_SUITE_END()