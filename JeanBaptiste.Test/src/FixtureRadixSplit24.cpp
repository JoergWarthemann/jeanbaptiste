#include <boost/test/unit_test.hpp>
#include <complex>
#include "../../JeanBaptiste/include/AlgorithmFactory.h"
#include "../include/AlgorithmFixture.h"
#include <string>

namespace ut = boost::unit_test;
namespace jb = jeanbaptiste;
namespace jbo = jeanbaptiste::options;

class RadixSplit24Fixture
    : public AlgorithmFixture
{
protected:
    bool initialized_;

public:
    RadixSplit24Fixture()
        : AlgorithmFixture()
    {
        BOOST_TEST_MESSAGE("Setup fixture: square pulse of 64 samples.");
        BOOST_TEST((initialized_ = algorithmResult_.initialize("../../test cases/square pulse (n=64).xml", "fft.in", workingSet_, expectedOutIFFT_, "fft.out", expectedOutFFT_)), "Loading test data failed.");
    }

    ~RadixSplit24Fixture()
    {}
};


BOOST_FIXTURE_TEST_SUITE(RadixSplit24TestSuite, RadixSplit24Fixture)

    BOOST_AUTO_TEST_CASE(fft_radix_split_2_4_dif)
    {
        if (!initialized_)
            return;

        BOOST_TEST_MESSAGE("Running split radix 2-4 DIF FFT and IFFT.");

        // Create Split-Radix-2-4 DIF FFT algorithms for sample counts 2 ... 256.
        jb::AlgorithmFactory<1, 8, jbo::Radix_Split_2_4, jbo::Decimation_In_Frequency, jbo::Direction_Forward, jbo::Window_None,
            jbo::Normalization_Square_Root, std::complex<double>> fftFactory;

        // Create Split-Radix-2-4 DIF IFFT algorithms for sample counts 2 ... 256.
        jb::AlgorithmFactory<1, 8, jbo::Radix_Split_2_4, jbo::Decimation_In_Frequency, jbo::Direction_Backward, jbo::Window_None,
            jbo::Normalization_Square_Root, std::complex<double>> ifftFactory;

        runAlgorithms(fftFactory.getAlgorithm(6), ifftFactory.getAlgorithm(6));
    }

    BOOST_AUTO_TEST_CASE(fft_radix_split_2_4_dit)
    {
        if (!initialized_)
            return;

        BOOST_TEST_MESSAGE("Running split radix 2-4 DIT FFT and IFFT.");

        // Create Split-Radix-2-4 DIT FFT algorithms for sample counts 2 ... 256.
        jb::AlgorithmFactory<1, 8, jbo::Radix_Split_2_4, jbo::Decimation_In_Time, jbo::Direction_Forward, jbo::Window_None,
            jbo::Normalization_Square_Root, std::complex<double>> fftFactory;

        // Create Split-Radix-2-4 DIT IFFT algorithms for sample counts 2 ... 256.
        jb::AlgorithmFactory<1, 8, jbo::Radix_Split_2_4, jbo::Decimation_In_Time, jbo::Direction_Backward, jbo::Window_None,
            jbo::Normalization_Square_Root, std::complex<double>> ifftFactory;

        runAlgorithms(fftFactory.getAlgorithm(6), ifftFactory.getAlgorithm(6));
    }

BOOST_AUTO_TEST_SUITE_END()