#include <boost/test/unit_test.hpp>
#include <complex>
#include "../../JeanBaptiste/include/AlgorithmFactory.h"
#include "../include/AlgorithmFixture.h"
#include <string>

namespace ut = boost::unit_test;
namespace jb = jeanbaptiste;
namespace jbo = jeanbaptiste::options;

class Radix4Fixture
    : public AlgorithmFixture
{
public:
    Radix4Fixture()
        : AlgorithmFixture()
    {
        BOOST_TEST_MESSAGE("Setup fixture: square pulse of 64 samples.");

		algorithmResult_.initialize("././test cases/square pulse (n=64).xml", "fft.in", workingSet_, expectedOutIFFT_, "fft.out", expectedOutFFT_);
    }

    ~Radix4Fixture()
    {}
};


BOOST_FIXTURE_TEST_SUITE(Radix4TestSuite, Radix4Fixture)
 
    BOOST_AUTO_TEST_CASE(fft_radix_4_dif)
    {
        BOOST_TEST_MESSAGE("Running radix 4 DIF FFT and IFFT.");

        // Create Radix-4 DIF FFT algorithms for sample counts 4 ... 1024.
        jb::AlgorithmFactory<1, 5, jbo::Radix_4, jbo::Decimation_In_Frequency, jbo::Direction_Forward, jbo::Window_None,
        std::complex<double>> fftFactory;

        // Create Radix-4 DIF IFFT algorithms for sample counts 4 ... 1024.
        jb::AlgorithmFactory<1, 5, jbo::Radix_4, jbo::Decimation_In_Frequency, jbo::Direction_Backward, jbo::Window_None,
        std::complex<double>> ifftFactory;

        runAlgorithms(fftFactory.getAlgorithm(3), ifftFactory.getAlgorithm(3));
    }

    BOOST_AUTO_TEST_CASE(fft_radix_4_dit)
    {
        BOOST_TEST_MESSAGE("Running radix 4 DIT FFT and IFFT.");

        // Create Radix-4 DIT FFT algorithms for sample counts 4 ... 1024.
        jb::AlgorithmFactory<1, 5, jbo::Radix_4, jbo::Decimation_In_Time, jbo::Direction_Forward, jbo::Window_None,
        std::complex<double>> fftFactory;

        // Create Radix-4 DIT IFFT algorithms for sample counts 4 ... 1024.
        jb::AlgorithmFactory<1, 5, jbo::Radix_4, jbo::Decimation_In_Time, jbo::Direction_Backward, jbo::Window_None,
        std::complex<double>> ifftFactory;

        runAlgorithms(fftFactory.getAlgorithm(3), ifftFactory.getAlgorithm(3));
    }

    BOOST_AUTO_TEST_CASE(wait)
    {
        std::cin.get();
    }

BOOST_AUTO_TEST_SUITE_END()