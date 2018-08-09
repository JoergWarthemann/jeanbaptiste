#include <boost/hana.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <complex>
#include "../include/AlgorithmFactory.h"
#include <iostream>
#include <memory>

#include "../include/basic/SineCosine.h"
#include "../include/basic/Abs.h"

#include "../include/windowing/BartlettWindow.h"

namespace constants = boost::math::constants;
namespace hana = boost::hana;
namespace jb = jeanbaptiste;
namespace jbo = jeanbaptiste::options;
namespace jw = jeanbaptiste::windowing;

int main()
{
    std::array<std::complex<double>, 16> complexData =
    {
        std::complex<double>(-1.0f, 0.0f),
		std::complex<double>(-1.0f, 0.0f),
		std::complex<double>(-1.0f, 0.0f),
		std::complex<double>(-1.0f, 0.0f),
		std::complex<double>( 1.0f, 0.0f),
		std::complex<double>( 1.0f, 0.0f),
		std::complex<double>( 1.0f, 0.0f),
		std::complex<double>( 1.0f, 0.0f),
		std::complex<double>(-1.0f, 0.0f),
		std::complex<double>(-1.0f, 0.0f),
		std::complex<double>(-1.0f, 0.0f),
		std::complex<double>(-1.0f, 0.0f),
		std::complex<double>( 1.0f, 0.0f),
		std::complex<double>( 1.0f, 0.0f),
		std::complex<double>( 1.0f, 0.0f),
		std::complex<double>( 1.0f, 0.0f)
    };

    //std::cout << "FFT:\n";
    //jb::AlgorithmFactory<1, 6, jbo::Radix_2, jbo::Decimation_In_Frequency, jbo::Direction_Forward, 
	////jb::AlgorithmFactory<1, 6, jbo::Radix_4, jbo::Decimation_In_Time, jbo::Direction_Forward,
    ////jb::AlgorithmFactory<1, 6, jbo::Radix_Split_2_4, jbo::Decimation_In_Time, jbo::Direction_Forward,  
    //    jbo::Window_None, std::complex<double>> algorithmFactory;
    //auto algorithm = algorithmFactory.getAlgorithm(4);
	////auto algorithm = algorithmFactory.getAlgorithm(2);
    ////auto algorithm = algorithmFactory.getAlgorithm(4);
    //algorithm->operator()(&complexData[0]);
//
    //for (auto value : complexData)
    //    std::cout << std::setw(10) << std::setprecision(5) << value.real() << "\t" << value.imag() << "I\n";
//
    //std::cout << "IFFT:\n";
    //jb::AlgorithmFactory<1, 6, jbo::Radix_2, jbo::Decimation_In_Frequency, jbo::Direction_Backward, 
    ////jb::AlgorithmFactory<1, 6, jbo::Radix_4, jbo::Decimation_In_Time, jbo::Direction_Backward,
    ////jb::AlgorithmFactory<1, 6, jbo::Radix_Split_2_4, jbo::Decimation_In_Time, jbo::Direction_Backward,
    //    jbo::Window_None, std::complex<double>> inverseAlgorithmFactory;
    //auto inverseAlgorithm = inverseAlgorithmFactory.getAlgorithm(4);
    ////auto inverseAlgorithm = inverseAlgorithmFactory.getAlgorithm(2);
    ////auto inverseAlgorithm = inverseAlgorithmFactory.getAlgorithm(4);
    //inverseAlgorithm->operator()(&complexData[0]);
//
    //for (auto value : complexData)
    //    std::cout << std::setw(10) << std::setprecision(5) << value.real() << "\t" << value.imag() << "I\n";

    ///////////////////////////////////
    //jeanbaptiste::windowing::BartlettWindow<std::integral_constant<unsigned, 16>, std::complex<double>> bartlett;

    const double pi = std::acos(-1);
    const double s11 = std::cos(pi/6);
    const double s12 = std::cos(-3*pi/4);

    constexpr double d1 = 1.0 / 6.0 * constants::pi<double>();
    constexpr double d2 = constants::pi<double>();
    constexpr double s31 = jeanbaptiste::basic::cosine<double>(1.0 / 6.0 * constants::pi<double>());
    constexpr double s32 = jeanbaptiste::basic::cosine<double>(-3.0 / 4.0 * constants::pi<double>());

    // pi*n   pi   2*pi*n - N*pi        2*n - N
    // ---- - -- = ------------- = pi * --------
    //   N     2         2*N               2*N

    //constexpr double s23 = jeanbaptiste::basic::sine<double>(2*n - N, 2*N);

    double f1 = std::abs(-1.6);
    constexpr double f2 = jb::basic::abs(-1.6);

    //jw::BartlettWindow<hana::int_c<128>, jbo::Direction_Forward, std::complex<double>> bartlett;
    jw::BartlettWindow<std::integral_constant<int, 16>, jbo::Direction_Forward, std::complex<double>> bartlett;
    bartlett(&complexData[0]);

    int i = 0;

    // std::vector<std::complex<double>> complexData =
    // {
    //     std::complex<double>(-1.0f, 0.0f),
    //     std::complex<double>(-1.0f, 0.0f),
    //     std::complex<double>(-1.0f, 0.0f),
    //     std::complex<double>(-1.0f, 0.0f),
    //     std::complex<double>(1.0f, 0.0f),
    //     std::complex<double>(1.0f, 0.0f),
    //     std::complex<double>(1.0f, 0.0f),
    //     std::complex<double>(1.0f, 0.0f),
    //     std::complex<double>(-1.0f, 0.0f),
    //     std::complex<double>(-1.0f, 0.0f),
    //     std::complex<double>(-1.0f, 0.0f),
    //     std::complex<double>(-1.0f, 0.0f),
    //     std::complex<double>(1.0f, 0.0f),
    //     std::complex<double>(1.0f, 0.0f),
    //     std::complex<double>(1.0f, 0.0f),
    //     std::complex<double>(1.0f, 0.0f)
    // };

    // // Radix-2 DIF TMP
    // std::cout << "fwd radix 2 DIF (TMP):" << std::endl;

    // jeanbaptiste::core::Radix2DIF<16, 1, double, std::complex> test;
    // test.apply(&complexData[0]);

    // jeanbaptiste::basic::IndexSwapper<16, 0U, double, std::complex> indexSwapper;
    // indexSwapper.apply(&complexData[0]);
    
    // jeanbaptiste::Normalization<16, 0, double, std::complex> norm;
    // norm.apply(&complexData[0]);

    // for (auto value : complexData)
    //     std::cout << std::setw(10) << std::setprecision(5) << value.real() << "\t" << value.imag() << "I\n";

    // std::cin.get();

    return 0;
}