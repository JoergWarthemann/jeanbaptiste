#include <boost/hana.hpp>
#include <complex>
#include "../include/AlgorithmFactory.h"
#include <iostream>
#include <memory>

namespace hana = boost::hana;
namespace jb = jeanbaptiste;
namespace jbo = jeanbaptiste::options;

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

    std::cout << "FFT:\n";
    //jb::AlgorithmFactory<1, 6, jbo::Radix_2, jbo::Decimation_In_Frequency, jbo::Direction_Forward, 
	//jb::AlgorithmFactory<1, 6, jbo::Radix_4, jbo::Decimation_In_Time, jbo::Direction_Forward,
    jb::AlgorithmFactory<1, 6, jbo::Radix_Split_2_4, jbo::Decimation_In_Time, jbo::Direction_Forward,  
        jbo::Window_None, std::complex<double>> algorithmFactory;
    //auto algorithm = algorithmFactory.getAlgorithm(4);
	//auto algorithm = algorithmFactory.getAlgorithm(2);
    auto algorithm = algorithmFactory.getAlgorithm(4);
    algorithm->operator()(&complexData[0]);

    for (auto value : complexData)
        std::cout << std::setw(10) << std::setprecision(5) << value.real() << "\t" << value.imag() << "I\n";

    std::cout << "IFFT:\n";
    //jb::AlgorithmFactory<1, 6, jbo::Radix_2, jbo::Decimation_In_Frequency, jbo::Direction_Backward, 
    //jb::AlgorithmFactory<1, 6, jbo::Radix_4, jbo::Decimation_In_Time, jbo::Direction_Backward,
    jb::AlgorithmFactory<1, 6, jbo::Radix_Split_2_4, jbo::Decimation_In_Time, jbo::Direction_Backward,
        jbo::Window_None, std::complex<double>> inverseAlgorithmFactory;
    //auto inverseAlgorithm = inverseAlgorithmFactory.getAlgorithm(4);
    //auto inverseAlgorithm = inverseAlgorithmFactory.getAlgorithm(2);
    auto inverseAlgorithm = inverseAlgorithmFactory.getAlgorithm(4);
    inverseAlgorithm->operator()(&complexData[0]);

    for (auto value : complexData)
        std::cout << std::setw(10) << std::setprecision(5) << value.real() << "\t" << value.imag() << "I\n";

    std::cin.get();


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