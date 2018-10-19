#include <boost/hana.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <complex>
#include "../include/AlgorithmFactory.h"
#include <iostream>
#include <memory>

#include "../include/basic/SineCosine.h"
#include "../include/basic/Ceil.h"
#include "../include/basic/Floor.h"
#include "../include/basic/Round.h"

namespace constants = boost::math::constants;
namespace hana = boost::hana;
namespace jb = jeanbaptiste;
namespace jbo = jeanbaptiste::options;
namespace jw = jeanbaptiste::windowing;

int main()
{
    // std::array<std::complex<double>, 16> complexData =
    // {
    //     std::complex<double>(-1.0f, 0.0f),
	// 	std::complex<double>(-1.0f, 0.0f),
	// 	std::complex<double>(-1.0f, 0.0f),
	// 	std::complex<double>(-1.0f, 0.0f),
	// 	std::complex<double>( 1.0f, 0.0f),
	// 	std::complex<double>( 1.0f, 0.0f),
	// 	std::complex<double>( 1.0f, 0.0f),
	// 	std::complex<double>( 1.0f, 0.0f),
	// 	std::complex<double>(-1.0f, 0.0f),
	// 	std::complex<double>(-1.0f, 0.0f),
	// 	std::complex<double>(-1.0f, 0.0f),
	// 	std::complex<double>(-1.0f, 0.0f),
	// 	std::complex<double>( 1.0f, 0.0f),
	// 	std::complex<double>( 1.0f, 0.0f),
	// 	std::complex<double>( 1.0f, 0.0f),
	// 	std::complex<double>( 1.0f, 0.0f)
    // };

    // std::cout << "FFT:\n";

    // jb::AlgorithmFactory<1, 6, jbo::Radix_2, jbo::Decimation_In_Frequency, jbo::Direction_Forward, 
    //    jbo::Window_None, jbo::Normalization_No, std::complex<double>> algorithmFactory;

    // auto algorithm = algorithmFactory.getAlgorithm(4);
    // algorithm->operator()(&complexData[0]);

    // for (auto value : complexData)
    //    std::cout << std::setw(10) << std::setprecision(5) << value.real() << "\t" << value.imag() << "I\n";

    // std::cin.get();

	/////////////////

	constexpr auto n = 40000;
	auto nc = jeanbaptiste::basic::ceil<double>(n);
	auto nf = jeanbaptiste::basic::floor<double>(n);
	auto nr = jeanbaptiste::basic::round<double>(n);
	
	//int32_t quadrant{0};
	//auto nRed = jeanbaptiste::basic::reduceRange<double>(n, quadrant);

	auto nSin = jeanbaptiste::basic::sine<double>(n);
	auto nCos = jeanbaptiste::basic::cosine<double>(n);

	//int8_t quadrant{0};

	for (float x = -3.1415926535897932; x <= 3.1415926535897932; x += 0.3926990816987241)
	{
		//auto r1 = jeanbaptiste::basic::reduceRange(x, quadrant);
		auto ns1 = jeanbaptiste::basic::sine<double>(x);
		auto nc1 = jeanbaptiste::basic::cosine<double>(x);

		auto ns2 = std::sin(x);
		auto nc2 = std::cos(x);

		int i = 0;
	}

     return 0;
}