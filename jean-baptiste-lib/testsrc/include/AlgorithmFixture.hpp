#ifndef JEANBAPTISTE_TESTING_ALGORITHM_FIXTURE_HPP_
#define JEANBAPTISTE_TESTING_ALGORITHM_FIXTURE_HPP_

#include <complex>
#include <memory>

#include "IExecutableAlgorithm.hpp"
#include "AlgorithmResultAnalysis.hpp"

//#include "../include/AlgorithmResultAnalysis.h"
//#include "../../JeanBaptiste/include/ExecutableAlgorithm.h"

namespace jeanbaptiste::testing {

class AlgorithmFixture
{
protected:
    using AlgorithmType = std::unique_ptr<IExecutableAlgorithm<std::complex<double>>>;

	std::vector<std::complex<double>> workingSet_;
	std::vector<std::complex<double>> expectedOutFFT_;
	std::vector<std::complex<double>> expectedOutIFFT_;

	Utilities::AlgorithmResultAnalysis<double> algorithmResult_;

public:
    AlgorithmFixture(void) = default;
    virtual ~AlgorithmFixture(void) = default;

    // TODO: Keep runAlgorithm to execute and check a single algorithm
    void runAlgorithm(AlgorithmType algorithm)
    {
        // TODO: do only use a span on workingSet_.
        algorithm->operator()(&workingSet_[0]);
        algorithmResult_.checkOutput(workingSet_, expectedOutFFT_);
    }

    template <typename ...Algorithm>
    void runAlgorithms(Algorithm ... algorithm)
    {
        runAlgorithm(algorithm ...);
    }
    // void runAlgorithms(AlgorithmType fft, AlgorithmType ifft)
    // {
    //     fft->operator()(&workingSet_[0]);
    //     algorithmResult_.checkOutput(workingSet_, expectedOutFFT_);

    //     ifft->operator()(&workingSet_[0]);
    //     algorithmResult_.checkOutput(workingSet_, expectedOutIFFT_);
    // }
};

} // namespace jeanbaptiste::testing

#endif // JEANBAPTISTE_TESTING_ALGORITHM_FIXTURE_HPP_
