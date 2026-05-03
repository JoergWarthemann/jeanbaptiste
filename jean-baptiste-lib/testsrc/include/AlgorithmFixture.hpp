#ifndef JEANBAPTISTE_TESTING_ALGORITHM_FIXTURE_HPP_
#define JEANBAPTISTE_TESTING_ALGORITHM_FIXTURE_HPP_

#include <complex>
#include <memory>
#include <vector>

#include "IExecutableAlgorithm.hpp"
#include "AlgorithmResultAnalysis.hpp"

namespace jeanbaptiste::testing {

/** Fixture for executable algorithms.
    Provides common data and functions for algorithm tests.
*/
class AlgorithmFixture {
protected:
    using TAlgorithmType = std::unique_ptr<IExecutableAlgorithm<std::complex<double>>>;
    using TDataSetType = std::vector<double>;
    using TComplexDataSetType = std::vector<std::complex<double>>;

	TComplexDataSetType mDataSetA;//mWorkingSet;
	TComplexDataSetType mDataSetB;//mExpectedOutFFT;
	//std::vector<std::complex<double>> mExpectedOutIFFT;

	//AlgorithmResultAnalysis<double> mAlgorithmResult;

public:
    AlgorithmFixture(void) = default;
    virtual ~AlgorithmFixture(void) = default;

    // TODO: Keep verifyAlgorithm to execute and check a single algorithm
    void verifyAlgorithm(TAlgorithmType algorithm)
    {
        algorithm->operator()(mDataSetA);
        //mAlgorithmResult.checkOutput(mWorkingSet, mExpectedOutFFT);
        AlgorithmResultAnalysis::compareComplexDataSets(mDataSetA, mDataSetB);
    }

    template <typename ...Algorithm>
    void verifyAlgorithms(Algorithm ... algorithm)
    {
        verifyAlgorithm(algorithm ...);
    }
    // void runAlgorithms(TAlgorithmType fft, TAlgorithmType ifft)
    // {
    //     fft->operator()(&workingSet_[0]);
    //     algorithmResult_.checkOutput(workingSet_, expectedOutFFT_);

    //     ifft->operator()(&workingSet_[0]);
    //     algorithmResult_.checkOutput(workingSet_, expectedOutIFFT_);
    // }
};

} // namespace jeanbaptiste::testing

#endif // JEANBAPTISTE_TESTING_ALGORITHM_FIXTURE_HPP_
