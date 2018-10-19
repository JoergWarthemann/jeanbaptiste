#pragma once

#include <boost/test/unit_test.hpp>
#include <complex>
#include "../include/AlgorithmResultAnalysis.h"
#include <memory>
//#include "../../JeanBaptiste/include/ExecutableAlgorithm.h"

namespace jb = jeanbaptiste;

class AlgorithmFixture
{
protected:
    using AlgorithmType = std::unique_ptr<jb::ExecutableAlgorithm<std::complex<double>>>;

	std::vector<std::complex<double>> workingSet_;
	std::vector<std::complex<double>> expectedOutFFT_;
	std::vector<std::complex<double>> expectedOutIFFT_;

	Utilities::AlgorithmResultAnalysis<double> algorithmResult_;

public:
    AlgorithmFixture(void) = default;
    virtual ~AlgorithmFixture(void) = default;

    void runAlgorithm(AlgorithmType fft)
    {
        fft->operator()(&workingSet_[0]);
        algorithmResult_.checkOutput(workingSet_, expectedOutFFT_);
    }

    void runAlgorithms(AlgorithmType fft, AlgorithmType ifft)
    {
        fft->operator()(&workingSet_[0]);
        algorithmResult_.checkOutput(workingSet_, expectedOutFFT_);

        ifft->operator()(&workingSet_[0]);
        algorithmResult_.checkOutput(workingSet_, expectedOutIFFT_);
    }
};