#pragma once

#include "basic/HeronSquareRoot.h"
#include <complex>

namespace jeanbaptiste
{
    /** Normalizes the FFT results.
        \param SampleCnt ... The count of samples to deal with.
        \param DenominatorShiftFactor ... An additional shift factor that is to be applied on the normalization factors denominator
                                          in real FFT backward mode.
        \param T ... Data type of the samples.
        \param Complex ... Complex data type.
    */
    template<typename SampleCnt,
             unsigned DenominatorShiftFactor,
             typename T = double,
             template<typename...> class Complex = std::complex>
    class Normalization
    {
        enum { kDenominator_ = SampleCnt::value >> DenominatorShiftFactor };

    public:
        // Normalize each sample with 1/sqrt(N) dependent on FFT mode (real, complex) and direction (forward, backward).
        void apply(Complex<T>* data)
        {
            for (auto i = 0; i < SampleCnt::value; ++i)
            {
                data[i] *= static_cast<T>(1 / basic::squareRoot<T>(0, 8, kDenominator_));
            }
        }
    };

    /** Does nothing.
        @param T ... Data type of the samples.
        @param Complex ... Complex data type.
    */
    template<typename T = double,
             template<typename...> class Complex = std::complex>
    class NoNormalization
    {
    public:
        void apply(Complex<T>* data)
        { }
    };
}