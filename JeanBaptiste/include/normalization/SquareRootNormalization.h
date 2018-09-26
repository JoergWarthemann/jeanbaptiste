#pragma once

#include "../basic/HeronSquareRoot.h"
#include <complex>
#include "../SubTask.h"

namespace jeanbaptiste::normalization
{
    /** Normalizes FFT results. Divides each frequency value by the square root of the sample length.
        \param SampleCnt ... The count of samples to deal with.
        \param DenominatorShiftFactor ... An additional shift factor that is to be applied on the normalization factors denominator
                                          in real FFT backward mode.
        \param Complex ... Complex data type.
    */
    template<typename SampleCnt,
             typename DenominatorShiftFactor,
             typename Complex>
    class SquareRootNormalization
        : public SubTask<SquareRootNormalization<SampleCnt, DenominatorShiftFactor, Complex>,
                         Complex>
    {
        /** Calculates the value used as denominator when normalizing data.
            \return std::size_t ... The denominator
        */
        static constexpr std::size_t getDenominator(void)
        {
            return SampleCnt::value >> DenominatorShiftFactor::value;
        }

    public:
        /** Normalizes each element of data with 1/sqrt(N).
            \param[in] data ... Pointer to an array of SampleCnt elements of type Complex.
        */
        void operator()(Complex* data) const
        {
            for (std::size_t i = 0; i < SampleCnt::value; ++i)
            {
                data[i] *= static_cast<typename Complex::value_type>(1 / basic::squareRoot<typename Complex::value_type>(0, 8, getDenominator()));
            }
        }
    };
}