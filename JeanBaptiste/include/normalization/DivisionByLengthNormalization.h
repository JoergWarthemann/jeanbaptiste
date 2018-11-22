#pragma once

#include <complex>
#include "../SubTask.h"

namespace jeanbaptiste::normalization
{
    /** Normalizes FFT results by dividing them by the length of the original signal (applying the factor 1/N to all samples).
        Applying the factor 1/N to all samples in one domain keeps the energy of the signal. When going back into the other domain
        no normalization should be used to get the original signal correctly normalized again.
        \param SampleCnt ... The count of samples to deal with.
        \param DenominatorShiftFactor ... An additional shift factor that is to be applied on the normalization factors denominator
                                          in real FFT backward mode.
        \param Complex ... Complex data type.
    */
    template<typename SampleCnt,
             typename DenominatorShiftFactor,
             typename Complex>
    class DivisionByLengthNormalization
        : public SubTask<DivisionByLengthNormalization<SampleCnt, DenominatorShiftFactor, Complex>,
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
        /** Normalizes each element of data with 1/N.
            \param[in] data ... Pointer to an array of SampleCnt elements of type Complex.
        */
        void operator()(Complex* data) const
        {
            for (std::size_t i = 0; i < SampleCnt::value; ++i)
            {
                data[i] *= static_cast<typename Complex::value_type>(1 / static_cast<typename Complex::value_type>(getDenominator()));
            }
        }
    };
}