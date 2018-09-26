#pragma once

#include <complex>
#include "../SubTask.h"

namespace jeanbaptiste::normalization
{
    /** Does not normalize FFT results.
        \param SampleCnt ... The count of samples to deal with.
        \param Complex ... Complex data type.
    */
    template<typename SampleCnt,
             typename Complex>
    class NoNormalization
        : public SubTask<NoNormalization<SampleCnt, Complex>,
                         Complex>
    {
    public:
        /** Does nothing.
            \param[in] Complex* ... Pointer to an array of SampleCnt elements of type Complex.
        */
        void operator()(Complex*) const
        {}
    };
}