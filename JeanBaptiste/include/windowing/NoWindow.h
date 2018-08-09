 #pragma once

#include "../SubTask.h"

namespace jeanbaptiste::windowing
{
    /** Creates an empty window (rectangular) for a specified sample count.
        \param SampleCnt ... The count of samples to be processed in this recursion level (stage)
        \param Complex ... The complex type.
    */
    template <typename SampleCnt,
              typename DirectionFactor,
              typename Complex>
    class NoWindow
        : public SubTask<NoWindow<SampleCnt, Complex>,
                         Complex>
    {
    public:
        void operator()(Complex* data) const
        {}
    };
}