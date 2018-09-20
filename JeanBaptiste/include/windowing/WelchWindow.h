#pragma once

#include "../basic/Abs.h"
#include "../basic/SineCosine.h"
#include "ExecuteWindowOnComplexData.h"
#include <functional>
#include <iostream>
#include "../SubTask.h"

namespace jeanbaptiste::windowing
{
    template <typename SampleCnt,
              typename Complex>
    class WelchWindow
        : public SubTask<WelchWindow<SampleCnt, Complex>,
                         Complex>
    {
        using ValueType = typename Complex::value_type;

        static constexpr ValueType createSample(const std::size_t index)
        {
            auto temp = [&index]() constexpr
            {
                return (index - kHalfSampleCntMinusOne_) * kHalfSampleCntPlusOneReciprocal_;
            };

            return 1.0 - temp() * temp();
        }

        template<std::size_t... Indices>
        static constexpr auto createWindowSamples(std::index_sequence<Indices...>)
        {
            return std::array<ValueType, sizeof...(Indices)>
            {
                createSample(Indices)...
            };
        }

        static constexpr auto getWindowSamples(void)
        {
            return createWindowSamples(std::make_index_sequence<SampleCnt::value>{});
        }

        static constexpr double kHalfSampleCntMinusOne_ = (SampleCnt::value - 1) / 2.0;
        static constexpr double kHalfSampleCntPlusOneReciprocal_ = 1.0 / ((SampleCnt::value + 1) / 2.0);
        static constexpr auto windowSamples_ = getWindowSamples();

    public:
        /** Fills the internal vector with values that represent a Welch window within SampleCnt samples.
            
            1                                       /     N-1  \ 2
                      .....                        |  n - ———   |
                    .........                      |       2    |
                  .............         w(n) = 1 - |————————————|
                .................                  |    N+1     |
               ...................                 |    ———     |
              .....................                 \    2     /
            +———————————————————————
            0                      N-1
		*/
        void operator()(Complex* data) const
        {
            std::transform(data, data + SampleCnt::value, windowSamples_.begin(),data, ExecuteWindowOnComplexData<Complex>());
        }
    };
}