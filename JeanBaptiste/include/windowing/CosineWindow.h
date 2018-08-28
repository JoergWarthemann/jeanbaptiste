#pragma once

#include "../basic/Abs.h"
#include "../basic/SineCosine.h"
#include <boost/math/constants/constants.hpp>
#include "ExecuteWindowOnComplexData.h"
#include <functional>
#include <iostream>
#include "../SubTask.h"

namespace constants = boost::math::constants;

namespace jeanbaptiste::windowing
{
    template <typename SampleCnt,
              typename Complex>
    class CosineWindow
        : public SubTask<CosineWindow<SampleCnt, Complex>,
                         Complex>
    {
        using ValueType = typename Complex::value_type;

        static constexpr ValueType createSample(const std::size_t index)
        {
            return jeanbaptiste::basic::cosine<double>(kPiDividedBySampleCnt  * index - constants::half_pi<double>());
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

        static constexpr double kPiDividedBySampleCnt = constants::pi<double>() / SampleCnt::value;
        static constexpr auto windowSamples_ = getWindowSamples();

    public:
        /** Fills the internal vector with values that represent a cosine window within SampleCnt samples.
            
            1
                       ...
                    .........                      / Pi * n    Pi \
                  .............         w(n) = cos|  ——————— - —— |
                .................                  \    N       2 /
               ...................
              .....................
            +———————————————————————
            0                      N-1
		*/
        void operator()(Complex* data) const
        {
            std::transform(data, data + SampleCnt::value, windowSamples_.begin(),data, ExecuteWindowOnComplexData<Complex>());
        }
    };
}