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
    class FlatTopWindow
        : public SubTask<FlatTopWindow<SampleCnt, Complex>,
                         Complex>
    {
        using ValueType = typename Complex::value_type;

        static constexpr ValueType createSample(const std::size_t index)
        {
            return 1.0
                 - 1.93  * jeanbaptiste::basic::cosine<double>(kTwoPiDividedBySampleCnt   * index)
                 + 1.29  * jeanbaptiste::basic::cosine<double>(kFourPiDividedBySampleCnt  * index)
                 - 0.388 * jeanbaptiste::basic::cosine<double>(kSixPiDividedBySampleCnt   * index)
                 + 0.028 * jeanbaptiste::basic::cosine<double>(kEightPiDividedBySampleCnt * index)
                 ;
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

        static constexpr double kTwoPiDividedBySampleCnt   = 2.0 * constants::pi<double>() / SampleCnt::value;
        static constexpr double kFourPiDividedBySampleCnt  = 2.0 * kTwoPiDividedBySampleCnt;
        static constexpr double kSixPiDividedBySampleCnt   = 3.0 * kTwoPiDividedBySampleCnt;
        static constexpr double kEightPiDividedBySampleCnt = 4.0 * kTwoPiDividedBySampleCnt;
        static constexpr auto windowSamples_ = getWindowSamples();

    public:
        /** Fills the internal vector with values that represent a FlatTop window within SampleCnt samples.
            
            5
                           .
                         .....                                / 2Pi * n \               / 4Pi * n \               / 6Pi * n \               / 8Pi * n \
                        .......         w(n) = 1 - 1.93 * cos|  ——————— | + 1.29 * cos |  ——————— | - 0.388 * cos|  ——————— | + 0.028 * cos|  ——————— |
                       .........                              \    N    /               \    N    /               \    N    /               \    N    /
                      ...........
                     .............
            +————————————————————————————————
			 ........             ........
			    ...                 ...
            0                                N-1
		*/
        void operator()(Complex* data) const
        {
            std::transform(data, data + SampleCnt::value, windowSamples_.begin(),data, ExecuteWindowOnComplexData<Complex>());
        }
    };
}