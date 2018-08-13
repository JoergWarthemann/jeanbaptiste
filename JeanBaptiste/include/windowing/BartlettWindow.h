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
    class BartlettWindow
        : public SubTask<BartlettWindow<SampleCnt, Complex>,
                         Complex>
    {
        using ValueType = typename Complex::value_type;

        static constexpr ValueType createSample(const std::size_t index)
        {
            return 1.0 - basic::abs<ValueType>(index - static_cast<ValueType>(kHalfSampleCnt_)) / SampleCnt::value;
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

        static constexpr unsigned kHalfSampleCnt_ = SampleCnt::value >> 1;
        static constexpr auto windowSamples_ = getWindowSamples();

    public:
        void operator()(Complex* data) const
        {
            std::transform(data, data + SampleCnt::value, windowSamples_.begin(),data, ExecuteWindowOnComplexData<Complex>());
        }
    };
}