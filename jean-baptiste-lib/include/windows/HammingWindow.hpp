#ifndef JEANBAPTISTE_WINDOWING_HAMMINGWINDOW_HPP_
#define JEANBAPTISTE_WINDOWING_HAMMINGWINDOW_HPP_

#include <algorithm>
#include <array>
#include <numbers>
#include <span>

#include "tools/SineCosine.hpp"
#include "tools/SubTask.hpp"
#include "windows/ExecuteWindowOnComplexData.hpp"

namespace jeanbaptiste::windowing {

template <typename SampleCnt, typename Complex>
class HammingWindow
    : public SubTask<HammingWindow<SampleCnt, Complex>, Complex> {
public:
    /**
     * Fills the internal vector with values that represent a Hamming window within SampleCnt samples.
     *
     *   1
     *              ...
     *           .........                                   / 2Pi * n \
     *          ...........         w(n) = 0.54 + 0.46 * cos|  ———————  |
     *         .............                                 \    N    /
     *        ...............
     *     .....................
     *   +———————————————————————
     *   0                      N-1
     */
    void operator()(std::span<Complex, SampleCnt::value> data) const
    {
        std::transform(data.begin(), data.end(), mWindowSamples.begin(), data.begin(), ExecuteWindowOnComplexData<Complex>());
    }

private:
    using ValueType = typename Complex::value_type;

    static consteval ValueType createSample(const std::size_t index)
    {
        return 0.54 +
            0.46 * tools::cosine<double>(kTwoPiDividedBySampleCnt * (index - static_cast<ValueType>(kHalfSampleCnt_)));
    }

    template <std::size_t... Indices>
    static consteval auto createWindowSamples(std::index_sequence<Indices...>)
    {
        return std::array<ValueType, sizeof...(Indices)>{
            createSample(Indices)...};
    }

    static consteval auto getWindowSamples(void)
    {
        return createWindowSamples(std::make_index_sequence<SampleCnt::value>{});
    }

    static constexpr double kTwoPiDividedBySampleCnt = 2.0 * std::numbers::pi / SampleCnt::value;
    static constexpr unsigned kHalfSampleCnt_ = SampleCnt::value >> 1;
    static constexpr auto mWindowSamples = getWindowSamples();
};

} // namespace jeanbaptiste::windowing

#endif // JEANBAPTISTE_WINDOWING_HAMMINGWINDOW_HPP_