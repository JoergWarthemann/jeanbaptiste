#ifndef JEANBAPTISTE_WINDOWING_BLACKMANWINDOW_HPP_
#define JEANBAPTISTE_WINDOWING_BLACKMANWINDOW_HPP_

#include <algorithm>
#include <array>
#include <numbers>
#include <span>

#include "tools/SineCosine.hpp"
#include "tools/SubTask.hpp"
#include "windows/ExecuteWindowOnComplexData.hpp"

namespace jeanbaptiste::windowing {

template <typename SampleCnt, typename Complex>
class BlackmanWindow
    : public SubTask<BlackmanWindow<SampleCnt, Complex>, Complex> {
public:
    /**
     * Fills the internal vector with values that represent a Blackman window within SampleCnt samples.
     *
     *  1
     *             .
     *          .......                                  / 2Pi * n \               / 4Pi * n \
     *         .........         w(n) = 0.42 - 0.5 * cos|  ——————— | + 0.08 * cos |  ——————— |
     *        ...........                                \    N    /               \    N    /
     *       .............
     *     .................
     *  +————————————————————
     *  0                    N-1
     */
    void operator()(std::span<Complex, SampleCnt::value> data) const
    {
        std::transform(data.begin(), data.end(), mWindowSamples.begin(), data.begin(), ExecuteWindowOnComplexData<Complex>());
    }

private:
    using ValueType = typename Complex::value_type;

    static consteval ValueType createSample(const std::size_t index)
    {
        return 0.42 - 0.5 * tools::cosine<double>(kTwoPiDividedBySampleCnt * index) + 0.08 * tools::cosine<double>(kFourPiDividedBySampleCnt * index);
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
    static constexpr double kFourPiDividedBySampleCnt = 2.0 * kTwoPiDividedBySampleCnt;
    static constexpr auto mWindowSamples = getWindowSamples();
};

} // namespace jeanbaptiste::windowing

#endif // JEANBAPTISTE_WINDOWING_BLACKMANWINDOW_HPP_