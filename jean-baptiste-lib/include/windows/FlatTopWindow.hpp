#ifndef JB_WINDOWS_FLATTOPWINDOW_HPP_
#define JB_WINDOWS_FLATTOPWINDOW_HPP_

#include <algorithm>
#include <array>
#include <numbers>
#include <span>

#include "tools/SineCosine.hpp"
#include "tools/SubTask.hpp"
#include "windows/ExecuteWindowOnComplexData.hpp"

namespace jb::windows {

template <typename SampleCnt, typename Complex>
    requires std::is_integral_v<typename SampleCnt::value_type> &&
    std::is_floating_point_v<typename Complex::value_type>
class FlatTopWindow
    : public tools::SubTask<FlatTopWindow<SampleCnt, Complex>, Complex> {
public:
    /**
     * Fills the internal vector with values that represent a FlatTop window within SampleCnt samples.
     *
     *     5
     *                     .
     *                   .....                                / 2Pi * n \               / 4Pi * n \               / 6Pi * n \               / 8Pi * n \
     *                  .......         w(n) = 1 - 1.93 * cos|  ——————— | + 1.29 * cos |  ——————— | - 0.388 * cos|  ——————— | + 0.028 * cos|  ——————— |
     *                 .........                              \    N    /               \    N    /               \    N    /               \    N    /
     *                ...........
     *               .............
     *     +————————————————————————————————
     *       ........             ........
     *          ...                 ...
     *     0                                N-1
     */
    void operator()(std::span<Complex, SampleCnt::value> data) const
    {
        std::transform(data.begin(), data.end(), mWindowSamples.begin(), data.begin(), ExecuteWindowOnComplexData<Complex>());
    }

private:
    using ValueType = typename Complex::value_type;

    static consteval ValueType createSample(const std::size_t index)
    {
        return 1.0 -
            1.93 * tools::cosine<double>(kTwoPiDividedBySampleCnt * index) +
            1.29 * tools::cosine<double>(kFourPiDividedBySampleCnt * index) -
            0.388 * tools::cosine<double>(kSixPiDividedBySampleCnt * index) +
            0.028 * tools::cosine<double>(kEightPiDividedBySampleCnt * index);
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
    static constexpr double kSixPiDividedBySampleCnt = 3.0 * kTwoPiDividedBySampleCnt;
    static constexpr double kEightPiDividedBySampleCnt = 4.0 * kTwoPiDividedBySampleCnt;
    static constexpr auto mWindowSamples = getWindowSamples();
};

} // namespace jb::windows

#endif // JB_WINDOWS_FLATTOPWINDOW_HPP_