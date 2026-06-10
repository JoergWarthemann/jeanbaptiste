#ifndef JB_WINDOWS_BLACKMANHARRISWINDOW_HPP_
#define JB_WINDOWS_BLACKMANHARRISWINDOW_HPP_

#include <algorithm>
#include <array>
#include <numbers>
#include <span>


#include "tools/Abs.hpp"
#include "tools/SineCosine.hpp"
#include "tools/SubTask.hpp"
#include "windows/ExecuteWindowOnComplexData.hpp"

namespace jb::windows {

template <typename SampleCnt, typename Complex>
class BlackmanHarrisWindow
    : public tools::SubTask<BlackmanHarrisWindow<SampleCnt, Complex>, Complex> {
public:
    /**
     * Fills the internal vector with values that represent a Blackman-Harris window within SampleCnt samples.
     *
     *  1
     *             .
     *          .......                                         / 2Pi * n \                  / 4Pi * n \                  / 6Pi * n \
     *         .........         w(n) = 0.35875 + 0.48829 * cos|  ——————— | + 0.14128 * cos |  ——————— | + 0.01168 * cos |  ——————— |
     *        ...........                                       \    N    /                  \    N    /                  \    N    /
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

    static consteval long modifyIndex(const std::size_t index)
    {
        return index - static_cast<ValueType>(kHalfSampleCnt);
    }

    static consteval ValueType createSample(const std::size_t index)
    {
        return 0.35875 + 0.48829 * tools::cosine<double>(kTwoPiDividedBySampleCnt * modifyIndex(index)) + 0.14128 * tools::cosine<double>(kFourPiDividedBySampleCnt * modifyIndex(index)) + 0.01168 * tools::cosine<double>(kSixPiDividedBySampleCnt * modifyIndex(index));
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

    static constexpr unsigned kHalfSampleCnt = SampleCnt::value >> 1;
    static constexpr double kTwoPiDividedBySampleCnt = 2.0 * std::numbers::pi / SampleCnt::value;
    static constexpr double kFourPiDividedBySampleCnt = 2.0 * kTwoPiDividedBySampleCnt;
    static constexpr double kSixPiDividedBySampleCnt = 3.0 * kTwoPiDividedBySampleCnt;
    static constexpr auto mWindowSamples = getWindowSamples();
};

} // namespace jb::windows

#endif // JB_WINDOWS_BLACKMANHARRISWINDOW_HPP_