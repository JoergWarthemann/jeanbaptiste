#ifndef JB_WINDOWS_WELCHWINDOW_HPP_
#define JB_WINDOWS_WELCHWINDOW_HPP_

#include <algorithm>
#include <array>
#include <span>

#include "tools/SubTask.hpp"
#include "windows/ExecuteWindowOnComplexData.hpp"

namespace jb::windows {

template <typename SampleCnt, typename Complex>
class WelchWindow
    : public tools::SubTask<WelchWindow<SampleCnt, Complex>, Complex> {
public:
    /**
     * Fills the internal vector with values that represent a Welch window within SampleCnt samples.
     *
     *   1                                       /     N-1  \ 2
     *             .....                        |  n - ———   |
     *           .........                      |       2    |
     *         .............         w(n) = 1 - |————————————|
     *       .................                  |    N+1     |
     *      ...................                 |    ———     |
     *     .....................                 \    2     /
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
        auto temp = [&index]() constexpr {
            return (index - kHalfSampleCntMinusOne_) * kHalfSampleCntPlusOneReciprocal_;
        };

        return 1.0 - temp() * temp();
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

    static constexpr double kHalfSampleCntMinusOne_ = (SampleCnt::value - 1) / 2.0;
    static constexpr double kHalfSampleCntPlusOneReciprocal_ = 1.0 / ((SampleCnt::value + 1) / 2.0);
    static constexpr auto mWindowSamples = getWindowSamples();
};

} // namespace jb::windows

#endif // JB_WINDOWS_WELCHWINDOW_HPP_