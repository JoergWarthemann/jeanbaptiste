#ifndef JB_WINDOWS_BARTLETTWINDOW_HPP_
#define JB_WINDOWS_BARTLETTWINDOW_HPP_

#include <algorithm>
#include <array>
#include <span>

#include "tools/Abs.hpp"
#include "tools/SubTask.hpp"
#include "windows/ExecuteWindowOnComplexData.hpp"

namespace jb::windows {

template <typename SampleCnt, typename Complex>
class BartlettWindow
    : public tools::SubTask<BartlettWindow<SampleCnt, Complex>, Complex> {
public:
    /**
     * Fills the internal vector with values that represent a Bartlett window within SampleCnt samples.
     *
     *  1
     *                .                                 N
     *              .....                         | n - — |
     *            .........                             2
     *          .............          w(n) = 1 - —————————
     *        .................                       N
     *      .....................                     —
     *    .........................                   2
     *  +———————————————————————————
     *  0                        N-1
     */
    void operator()(std::span<Complex, SampleCnt::value> data) const
    {
        std::transform(data.begin(), data.end(), mWindowSamples.begin(), data.begin(), ExecuteWindowOnComplexData<Complex>());
    }

private:
    using ValueType = typename Complex::value_type;

    static consteval ValueType createSample(const std::size_t index)
    {
        return 1.0 - tools::abs<ValueType>(index - static_cast<ValueType>(kHalfSampleCnt)) / kHalfSampleCnt;
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
    static constexpr auto mWindowSamples = getWindowSamples();
};

} // namespace jb::windows

#endif // JB_WINDOWS_BARTLETTWINDOW_HPP_