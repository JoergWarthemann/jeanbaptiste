#ifndef JEANBAPTISTE_WINDOWING_COSINEWINDOW_HPP_
#define JEANBAPTISTE_WINDOWING_COSINEWINDOW_HPP_

#include <algorithm>
#include <array>
#include <span>

#include <boost/math/constants/constants.hpp>

#include "tools/SineCosine.hpp"
#include "tools/SubTask.hpp"
#include "windows/ExecuteWindowOnComplexData.hpp"

namespace jeanbaptiste::windowing {

template <typename SampleCnt, typename Complex>
class CosineWindow
    : public SubTask<CosineWindow<SampleCnt, Complex>, Complex> {
public:
    /**
     * Fills the internal vector with values that represent a cosine window within SampleCnt samples.
     *
     *
     *  1
     *             ...
     *          .........                      / Pi * n    Pi \
     *        .............         w(n) = cos|  ——————— - —— |
     *      .................                  \    N       2 /
     *     ...................
     *    .....................
     *  +———————————————————————
     *  0                      N-1
     */
    void operator()(std::span<Complex, SampleCnt::value> data) const
    {
        std::transform(data.begin(), data.end(), mWindowSamples.begin(), data.begin(), ExecuteWindowOnComplexData<Complex>());
    }

private:
    using ValueType = typename Complex::value_type;

    static consteval ValueType createSample(const std::size_t index)
    {
        return tools::cosine<double>(kPiDividedBySampleCnt * index - boost::math::constants::half_pi<double>());
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

    static constexpr double kPiDividedBySampleCnt = boost::math::constants::pi<double>() / SampleCnt::value;
    static constexpr auto mWindowSamples = getWindowSamples();
};

} // namespace jeanbaptiste::windowing

#endif // JEANBAPTISTE_WINDOWING_COSINEWINDOW_HPP_