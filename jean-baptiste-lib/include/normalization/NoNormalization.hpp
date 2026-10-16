#ifndef JB_NORMALIZATION_NONORMALIZATION_HPP_
#define JB_NORMALIZATION_NONORMALIZATION_HPP_

#include <span>
#include <type_traits>

#include "tools/SubTask.hpp"

namespace jb::normalization {

/**
 * Does not normalize FFT results.
 * \tparam SampleCnt ... The count of samples to deal with.
 * \tparam Complex ... Complex data type.
 */
template <typename SampleCnt, typename Complex>
    requires std::is_integral_v<typename SampleCnt::value_type> &&
    std::is_floating_point_v<typename Complex::value_type>
class NoNormalization
    : public tools::SubTask<NoNormalization<SampleCnt, Complex>, Complex> {
public:
    /**
     * Does nothing.
     * \param[in] data ... Span of SampleCnt elements of type Complex.
     */
    void operator()(std::span<Complex, SampleCnt::value> data) const
    {}
};

} // namespace jb::normalization
#endif // JB_NORMALIZATION_NONORMALIZATION_HPP_