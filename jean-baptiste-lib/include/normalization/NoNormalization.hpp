#ifndef JEANBAPTISTE_NORMALIZATION_NONORMALIZATION_HPP_
#define JEANBAPTISTE_NORMALIZATION_NONORMALIZATION_HPP_

#include <span>
#include <type_traits>

#include "tools/SubTask.hpp"

namespace jeanbaptiste::normalization {

/**
 * Does not normalize FFT results.
 * \param SampleCnt ... The count of samples to deal with.
 * \param Complex ... Complex data type.
 */
template<typename SampleCnt, typename Complex>
requires std::is_integral_v<SampleCnt> &&
    std::is_floating_point_v<typename Complex::value_type>
class NoNormalization
    : public SubTask<NoNormalization<SampleCnt, Complex>, Complex> {
public:
    /**
     * Does nothing.
     * \param[in] data ... Span of SampleCnt elements of type Complex.
     */
    void operator()(std::span<Complex, SampleCnt::value> data) const
    {}
};

}
#endif // JEANBAPTISTE_NORMALIZATION_NONORMALIZATION_HPP_