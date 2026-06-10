#ifndef JB_NORMALIZATION_DIVISIONBYLENGTHNORMALIZATION_HPP_
#define JB_NORMALIZATION_DIVISIONBYLENGTHNORMALIZATION_HPP_

#include <span>
#include <type_traits>

#include "tools/SubTask.hpp"

namespace jb::normalization {

/**
 * Normalizes FFT results by dividing them by the length of the original signal (applying the factor 1/N to all samples).
 * Applying the factor 1/N to all samples in one domain keeps the energy of the signal. When going back into the other domain
 * no normalization should be used to get the original signal correctly normalized again.
 * \param SampleCnt ... The count of samples to deal with.
 * \param DenominatorShiftFactor ... Additional shift factor applied to the normalization factors denominator in real FFT backward mode.
 * \param Complex ... Complex data type.
 */
template <typename SampleCnt, typename DenominatorShiftFactor, typename Complex>
    requires std::is_integral_v<typename SampleCnt::value_type> &&
    std::is_integral_v<typename DenominatorShiftFactor::value_type> &&
    std::is_floating_point_v<typename Complex::value_type>
class DivisionByLengthNormalization
    : public tools::SubTask<DivisionByLengthNormalization<SampleCnt, DenominatorShiftFactor, Complex>, Complex> {
public:
    /**
     * Normalizes each element of data with 1/N.
     * \param[in] data ... Span of SampleCnt elements of type Complex.
     */
    void operator()(std::span<Complex, SampleCnt::value> data) const
    {
        for (auto& element : data) {
            element *= static_cast<typename Complex::value_type>(
                1 / static_cast<typename Complex::value_type>(getDenominator()));
        }
    }

private:
    /**
     * Calculates the value used as denominator when normalizing data.
     * \return std::size_t ... The denominator
     */
    static constexpr std::size_t getDenominator(void)
    {
        return SampleCnt::value >> DenominatorShiftFactor::value;
    }
};

} // namespace jb::normalization

#endif // JB_NORMALIZATION_DIVISIONBYLENGTHNORMALIZATION_HPP_