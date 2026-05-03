#ifndef JEANBAPTISTE_NORMALIZATION_SQUAREROOTNORMALIZATION_HPP_
#define JEANBAPTISTE_NORMALIZATION_SQUAREROOTNORMALIZATION_HPP_

#include <span>
#include <type_traits>

#include "tools/HeronSquareRoot.hpp"
#include "tools/SubTask.hpp"

namespace jeanbaptiste::normalization {

/**
 * Normalizes FFT results with respect to Parseval's identity.
 * Normalized transforms have the property that energy computed in one domain equals energy computed in the transform domain.
 * Applying the factor 1/sqrt(N) to all samples in both domains enables the norms in both domains to be equivalent.
 * \param SampleCnt ... The count of samples to deal with.
 * \param DenominatorShiftFactor ... Additional shift factor applied to the normalization factors denominator in real FFT backward mode.
 * \param Complex ... Complex data type.
*/
template<typename SampleCnt, typename DenominatorShiftFactor, typename Complex>
requires std::is_integral_v<SampleCnt> &&
    std::is_integral_v<DenominatorShiftFactor> &&
    std::is_floating_point_v<typename Complex::value_type>
class SquareRootNormalization
    : public SubTask<SquareRootNormalization<SampleCnt, DenominatorShiftFactor, Complex>, Complex> {
public:
    /**
     * Normalizes each element of data with 1/sqrt(N).
     * \param[in] data ... Span of SampleCnt elements of type Complex.
     */
    void operator()(std::span<Complex, SampleCnt::value> data) const
    {
        for (auto& element : data) {
            element *= static_cast<typename Complex::value_type>(
                1 / tools::squareRoot<typename Complex::value_type>(0, 8, getDenominator()));
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

}

#endif // JEANBAPTISTE_NORMALIZATION_SQUAREROOTNORMALIZATION_HPP_