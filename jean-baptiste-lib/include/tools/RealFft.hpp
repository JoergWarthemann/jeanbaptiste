#ifndef JB_TOOLS_REALFFT_HPP_
#define JB_TOOLS_REALFFT_HPP_

#include <span>
#include <type_traits>

#include <boost/math/constants/constants.hpp>

#include "tools/SineCosine.hpp"

namespace constants = boost::math::constants;

namespace jb::tools {

/**
 * Prepares the input data for the real FFT by reorganizing the samples.
 * \param InternalSampleCnt ... The internal sample count for the real FFT.
 * \param DirectionFactor ... The direction of the FFT (forward or inverse).
 * \param Complex ... The complex type used for the FFT data.
 */
template <typename InternalSampleCnt, typename DirectionFactor, typename Complex>
    requires std::is_integral_v<typename InternalSampleCnt::value_type> &&
    std::is_integral_v<typename DirectionFactor::value_type> &&
    std::is_floating_point_v<typename Complex::value_type>
class RealFftInputPacking {
public:
    void operator()(std::span<Complex, 2 * InternalSampleCnt::value> data) const
    {
        if constexpr (DirectionFactor::value == 1) {
            // Go through samples and reorganize them.
            for (auto idxNode = 0U; idxNode < InternalSampleCnt::value; ++idxNode) {
                data[idxNode].real(data[2 * idxNode].real());
                data[idxNode].imag(data[2 * idxNode + 1].real());
            }
        }
        else {
            data[0].imag(data[InternalSampleCnt::value].real());
        }
    }
};

/**
 * Recombines the even and odd frequency components of a packed real FFT.
 *
 * For a forward transform, converts the result of an N-sample complex FFT on
 * packed real input into the spectrum of a 2N-sample real FFT. For an inverse
 * transform, applies the inverse recombination before unpacking the real data.
 *
 * \param InternalSampleCnt ... The internal sample count for the real FFT.
 * \param DirectionFactor ... The direction of the FFT (forward or inverse).
 * \param Complex ... The complex type used for the FFT data.
 */
template <typename InternalSampleCnt, typename DirectionFactor, typename Complex>
    requires std::is_integral_v<typename InternalSampleCnt::value_type> &&
    std::is_integral_v<typename DirectionFactor::value_type> &&
    std::is_floating_point_v<typename Complex::value_type>
class RealFftEvenOddRecombination {
public:
    void operator()(std::span<Complex, InternalSampleCnt::value> data) const
    {
        using ValueType = typename Complex::value_type;

        constexpr auto halfSampleCnt = InternalSampleCnt::value >> 1;
        constexpr Complex twiddleMultiplier(
            static_cast<ValueType>(
                -2.0 *
                tools::sine<ValueType>(constants::pi<ValueType>() / (2.0 * InternalSampleCnt::value)) *
                tools::sine<ValueType>(constants::pi<ValueType>() / (2.0 * InternalSampleCnt::value))),
            static_cast<ValueType>(
                DirectionFactor::value *
                tools::sine<ValueType>(constants::pi<ValueType>() / InternalSampleCnt::value)));

        Complex twiddleFactor(1.0 + twiddleMultiplier.real(), twiddleMultiplier.imag());

        // Omit calculation of DC (sample 0) and Nyquist frequency (sample N/2)
        for (auto idxNode0 = 1U; idxNode0 < halfSampleCnt; ++idxNode0) {
            // idxNode1 = N - n
            const auto idxNode1 = InternalSampleCnt::value - idxNode0;

            // X.r(n) =               ((R(n) + R(N-n)) / 2)		-> temp1
            //        + cos(pi*n/N) * ((I(n) - I(N-n)) / 2)		-> temp2
            //        - sin(pi*n/N) * ((R(n) - R(N-n)) / 2)		-> temp2
            // X.i(n) =               ((I(n) - I(N-n)) / 2)		-> temp1
            //        - sin(pi*n/N) * ((I(n) + I(N-n)) / 2)		-> temp2
            //        - cos(pi*n/N) * ((R(n) - R(N-n)) / 2)		-> temp2

            Complex temp1(
                static_cast<ValueType>(0.5) * (data[idxNode0].real() + data[idxNode1].real()),
                static_cast<ValueType>(0.5) * (data[idxNode0].imag() - data[idxNode1].imag()));
            Complex temp2(
                DirectionFactor::value * static_cast<ValueType>(0.5) * (data[idxNode0].imag() + data[idxNode1].imag()),
                -DirectionFactor::value * static_cast<ValueType>(0.5) * (data[idxNode0].real() - data[idxNode1].real()));
            temp2 *= twiddleFactor;

            data[idxNode0] = temp1 + temp2;
            data[idxNode1] = temp1 - temp2;
            data[idxNode1].imag(-data[idxNode1].imag());

            // Calculate the next transform factor via trigonometric recurrence.
            twiddleFactor += twiddleMultiplier * twiddleFactor;
        }

        // Calculate DC and Nyquist and save them in sample 1.
        auto temp = data[0].real();
        if constexpr (DirectionFactor::value == 1) {
            data[0].real(temp + data[0].imag()); // DC
            data[0].imag(temp - data[0].imag()); // Nyquist
        }
        else {
            data[0].real(static_cast<ValueType>(0.5) * (temp + data[0].imag())); // DC
            data[0].imag(static_cast<ValueType>(0.5) * (temp - data[0].imag())); // Nyquist
        }
    }
};

/**
 * Unpacks packed real FFT data into the full 2N-sample output layout.
 *
 * For a forward transform, expands the packed half-spectrum into a
 * conjugate-symmetric 2N-sample complex spectrum. For an inverse transform,
 * unpacks N complex samples into 2N real-valued samples.
 *
 * \param InternalSampleCnt ... The internal sample count for the real FFT.
 * \param DirectionFactor ... The direction of the FFT (forward or inverse).
 * \param Complex ... The complex type used for the FFT data.
 */
template <typename InternalSampleCnt, typename DirectionFactor, typename Complex>
    requires std::is_integral_v<typename InternalSampleCnt::value_type> &&
    std::is_integral_v<typename DirectionFactor::value_type> &&
    std::is_floating_point_v<typename Complex::value_type>
class RealFftOutputUnpacking {
public:
    void operator()(std::span<Complex, 2 * InternalSampleCnt::value> data) const
    {
        if constexpr (DirectionFactor::value == 1) {
            // Go through samples and reorganize them.
            for (auto idxNode = 1U; idxNode < InternalSampleCnt::value; ++idxNode) {
                data[2 * InternalSampleCnt::value - idxNode].real(data[idxNode].real());
                data[2 * InternalSampleCnt::value - idxNode].imag(-data[idxNode].imag());
            }

            data[InternalSampleCnt::value].real(data[0].imag()); // Move the Nyquist component to the middle of the spectrum.
            data[InternalSampleCnt::value].imag(0.0);
            data[0].imag(0.0);
        }
        else {
            // Go through samples and reorganize them.
            for (auto idxNode = InternalSampleCnt::value; idxNode-- > 0;) {
                const auto packedSample = data[idxNode];
                data[2 * idxNode].real(packedSample.real());
                data[2 * idxNode].imag(0.0);
                data[2 * idxNode + 1].real(packedSample.imag());
                data[2 * idxNode + 1].imag(0.0);
            }
        }
    }
};

} // namespace jb::tools

#endif // JB_TOOLS_REALFFT_HPP_