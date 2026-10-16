/**
 * Real-valued FFT algorithm wrapper.
 *
 * Orchestrates the full real FFT pipeline around the core complex FFT for a
 * selected stage, radix, decimation strategy, direction, window, normalization,
 * and complex sample type. Forward transforms pack 2N real samples, run the
 * internal N-sample complex FFT, then recombine and unpack the spectrum.
 * Inverse transforms apply the inverse recombination before the internal FFT,
 * then unpack the complex samples back into the 2N real-valued layout.
 */

#ifndef JB_REALALGORITHM_HPP_
#define JB_REALALGORITHM_HPP_

#include <cassert>
#include <span>
#include <type_traits>

#include "Options.hpp"
#include "SampleCount.hpp"
#include "core/FftKernel.hpp"
#include "core/RealFft.hpp"
#include "normalization/DivisionByLengthNormalization.hpp"
#include "normalization/NoNormalization.hpp"
#include "normalization/SquareRootNormalization.hpp"
#include "windows/BartlettWindow.hpp"
#include "windows/BlackmanHarrisWindow.hpp"
#include "windows/BlackmanWindow.hpp"
#include "windows/CosineWindow.hpp"
#include "windows/FlatTopWindow.hpp"
#include "windows/HammingWindow.hpp"
#include "windows/NoWindow.hpp"
#include "windows/VonHannWindow.hpp"
#include "windows/WelchWindow.hpp"

namespace jb {

/**
 * Real-valued FFT algorithm wrapper.
 *
 * Defines a tuple of executable sub tasks which belong to an algorithm, e.g. FFT, normalization, bit reversal.
 * \tparam Stage ... The count of stages inside an algorithm. E.g. Stages = 4 -> sample count = 2^4
 * \tparam Radix ... The radix used for the FFT (e.g., Radix_2, Radix_4, Radix_Split_2_4).
 * \tparam Decimation ... The decimation strategy used for the FFT (e.g., Decimation_InTime, Decimation_InFrequency).
 * \tparam Direction ... The direction of the FFT (e.g., Direction_Forward, Direction_Backward).
 * \tparam Window ... The window function applied to the input data (e.g., Window_Hamming, Window_Blackman).
 * \tparam Normalization ... The normalization method applied to the FFT output (e.g., Normalization_No, Normalization_Division_By_Length).
 * \tparam Complex ... The complex data type.
 */
template <std::size_t Stage,
    typename Radix,
    typename Decimation,
    typename Direction,
    typename Window,
    typename Normalization,
    typename Complex>
class RealAlgorithm {
public:
    void operator()(std::span<Complex> data) const
    {
        assert(data.size() == ExternalSampleCnt::value && "The size of the real FFT input data must be equal to the external sample count for the selected radix and stage.");
        applyWithSampleCounts<ExternalSampleCnt, InternalSampleCnt>(data);
    }

private:
    using SampleCounts = SampleCount<Stage, Radix, options::Data_Real>;
    using ExternalSampleCnt = typename SampleCounts::External;
    using InternalSampleCnt = typename SampleCounts::Internal;

    static consteval auto getDirectionValue() noexcept
    {
        if constexpr (std::is_same_v<Direction, options::Direction_Forward>) {
            return std::integral_constant<int, 1>{};
        }
        else {
            return std::integral_constant<int, -1>{};
        }
    }

    template <typename SampleCnt>
    static consteval auto getWindowValue() noexcept
    {
        if constexpr (std::is_same_v<Window, options::Window_Bartlett>) {
            return windows::BartlettWindow<SampleCnt, Complex>{};
        }
        else if constexpr (std::is_same_v<Window, options::Window_Blackman>) {
            return windows::BlackmanWindow<SampleCnt, Complex>{};
        }
        else if constexpr (std::is_same_v<Window, options::Window_BlackmanHarris>) {
            return windows::BlackmanHarrisWindow<SampleCnt, Complex>{};
        }
        else if constexpr (std::is_same_v<Window, options::Window_Cosine>) {
            return windows::CosineWindow<SampleCnt, Complex>{};
        }
        else if constexpr (std::is_same_v<Window, options::Window_FlatTop>) {
            return windows::FlatTopWindow<SampleCnt, Complex>{};
        }
        else if constexpr (std::is_same_v<Window, options::Window_Hamming>) {
            return windows::HammingWindow<SampleCnt, Complex>{};
        }
        else if constexpr (std::is_same_v<Window, options::Window_vonHann>) {
            return windows::VonHannWindow<SampleCnt, Complex>{};
        }
        else if constexpr (std::is_same_v<Window, options::Window_Welch>) {
            return windows::WelchWindow<SampleCnt, Complex>{};
        }
        else {
            return windows::NoWindow<SampleCnt, Complex>{};
        }
    }

    template <typename SampleCnt>
    static consteval auto getNormalizationValue() noexcept
    {
        using DenominatorShiftFactor = std::integral_constant<int, decltype(getDirectionValue())::value == -1 ? 1 : 0>;

        if constexpr (std::is_same_v<Normalization, options::Normalization_Division_By_Length>) {
            return normalization::DivisionByLengthNormalization<SampleCnt, DenominatorShiftFactor, Complex>{};
        }
        else if constexpr (std::is_same_v<Normalization, options::Normalization_Square_Root>) {
            return normalization::SquareRootNormalization<SampleCnt, DenominatorShiftFactor, Complex>{};
        }
        else {
            return normalization::NoNormalization<SampleCnt, Complex>{};
        }
    }

    template <typename ExternalSampleCnt, typename InternalSampleCnt>
    void applyWithSampleCounts(std::span<Complex> data) const
    {
        auto externalData = std::span<Complex, ExternalSampleCnt::value>(data.data(), ExternalSampleCnt::value);
        auto internalData = std::span<Complex, InternalSampleCnt::value>(data.data(), InternalSampleCnt::value);

        decltype(getWindowValue<ExternalSampleCnt>()){}(externalData);
        core::RealFftInputPacking<
            InternalSampleCnt,
            decltype(getDirectionValue()),
            Complex>{}(externalData);

        if constexpr (decltype(getDirectionValue())::value == 1) {
            core::FftKernel<
                Radix,
                Decimation,
                decltype(getDirectionValue()),
                Complex>{}
                .template apply<InternalSampleCnt>(internalData);
            core::RealFftEvenOddRecombination<
                InternalSampleCnt,
                decltype(getDirectionValue()),
                Complex>{}(internalData);
        }
        else {
            core::RealFftEvenOddRecombination<
                InternalSampleCnt,
                decltype(getDirectionValue()),
                Complex>{}(internalData);
            core::FftKernel<
                Radix,
                Decimation,
                decltype(getDirectionValue()),
                Complex>{}
                .template apply<InternalSampleCnt>(internalData);
        }

        core::RealFftOutputUnpacking<
            InternalSampleCnt,
            decltype(getDirectionValue()),
            Complex>{}(externalData);
        decltype(getNormalizationValue<ExternalSampleCnt>()){}(externalData);
    }
};

} // namespace jb

#endif // JB_REALALGORITHM_HPP_