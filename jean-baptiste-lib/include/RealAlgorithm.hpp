#ifndef JB_REALALGORITHM_HPP_
#define JB_REALALGORITHM_HPP_

#include <cassert>
#include <span>
#include <type_traits>

#include "Options.hpp"
#include "SampleCount.hpp"
#include "core/FftKernel.hpp"
#include "normalization/DivisionByLengthNormalization.hpp"
#include "normalization/NoNormalization.hpp"
#include "normalization/SquareRootNormalization.hpp"
#include "tools/RealFft.hpp"
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
        tools::RealFftInputPacking<InternalSampleCnt, decltype(getDirectionValue()), Complex>{}(externalData);

        if constexpr (decltype(getDirectionValue())::value == 1) {
            core::FftKernel<Radix, Decimation, decltype(getDirectionValue()), Complex>{}.template apply<InternalSampleCnt>(internalData);
            tools::RealFftEvenOddRecombination<InternalSampleCnt, decltype(getDirectionValue()), Complex>{}(internalData);
        }
        else {
            tools::RealFftEvenOddRecombination<InternalSampleCnt, decltype(getDirectionValue()), Complex>{}(internalData);
            core::FftKernel<Radix, Decimation, decltype(getDirectionValue()), Complex>{}.template apply<InternalSampleCnt>(internalData);
        }

        tools::RealFftOutputUnpacking<InternalSampleCnt, decltype(getDirectionValue()), Complex>{}(externalData);
        decltype(getNormalizationValue<ExternalSampleCnt>()){}(externalData);
    }
};

} // namespace jb

#endif // JB_REALALGORITHM_HPP_