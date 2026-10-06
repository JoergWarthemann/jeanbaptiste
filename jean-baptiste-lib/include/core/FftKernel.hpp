/**
 * Compile-time dispatcher for the core complex FFT implementations.
 *
 * Selects and executes the radix-2, radix-4, or split-radix 2/4 kernel for a
 * fixed decimation strategy, transform direction, and complex sample type. The
 * wrapper also places the bit-reversal permutation on the correct side of the
 * kernel: before decimation-in-time transforms and after decimation-in-frequency
 * transforms.
 */

#ifndef JB_CORE_FFTKERNEL_HPP_
#define JB_CORE_FFTKERNEL_HPP_

#include <cstddef>
#include <span>
#include <type_traits>

#include "Options.hpp"
#include "core/Radix2.hpp"
#include "core/Radix4.hpp"
#include "core/RadixSplit24.hpp"
#include "tools/BitReversalIndexSwapping.hpp"

namespace jb::core {

template <typename>
struct AlwaysFalse : std::false_type {};

/**
 * Compile-time FFT kernel dispatcher.
 *
 * Selects the appropriate FFT kernel (Radix-2, Radix-4, or Split-Radix 2/4)
 * based on the template parameters and applies it to the given data span.
 *
 * \tparam Radix ... The FFT radix to use (options::Radix_2, options::Radix_4, or options::Radix_Split_2_4).
 * \tparam Decimation ... The decimation strategy (options::Decimation_In_Time or options::Decimation_In_Frequency).
 * \tparam DirectionFactor ... The direction of the FFT (options::Direction_Forward or options::Direction_Backward).
 * \tparam Complex ... The complex sample type.
 */
template <typename Radix, typename Decimation, typename DirectionFactor, typename Complex>
class FftKernel {
public:
    /**
     * Applies the FFT kernel to the given data span.
     *
     * Allows hana::tuple execution to call it like other subtasks.
     * \tparam Extent ... The extent of the data span for the operator() overload.
     * \param data ... The data span to process.
     */
    template <std::size_t Extent>
    void operator()(std::span<Complex, Extent> data) const
    {
        apply<std::integral_constant<int, static_cast<int>(Extent)>>(data);
    }

    /**
     * Applies the FFT kernel to the given data span.
     *
     * \tparam SampleCnt ... The number of samples in the data span.
     * \param data ... The data span to process.
     */
    template <typename SampleCnt>
    void operator()(std::span<Complex, SampleCnt::value> data) const
    {
        apply<SampleCnt>(data);
    }

    /**
     * Applies the FFT kernel to the given data span.
     *
     * \tparam SampleCnt ... The number of samples in the data span.
     * \param data ... The data span to process.
     */
    template <typename SampleCnt>
    void apply(std::span<Complex, SampleCnt::value> data) const
    {
        if constexpr (std::is_same_v<Radix, options::Radix_2>) {
            applyRadix2<SampleCnt>(data);
        }
        else if constexpr (std::is_same_v<Radix, options::Radix_4>) {
            applyRadix4<SampleCnt>(data);
        }
        else if constexpr (std::is_same_v<Radix, options::Radix_Split_2_4>) {
            applySplitRadix24<SampleCnt>(data);
        }
        else {
            static_assert(AlwaysFalse<Radix>::value, "Unsupported FFT radix option.");
        }
    }

private:
    template <typename SampleCnt>
    void applyRadix2(std::span<Complex, SampleCnt::value> data) const
    {
        if constexpr (std::is_same_v<Decimation, options::Decimation_In_Time>) {
            tools::BitReversalIndexSwapping<SampleCnt, Complex>{}(data);
            Radix2DIT<SampleCnt, DirectionFactor, Complex>{}(data);
        }
        else {
            Radix2DIF<SampleCnt, DirectionFactor, Complex>{}(data);
            tools::BitReversalIndexSwapping<SampleCnt, Complex>{}(data);
        }
    }

    template <typename SampleCnt>
    void applyRadix4(std::span<Complex, SampleCnt::value> data) const
    {
        if constexpr (std::is_same_v<Decimation, options::Decimation_In_Time>) {
            tools::BitReversalIndexSwapping<SampleCnt, Complex>{}(data);
            Radix4DIT<SampleCnt, DirectionFactor, Complex>{}(data);
        }
        else {
            Radix4DIF<SampleCnt, DirectionFactor, Complex>{}(data);
            tools::BitReversalIndexSwapping<SampleCnt, Complex>{}(data);
        }
    }

    template <typename SampleCnt>
    void applySplitRadix24(std::span<Complex, SampleCnt::value> data) const
    {
        if constexpr (std::is_same_v<Decimation, options::Decimation_In_Time>) {
            tools::BitReversalIndexSwapping<SampleCnt, Complex>{}(data);
            RadixSplit24DIT<SampleCnt, DirectionFactor, Complex>{}(data);
        }
        else {
            RadixSplit24DIF<SampleCnt, DirectionFactor, Complex>{}(data);
            tools::BitReversalIndexSwapping<SampleCnt, Complex>{}(data);
        }
    }
};

} // namespace jb::core

#endif // JB_CORE_FFTKERNEL_HPP_