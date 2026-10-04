#ifndef JB_CORE_FFTKERNEL_HPP_
#define JB_CORE_FFTKERNEL_HPP_

#include <span>
#include <type_traits>

#include "Options.hpp"
#include "core/Radix2.hpp"
#include "core/Radix4.hpp"
#include "core/RadixSplit24.hpp"
#include "tools/BitReversalIndexSwapping.hpp"

namespace jb::core {

template <typename>
inline constexpr bool alwaysFalse = false;

template <typename Radix, typename Decimation, typename DirectionFactor, typename Complex>
class FftKernel {
public:
    template <typename SampleCnt>
    void operator()(std::span<Complex, SampleCnt::value> data) const
    {
        apply<SampleCnt>(data);
    }

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
            static_assert(alwaysFalse<Radix>, "Unsupported FFT radix option.");
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