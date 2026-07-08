#ifndef JB_ALGORITHM_HPP_
#define JB_ALGORITHM_HPP_

#include <cassert>
#include <span>
#include <tuple>
#include <type_traits>

#include <boost/hana.hpp>

#include "core/Radix2.hpp"
#include "core/Radix4.hpp"
#include "tools/BitReversalIndexSwapping.hpp"

// #include "core/RadixSplit24.h"
#include "IExecutableAlgorithm.hpp"
#include "Options.hpp"

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

namespace hana = boost::hana;

namespace jb {

/**
 * Defines a tuple of executable sub tasks which belong to an algorithm, e.g. FFT, normalization, bit reversal.
 * \param Stage ... The count of stages inside an algorithm. E.g. Stages = 4 -> sample count = 2^4
 * \param Complex ... The complex data type.
 */
template <std::size_t Stage,
    typename Radix,
    typename Decimation,
    typename Direction,
    typename Window,
    typename Normalization,
    typename Complex>
class Algorithm : public IExecutableAlgorithm<Complex> {
public:
    /**
     * Executes all sub tasks sequentially.
     * \param[in] data ... View to an array of SampleCnt elements of type Complex.
     */
    void operator()(std::span<Complex> data) const override
    {
        if constexpr (std::is_same_v<Radix, options::Radix_4>) {
            assert(data.size() == R4SampleCnt::value && "The size of the input data must be equal to 4^Stage.");
            hana::for_each(tupleOfSubTasks_, [&](const auto& subTask) {
                subTask(std::span<Complex, R4SampleCnt::value>(data.data(), R4SampleCnt::value));
            });
        }
        else {
            assert(data.size() == R2SampleCnt::value && "The size of the input data must be equal to 2^Stage.");
            hana::for_each(tupleOfSubTasks_, [&](const auto& subTask) {
                subTask(std::span<Complex, R2SampleCnt::value>(data.data(), R2SampleCnt::value));
            });
        }
    }

private:
    using R2SampleCnt = std::integral_constant<int, 1 << Stage>;
    using R4SampleCnt = std::integral_constant<int, 1 << (Stage << 1)>;

    /**
     * Creates a value of the selected direction type (FFT or iFFT) at compilation time.
     * \return value ... The selected value.
     */
    static consteval auto getDirectionValue() noexcept
    {
        if constexpr (std::is_same_v<Direction, options::Direction_Forward>) {
            return std::integral_constant<int, 1>{};
        }
        else {
            return std::integral_constant<int, -1>{};
        }
    }

    /**
     * Creates a value of the selected window type at compilation time.
     * \return value ... The selected value.
     */
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

    /**
     * Creates a value of the selected normalization type for radix-2 FFT results at compilation time.
     * \return value ... The selected value.
     */
    static consteval auto getRadix2NormalizationValue() noexcept
    {
        if constexpr (std::is_same_v<Normalization, options::Normalization_Division_By_Length>) {
            return normalization::DivisionByLengthNormalization<R2SampleCnt, std::integral_constant<int, 0>, Complex>{};
        }
        else if constexpr (std::is_same_v<Normalization, options::Normalization_Square_Root>) {
            return normalization::SquareRootNormalization<R2SampleCnt, std::integral_constant<int, 0>, Complex>{};
        }
        else {
            return normalization::NoNormalization<R2SampleCnt, Complex>{};
        }
    }

    /**
     * Creates a tupel of sub task type values belonging to a radix-2 FFT task at compilation time.
     * \return hana::tuple ... A tuple of sub task type values.
     */
    static consteval auto radix2SubTaskTypeValues() noexcept
    {
        if constexpr (std::is_same_v<Decimation, options::Decimation_In_Time>) {
            return hana::make_tuple(
                decltype(getWindowValue<R2SampleCnt>()){},
                tools::BitReversalIndexSwapping<R2SampleCnt, Complex>{},
                core::Radix2DIT<R2SampleCnt, decltype(getDirectionValue()), Complex>{},
                decltype(getRadix2NormalizationValue()){});
        }
        else {
            return hana::make_tuple(
                decltype(getWindowValue<R2SampleCnt>()){},
                core::Radix2DIF<R2SampleCnt, decltype(getDirectionValue()), Complex>{},
                tools::BitReversalIndexSwapping<R2SampleCnt, Complex>{},
                decltype(getRadix2NormalizationValue()){});
        }
    }

    /**
     * Creates a value of the selected normalization type at compilation time.
     * \return value ... The selected value.
     */
    static consteval auto getRadix4NormalizationValue() noexcept
    {
        if constexpr (std::is_same_v<Normalization, options::Normalization_Division_By_Length>) {
            return normalization::DivisionByLengthNormalization<R4SampleCnt, std::integral_constant<int, 0>, Complex>{};
        }
        else if constexpr (std::is_same_v<Normalization, options::Normalization_Square_Root>) {
            return normalization::SquareRootNormalization<R4SampleCnt, std::integral_constant<int, 0>, Complex>{};
        }
        else {
            return normalization::NoNormalization<R4SampleCnt, Complex>{};
        }
    }

    /**
     * Creates a tupel of sub task type values belonging to a radix-4 FFT task at compilation time.
     * \return hana::tuple ... A tuple of sub task type values.
     */
    static consteval auto radix4SubTaskTypeValues() noexcept
    {
        if constexpr (std::is_same_v<Decimation, options::Decimation_In_Time>) {
            return hana::make_tuple(
                decltype(getWindowValue<R4SampleCnt>()){},
                tools::BitReversalIndexSwapping<R4SampleCnt, Complex>{},
                core::Radix4DIT<R4SampleCnt, decltype(getDirectionValue()), Complex>{},
                decltype(getRadix4NormalizationValue()){});
        }
        else {
            return hana::make_tuple(
                decltype(getWindowValue<R4SampleCnt>()){},
                core::Radix4DIF<R4SampleCnt, decltype(getDirectionValue()), Complex>{},
                tools::BitReversalIndexSwapping<R4SampleCnt, Complex>{},
                decltype(getRadix4NormalizationValue()){});
        }
    }

    /** Creates a tupel of sub task type values belonging to a split radix 2-4 task at compilation time.
        \return hana::tuple ... A tuple of sub task type values.
    */
    // TODO: add radix split 2-4 case.
    // static consteval auto radixSplit24SubtaskTypeValues() noexcept
    // {
    //     if constexpr (std::is_same_v<Decimation, jbo::Decimation_In_Time>) {
    //         return hana::make_tuple(
    //             decltype(getWindowValue()){},
    //             basic::BitReversalIndexSwapping<R2SampleCnt, Complex>{},
    //             core::RadixSplit24DIT<R2SampleCnt, decltype(getDirectionValue()), Complex>{},
    //             decltype(getRadix2NormalizationValue()){});
    //     }

    //     return hana::make_tuple(
    //         decltype(getWindowValue()){},
    //         core::RadixSplit24DIF<R2SampleCnt, decltype(getDirectionValue()), Complex>{},
    //         basic::BitReversalIndexSwapping<R2SampleCnt, Complex>{},
    //         decltype(getRadix2NormalizationValue()){});
    // }

    /**
     * Creates a tupel of sub task type values at compilation time. Sub tasks belong to a FFT task.
     * \return hana::tuple ... A tuple of sub task type values.
     */
    static consteval auto createTupleOfSubTaskTypeValues() noexcept
    {
        if constexpr (std::is_same_v<Radix, options::Radix_2>) {
            return radix2SubTaskTypeValues();
        }
        else if constexpr (std::is_same_v<Radix, options::Radix_4>) {
            return radix4SubTaskTypeValues();
        }

        // TODO: add split radix 2-4 case.
        // return radixSplit24SubtaskTypeValues();
    }

    using SubTaskTypes = decltype(createTupleOfSubTaskTypeValues());
    SubTaskTypes tupleOfSubTasks_;
};

} // namespace jb

#endif // JB_ALGORITHM_HPP_