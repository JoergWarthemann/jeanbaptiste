#ifndef JEANBAPTISTE_ALGORITHM_HPP_
#define JEANBAPTISTE_ALGORITHM_HPP_

#include <cassert>
#include <span>
#include <tuple>
#include <type_traits>

#include <boost/hana.hpp>

#include "tools/BitReversalIndexSwapping.hpp"
#include "core/Radix2.hpp"
//#include "core/Radix4.h"
//#include "core/RadixSplit24.h"
#include "IExecutableAlgorithm.hpp"
#include "Options.hpp"

// #include "normalization/DivisionByLengthNormalization.h"
// #include "normalization/NoNormalization.h"
// #include "normalization/SquareRootNormalization.h"
// #include "windowing/BartlettWindow.h"
// #include "windowing/BlackmanHarrisWindow.h"
// #include "windowing/BlackmanWindow.h"
// #include "windowing/CosineWindow.h"
// #include "windowing/FlatTopWindow.h"
// #include "windowing/HammingWindow.h"
// #include "windowing/NoWindow.h"
// #include "windowing/VonHannWindow.h"
// #include "windowing/WelchWindow.h"

namespace hana = boost::hana;
namespace jbo = jeanbaptiste::options;

namespace jeanbaptiste {

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
    void operator()(std::span<Complex, SampleCnt::value> data) const override
    {
        hana::for_each(tupleOfSubTasks_, [&](const auto& subTask)
        {
            subTask(data);
        });
    }

private:
    using R2SampleCnt = std::integral_constant<int, 1 << Stage>;
    using R4SampleCnt = std::integral_constant<int, 1 << (Stage << 1)>;
    using SubTaskTypes = decltype(createTupleOfSubTaskTypeValues());

    /** Creates a value of the selected direction type at compilation time.
        \return value ... The selected value.
    */
    static consteval auto getDirectionValue() noexcept
    {
        if constexpr (std::is_same_v<Direction, jbo::Direction_Forward>)
            return std::integral_constant<int, 1>{};
        else
            return std::integral_constant<int, -1>{};
    }

    /** Creates a value of the selected window type at compilation time.
        \return value ... The selected value.
    */
    static consteval auto getWindowValue() noexcept
    {
        if constexpr (std::is_same_v<Window, jbo::Window_Bartlett>)
            return windowing::BartlettWindow<R2SampleCnt, Complex>{};
        else if constexpr (std::is_same_v<Window, jbo::Window_Blackman>)
            return windowing::BlackmanWindow<R2SampleCnt, Complex>{};
        else if constexpr (std::is_same_v<Window, jbo::Window_BlackmanHarris>)
            return windowing::BlackmanHarrisWindow<R2SampleCnt, Complex>{};
        else if constexpr (std::is_same_v<Window, jbo::Window_Cosine>)
            return windowing::CosineWindow<R2SampleCnt, Complex>{};
        else if constexpr (std::is_same_v<Window, jbo::Window_FlatTop>)
            return windowing::FlatTopWindow<R2SampleCnt, Complex>{};
        else if constexpr (std::is_same_v<Window, jbo::Window_Hamming>)
            return windowing::HammingWindow<R2SampleCnt, Complex>{};
        else if constexpr (std::is_same_v<Window, jbo::Window_vonHann>)
            return windowing::VonHannWindow<R2SampleCnt, Complex>{};
        else if constexpr (std::is_same_v<Window, jbo::Window_Welch>)
            return windowing::WelchWindow<R2SampleCnt, Complex>{};
        else
            return windowing::NoWindow<R2SampleCnt, Complex>{};
    }

    /** Creates a value of the selected normalization type at compilation time.
        \return value ... The selected value.
    */
    static consteval auto getRadix2NormalizationValue() noexcept
    {
        if constexpr (std::is_same_v<Normalization, jbo::Normalization_Division_By_Length>)
            return normalization::DivisionByLengthNormalization<R2SampleCnt, std::integral_constant<int, 0>, Complex>{};
        else if constexpr (std::is_same_v<Normalization, jbo::Normalization_Square_Root>)
            return normalization::SquareRootNormalization<R2SampleCnt, std::integral_constant<int, 0>, Complex>{};
        else
            return normalization::NoNormalization<R2SampleCnt, Complex>{};
        }

    /** Creates a tupel of sub task type values belonging to a radix 2 task at compilation time.
        \return hana::tuple ... A tuple of sub task type values.
    */
    static consteval auto radix2SubTaskTypeValues() noexcept
    {
        if constexpr (std::is_same_v<Decimation, jbo::Decimation_In_Time>)
            return hana::make_tuple(
                decltype(getWindowValue()){},
                basic::BitReversalIndexSwapping<R2SampleCnt, Complex>{},
                core::Radix2DIT<R2SampleCnt, decltype(getDirectionValue()), Complex>{},
                decltype(getRadix2NormalizationValue()){});
        else
            return hana::make_tuple(
                decltype(getWindowValue()){},
                core::Radix2DIF<R2SampleCnt, decltype(getDirectionValue()), Complex>{},
                basic::BitReversalIndexSwapping<R2SampleCnt, Complex>{},
                decltype(getRadix2NormalizationValue()){});
    }

    /** Creates a value of the selected normalization type at compilation time.
        \return value ... The selected value.
    */
    static consteval auto getRadix4NormalizationValue() noexcept
    {
        if constexpr (std::is_same_v<Normalization, jbo::Normalization_Division_By_Length>)
            return normalization::DivisionByLengthNormalization<R4SampleCnt, std::integral_constant<int, 0>, Complex>{};
        else if constexpr (std::is_same_v<Normalization, jbo::Normalization_Square_Root>)
            return normalization::SquareRootNormalization<R4SampleCnt, std::integral_constant<int, 0>, Complex>{};
        else
            return normalization::NoNormalization<R4SampleCnt, Complex>{};
    }

    /** Creates a tupel of sub task type values belonging to a radix 4 task at compilation time.
        \return hana::tuple ... A tuple of sub task type values.
    */
    static consteval auto radix4SubTaskTypeValues() noexcept
    {
        if constexpr (std::is_same_v<Decimation, jbo::Decimation_In_Time>)
            return hana::make_tuple(
                decltype(getWindowValue()){},
                basic::BitReversalIndexSwapping<R4SampleCnt, Complex>{},
                core::Radix4DIT<R4SampleCnt, decltype(getDirectionValue()), Complex>{},
                decltype(getRadix4NormalizationValue()){});
        else
            return hana::make_tuple(
                decltype(getWindowValue()){},
                core::Radix4DIF<R4SampleCnt, decltype(getDirectionValue()), Complex>{},
                basic::BitReversalIndexSwapping<R4SampleCnt, Complex>{},
                decltype(getRadix4NormalizationValue()){});
    }

    /** Creates a tupel of sub task type values belonging to a split radix 2-4 task at compilation time.
        \return hana::tuple ... A tuple of sub task type values.
    */
    static consteval auto radixSplit24SubtaskTypeValues() noexcept
    {
        if constexpr (std::is_same_v<Decimation, jbo::Decimation_In_Time>)
            return hana::make_tuple(
                decltype(getWindowValue()){},
                basic::BitReversalIndexSwapping<R2SampleCnt, Complex>{},
                core::RadixSplit24DIT<R2SampleCnt, decltype(getDirectionValue()), Complex>{},
                decltype(getRadix2NormalizationValue()){});
        else
            return hana::make_tuple(
                decltype(getWindowValue()){},
                core::RadixSplit24DIF<R2SampleCnt, decltype(getDirectionValue()), Complex>{},
                basic::BitReversalIndexSwapping<R2SampleCnt, Complex>{},
                decltype(getRadix2NormalizationValue()){});
    }

    /** Creates a tupel of sub task type values at compilation time. Sub tasks belong to a FFT task.
        \return hana::tuple ... A tuple of sub task type values.
    */
    static consteval auto createTupleOfSubTaskTypeValues() noexcept
    {
        if constexpr (std::is_same_v<Radix, jbo::Radix_2>)
            return radix2SubTaskTypeValues();
        else if constexpr (std::is_same_v<Radix, jbo::Radix_4>)
            return radix4SubTaskTypeValues();
        else
            return radixSplit24SubtaskTypeValues();
    }

    SubTaskTypes tupleOfSubTasks_;
};

}