#pragma once

#include "basic/BitReversalIndexSwapping.h"
#include <boost/hana.hpp>
#include <cassert>
#include "core/Radix2.h"
#include "core/Radix4.h"
#include "core/RadixSplit24.h"
#include "ExecutableAlgorithm.h"
#include "normalization/DivisionByLengthNormalization.h"
#include "normalization/NoNormalization.h"
#include "normalization/SquareRootNormalization.h"
#include "Options.h"
#include "windowing/BartlettWindow.h"
#include "windowing/BlackmanHarrisWindow.h"
#include "windowing/BlackmanWindow.h"
#include "windowing/CosineWindow.h"
#include "windowing/FlatTopWindow.h"
#include "windowing/HammingWindow.h"
#include "windowing/NoWindow.h"
#include "windowing/VonHannWindow.h"
#include "windowing/WelchWindow.h"

namespace hana = boost::hana;
namespace jbo = jeanbaptiste::options;

namespace jeanbaptiste
{
    /** Defines a tuple of executable sub tasks which belong to a FFT task, e.g. FFT, normalization, bit reversal.
        \param Stage ... The count of stages inside an FFT algorithm. E.g. Stages = 4 -> sample count = 2^4
        \param Complex ... The complex data type.
    */
    template <typename Stage,
              typename Radix,
              typename Decimation,
              typename Direction,
              typename Window,
              typename Normalization,
              typename Complex>
    class Algorithm
        : public ExecutableAlgorithm<Complex>
    {
        /** Calculates the number of samples that can be processed by this algorithm.
        */
        static constexpr auto getNumberOfSamples(void)
        {
            if constexpr (hana::typeid_(Radix{}) == hana::type<jbo::Radix_4>{})
                return std::integral_constant<int, 1 << (Stage::value << 1)>{};
            else
                return std::integral_constant<int, 1 << Stage::value>{};
        }

        /** Calculates the number of frequencies that this algorithm puts out.
            Shannon-Nyquist: half the sample count
        */
        static constexpr auto getNumberOfFrequencies(void)
        {
            return getNumberOfSamples() >> 1;
        }

        /** Creates a value of the selected direction type at compilation time.
            \return value ... The selected value.
        */
        static constexpr auto getDirectionValue(void)
        {
            return hana::if_(
                hana::typeid_(Direction{}) == hana::type<jbo::Direction_Forward>{},
                    std::integral_constant<int, 1>{},
                    std::integral_constant<int, -1>{});
        }

        /** Creates a value of the selected window type at compilation time.
            \return value ... The selected value.
        */
        static constexpr auto getWindowValue(void)
        {
            return
                hana::if_(hana::typeid_(Window{}) == hana::type<jbo::Window_Bartlett>{},
                    windowing::BartlettWindow<
                        typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                        Complex>{},
                hana::if_(hana::typeid_(Window{}) == hana::type<jbo::Window_Blackman>{},
                    windowing::BlackmanWindow<
                        typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                        Complex>{},
                hana::if_(hana::typeid_(Window{}) == hana::type<jbo::Window_BlackmanHarris>{},
                    windowing::BlackmanHarrisWindow<
                        typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                        Complex>{},
                hana::if_(hana::typeid_(Window{}) == hana::type<jbo::Window_Cosine>{},
                    windowing::CosineWindow<
                        typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                        Complex>{},
                hana::if_(hana::typeid_(Window{}) == hana::type<jbo::Window_FlatTop>{},
                    windowing::FlatTopWindow<
                        typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                        Complex>{},
                hana::if_(hana::typeid_(Window{}) == hana::type<jbo::Window_Hamming>{},
                    windowing::HammingWindow<
                        typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                        Complex>{},
                hana::if_(hana::typeid_(Window{}) == hana::type<jbo::Window_vonHann>{},
                    windowing::VonHannWindow<
                        typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                        Complex>{},
                hana::if_(hana::typeid_(Window{}) == hana::type<jbo::Window_Welch>{},
                    windowing::WelchWindow<
                        typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                        Complex>{},
                    windowing::NoWindow<
                        typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                        Complex>{}
                ))))))));
        }

        /** Creates a value of the selected normalization type at compilation time.
            \return value ... The selected value.
        */
        static constexpr auto getRadix2NormalizationValue(void)
        {
            return
                hana::if_(hana::typeid_(Normalization{}) == hana::type<jbo::Normalization_Division_By_Length>{},
                    normalization::DivisionByLengthNormalization<
                        typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                        typename decltype(std::integral_constant<int, 0>{})::type,
                        Complex>{},
                hana::if_(hana::typeid_(Normalization{}) == hana::type<jbo::Normalization_Square_Root>{},
                    normalization::SquareRootNormalization<
                        typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                        typename decltype(std::integral_constant<int, 0>{})::type,
                        Complex>{},
                    normalization::NoNormalization<
                        typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                        Complex>{}
                ));
         }

        /** Creates a tupel of sub task type values belonging to a radix 2 task at compilation time.
            \return hana::tuple_t ... A tuple of sub task type values.
        */
        static constexpr auto radix2SubTaskTypeValues(void)
        {
            return hana::if_(
                hana::typeid_(Decimation{}) == hana::type<jbo::Decimation_In_Time>{},
                    hana::tuple_t<
                        decltype(getWindowValue()),
                        basic::BitReversalIndexSwapping<
                            typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                            Complex>,
                        core::Radix2DIT<
                            typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                            typename decltype(getDirectionValue())::type,
                            Complex>,
                        decltype(getRadix2NormalizationValue())>,
                    hana::tuple_t<
                        decltype(getWindowValue()),
                        core::Radix2DIF<
                            typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                            typename decltype(getDirectionValue())::type,
                            Complex>,
                        basic::BitReversalIndexSwapping<
                            typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                            Complex>,
                        decltype(getRadix2NormalizationValue())>);
        }

        /** Creates a value of the selected normalization type at compilation time.
            \return value ... The selected value.
        */
        static constexpr auto getRadix4NormalizationValue(void)
        {
            return
                hana::if_(hana::typeid_(Normalization{}) == hana::type<jbo::Normalization_Division_By_Length>{},
                    normalization::DivisionByLengthNormalization<
                        typename decltype(std::integral_constant<int, 1 << (Stage::value << 1)>{})::type,
                        typename decltype(std::integral_constant<int, 0>{})::type,
                        Complex>{},
                hana::if_(hana::typeid_(Normalization{}) == hana::type<jbo::Normalization_Square_Root>{},
                    normalization::SquareRootNormalization<
                        typename decltype(std::integral_constant<int, 1 << (Stage::value << 1)>{})::type,
                        typename decltype(std::integral_constant<int, 0>{})::type,
                        Complex>{},
                    normalization::NoNormalization<
                        typename decltype(std::integral_constant<int, 1 << (Stage::value << 1)>{})::type,
                        Complex>{}
                ));
 
        }

        /** Creates a tupel of sub task type values belonging to a radix 4 task at compilation time.
            \return hana::tuple_t ... A tuple of sub task type values.
        */
        static constexpr auto radix4SubTaskTypeValues(void)
        {
            return hana::if_(
                hana::typeid_(Decimation{}) == hana::type<jbo::Decimation_In_Time>{},
                    hana::tuple_t<
                        decltype(getWindowValue()),
                        basic::BitReversalIndexSwapping<
                            typename decltype(std::integral_constant<int, 1 << (Stage::value << 1)>{})::type,
                            Complex>,
                        core::Radix4DIT<
                            typename decltype(std::integral_constant<int, 1 << (Stage::value << 1)>{})::type,
                            typename decltype(getDirectionValue())::type,
                            Complex>,
                        decltype(getRadix4NormalizationValue())>,
                    hana::tuple_t<
                        decltype(getWindowValue()),
                        core::Radix4DIF<
                            typename decltype(std::integral_constant<int, 1 << (Stage::value << 1)>{})::type,
                            typename decltype(getDirectionValue())::type,
                            Complex>,
                        basic::BitReversalIndexSwapping<
                            typename decltype(std::integral_constant<int, 1 << (Stage::value << 1)>{})::type,
                            Complex>,
                        decltype(getRadix4NormalizationValue())>);
        }

        /** Creates a tupel of sub task type values belonging to a split radix 2-4 task at compilation time.
            \return hana::tuple_t ... A tuple of sub task type values.
        */
        static constexpr auto radixSplit24SubtaskTypeValues(void)
        {
            return hana::if_(
                hana::typeid_(Decimation{}) == hana::type<jbo::Decimation_In_Time>{},
                    hana::tuple_t<
                        decltype(getWindowValue()),
                        basic::BitReversalIndexSwapping<
                            typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                            Complex>,
                        core::RadixSplit24DIT<
                            typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                            typename decltype(getDirectionValue())::type,
                            Complex>,
                        decltype(getRadix2NormalizationValue())>,
                    hana::tuple_t<
                        decltype(getWindowValue()),
                        core::RadixSplit24DIF<
                            typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                            typename decltype(getDirectionValue())::type,
                            Complex>,
                        basic::BitReversalIndexSwapping<
                            typename decltype(std::integral_constant<int, 1 << Stage::value>{})::type,
                            Complex>,
                        decltype(getRadix2NormalizationValue())>);
        }

        /** Creates a tupel of sub task type values at compilation time. Sub tasks belong to a FFT task.
            \return hana::tuple_t ... A tuple of sub task type values.
        */
        static constexpr auto createTupleOfSubTaskTypeValues(void)
        {
            if constexpr (hana::typeid_(Radix{}) == hana::type<jbo::Radix_2>{})
                return radix2SubTaskTypeValues();
            else if constexpr (hana::typeid_(Radix{}) == hana::type<jbo::Radix_4>{})
                return radix4SubTaskTypeValues();
            else
                return radixSplit24SubtaskTypeValues();
        }

        using SubTaskTypes = typename decltype(hana::unpack(createTupleOfSubTaskTypeValues(), hana::template_<hana::tuple>))::type;

        SubTaskTypes tupleOfSubTasks_;

    public:
        /** Executes all sub tasks sequentially.
            \param[in] data ... Pointer to an array of SampleCnt elements of type Complex.
        */
        void operator()(Complex* data) const override
        {
            hana::for_each(tupleOfSubTasks_, [&](const auto& subTask)
            {
                subTask(data);
            });
        }

        std::size_t numberOfSamples(void) const override
        {
            return getNumberOfSamples();
        }
        
        std::size_t numberOfFrequencies(void) const override
        {
            return getNumberOfFrequencies();
        }
    };
}