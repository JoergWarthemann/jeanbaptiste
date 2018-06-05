#pragma once

#include "basic/BitReversalIndexSwapping.h"
#include "basic/Normalization.h"
#include <boost/hana.hpp>
#include <cassert>
#include "core/Radix2.h"
#include "core/Radix4.h"
#include "ExecutableAlgorithm.h"
#include "Options.h"

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
              typename Complex>
    class Algorithm
        : public ExecutableAlgorithm<Complex>
    {
        static constexpr auto getDirectionType(void)
        {
            return hana::if_(hana::decltype_(Direction{}) == hana::type_c<jbo::Direction_Forward>, hana::int_c<1>, hana::int_c<-1>);
        }

        static constexpr auto getWindowType(void)
        {

        }

        static constexpr auto radix2SubTaskTypes(void)
        {
            return hana::if_(
                hana::decltype_(Decimation{}) == hana::type_c<jbo::Decimation_In_Time>,
                    hana::tuple_t<
                        basic::BitReversalIndexSwapping<
                            typename decltype(hana::int_c<1 << Stage::value>)::type,
                            Complex>,
                        core::Radix2DIT<
                            typename decltype(hana::int_c<1 << Stage::value>)::type,
                            typename decltype(getDirectionType())::type,
                            Complex>,
                        basic::Normalization<
                            typename decltype(hana::int_c<1 << Stage::value>)::type,
                            typename decltype(hana::int_c<0>)::type,
                            Complex>>,
                    hana::tuple_t<
                        core::Radix2DIF<
                            typename decltype(hana::int_c<1 << Stage::value>)::type,
                            typename decltype(getDirectionType())::type,
                            Complex>,
                        basic::BitReversalIndexSwapping<
                            typename decltype(hana::int_c<1 << Stage::value>)::type,
                            Complex>,
                        basic::Normalization<
                            typename decltype(hana::int_c<1 << Stage::value>)::type,
                            typename decltype(hana::int_c<0>)::type,
                            Complex>>);
        }

        static constexpr auto radix4SubTaskTypes(void)
        {
            return hana::if_(
                hana::decltype_(Decimation{}) == hana::type_c<jbo::Decimation_In_Time>,
                    hana::tuple_t<
                        basic::BitReversalIndexSwapping<
                            typename decltype(hana::int_c<1 << (Stage::value << 1)>)::type,
                            Complex>,
                        core::Radix4DIT<
                            typename decltype(hana::int_c<1 << (Stage::value << 1)>)::type,
                            typename decltype(getDirectionType())::type,
                            Complex>,
                        basic::Normalization<
                            typename decltype(hana::int_c<1 << (Stage::value << 1)>)::type,
                            typename decltype(hana::int_c<0>)::type,
                            Complex>>,
                    hana::tuple_t<
                        core::Radix4DIF<
                            typename decltype(hana::int_c<1 << (Stage::value << 1)>)::type,
                            typename decltype(getDirectionType())::type,
                            Complex>,
                        basic::BitReversalIndexSwapping<
                            typename decltype(hana::int_c<1 << (Stage::value << 1)>)::type,
                            Complex>,
                        basic::Normalization<
                            typename decltype(hana::int_c<1 << (Stage::value << 1)>)::type,
                            typename decltype(hana::int_c<0>)::type,
                            Complex>>);
        }

        static constexpr auto splitRadixSubtaskTypes(void)
        {
            return hana::tuple_t<
                basic::BitReversalIndexSwapping<
                    std::integral_constant<typename Stage::value_type, 1 << Stage::value>,
                    Complex>
                >;
        }

        /** Creates a tupel of sub task types at compilation time. Sub tasks belong to a FFT task.
            \return hana::tuple_t ... A tuple of sub tasks.
        */
        static constexpr auto createTupleOfSubTaskTypes(void)
        {
            if constexpr (hana::decltype_(Radix{}) == hana::type_c<jbo::Radix_2>)
                return radix2SubTaskTypes();
            else if constexpr (hana::decltype_(Radix{}) == hana::type_c<jbo::Radix_4>)
                return radix4SubTaskTypes();
            else
                return splitRadixSubtaskTypes();
        }

        using SubTasks = typename decltype(hana::unpack(createTupleOfSubTaskTypes(), hana::template_<hana::tuple>))::type;

        SubTasks tupleOfSubTasks_;

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
    };
}