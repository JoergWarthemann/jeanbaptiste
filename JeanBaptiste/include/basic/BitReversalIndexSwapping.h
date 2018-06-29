#pragma once

#include <memory>
#include <boost/hana.hpp>
#include <complex>
#include "../SubTask.h"
#include <type_traits>

#include <iostream>

namespace hana = boost::hana;

namespace jeanbaptiste::basic
{
    /** Applies bit reversal index swapping which brings calculations into the right order before/after having applied a radix 2 FFT.
        \param SampleCnt ... The count of samples used for the bit reversal.
                             This masks the significant sequence of bits in a 16 bit number.
        \param Complex ... The complex data type.
    */
    template<typename SampleCnt,
             typename Complex>
    class BitReversalIndexSwapping
        : public SubTask<BitReversalIndexSwapping<SampleCnt, Complex>,
                         Complex>
    {
        /** The reversal is based on 16 bits always. There is a significant sequence of bits within these 16 bits. For example having a 
            field of 256 numbers this significant sequence of bits is about 8 bits (2^8 = 256). The result of the reversal is always 
            based on 16 bits. In order to have it based on 8 bits again, we have to shift it 8 bits to the right.

            This creates a lookup table of 2^n cases at compilation time. A single lookup masks the significant bits of a 16 bit number.
            example 1, given: sample count: 256 -> 1 byte
                              x =  19 =      0001 0011
                              wanted:   the corresponding shift value to represent the reversal of x in the range [0 ... 256] as a 8 bit number
                              solution: 200     =      1100 1000
            example 2, given: sample count: 1024 -> 2 bytes
                              x = 614 = 0010 0110 0110
                              wanted:   the corresponding shift value to represent the reversal of x in the range [0 ... 1024] as a 12 bit number
                              solution: 1636    = 0110 0110 0100
        */
        static constexpr auto createSignificantBitShiftLookupTable(void)
        {
            auto shifts = hana::make_range(hana::int_c<1>, hana::int_c<sizeof(char16_t) * 8 + 1>);

            return hana::unpack(shifts, [](auto... shift)
            {
                return hana::make_map(
                    hana::make_pair(
                        hana::int_c<1> << (hana::int_c<sizeof(char16_t) * 8> - shift),
                        shift)...);
            });
        }

        /** Calculates the bit reversed counterpart of a given index.
            \param[in] index ... The input number.
            \return std::size_t ... The bit reversed counterpart of input.
        */
        static constexpr std::size_t getReverseIndex(const std::size_t index)
        {
            std::size_t value1 = (( index & 0xAAAA) >> 1) | (( index & 0x5555) << 1);
            std::size_t value2 = ((value1 & 0xCCCC) >> 2) | ((value1 & 0x3333) << 2);
            std::size_t value3 = ((value2 & 0xF0F0) >> 4) | ((value2 & 0x0F0F) << 4);

            return               ((value3 & 0xFF00) >> 8) | ((value3 & 0x00FF) << 8)
                              >> significantBitShiftLookupTable_[hana::int_c<SampleCnt::value>];
        }

        static constexpr std::size_t createReverseIndicesPair(const std::size_t index)
        {
            return getReverseIndex(index);
        }

        template<std::size_t... Indices>
        static constexpr auto createReverseIndicesLookupTable(std::index_sequence<Indices...>)
        {
            return std::array<std::size_t, sizeof...(Indices)>
            {
                createReverseIndicesPair(Indices)...
            };
        }

        /** Creates a lookup table of indices and their bit reversed counterparts at compilation time.
        */
        static constexpr auto getSwapIndicesLookupTable(void)
        {
            static_assert(sizeof(SampleCnt) <= 2, "Trying to use a sample count bigger than 2 bytes.");

            return createReverseIndicesLookupTable(std::make_index_sequence<SampleCnt::value>{});
        }

        static constexpr auto significantBitShiftLookupTable_ = createSignificantBitShiftLookupTable();
        static constexpr auto swapLookupTable_ = getSwapIndicesLookupTable();

    public:
        /** Swaps all elements of data based on their index and its bit reversed counterpart index.
            \param[in] data ... Pointer to an array of SampleCnt elements of type Complex.
        */
        void operator()(Complex* data) const
        {
            for (auto i = 0; i < SampleCnt::value; ++i)
            {
                if (i < swapLookupTable_[i])
                    std::swap(data[i], data[swapLookupTable_[i]]);
            }
        }
    };
}