/**
 * Compile-time sample-count selection for FFT algorithm variants.
 *
 * Defines the external sample count expected by the public algorithm interface
 * and the internal sample count processed by the core FFT implementation for a
 * given stage, radix, and data type.
 *
 * Complex FFTs use the same internal and external sample count. Real radix-2
 * and split-radix FFTs expose a 2N-sample real layout and process N complex
 * samples internally. Real radix-4 FFTs use the radix-4 complex sample count
 * internally and expose twice that count externally.
 *
 *
 *    ------------------------------------------------------------------------------------------
 *    | stages (id) | sample count          | sample count             | sample count          |
 *    |             | (internal = external) | (internal = external)    | (internal = external) |
 *    |             | radix 2               | radix 4                  | split radix           |
 *    ------------------------------------------------------------------------------------------
 *    | 1           | 2^1 =  2              | 4^1 =   4                | 2^1 =  2              |
 *    | 2           | 2^2 =  4              | 4^2 =  16                | 2^2 =  4              |
 *    | 3           | 2^3 =  8              | 4^3 =  64                | 2^3 =  8              |
 *    | 4           | 2^4 = 16              | 4^4 = 256                | 2^4 = 16              |
 *    | ...         | ...                   | ...                      | ...                   |
 *    ------------------------------------------------------------------------------------------
 *
 *
 *    ------------------------------------------------------------------------------------------
 *    | stages (id) | sample count          | sample count             | sample count          |
 *    |             | (internal, external)  | (internal, external)     | (internal, external)  |
 *    |             | real radix 2          | real radix 4             | real split radix      |
 *    ------------------------------------------------------------------------------------------
 *    | 1           | -                     | 4^1 =   4, 4^1 * 2 = 8   | -                     |
 *    | 2           | 2^2-1 = 2, 2^2 =  4   | 4^2 =  16, 4^2 * 2 = 32  | 2^2-1 = 2, 2^2 =  4   |
 *    | 3           | 2^3-1 = 4, 2^3 =  8   | 4^3 =  64, 4^3 * 2 = 128 | 2^3-1 = 4, 2^3 =  8   |
 *    | 4           | 2^4-1 = 8, 2^4 = 16   | 4^4 = 256, 4^4 * 2 = 512 | 2^4-1 = 8, 2^4 = 16   |
 *    | ...         | ...                   | ...                      | ...                   |
 *    ------------------------------------------------------------------------------------------
 */

#ifndef JB_SAMPLECOUNT_HPP_
#define JB_SAMPLECOUNT_HPP_

#include <cstddef>
#include <type_traits>

#include "Options.hpp"

namespace jb {

template <std::size_t Stage, typename Radix, typename Data>
struct SampleCount;

template <std::size_t Stage>
struct SampleCount<Stage, options::Radix_2, options::Data_Complex> {
    using External = std::integral_constant<int, 1 << Stage>;
    using Internal = External;
};

template <std::size_t Stage>
struct SampleCount<Stage, options::Radix_Split_2_4, options::Data_Complex> {
    using External = std::integral_constant<int, 1 << Stage>;
    using Internal = External;
};

template <std::size_t Stage>
struct SampleCount<Stage, options::Radix_4, options::Data_Complex> {
    using External = std::integral_constant<int, 1 << (Stage << 1)>;
    using Internal = External;
};

template <std::size_t Stage>
struct SampleCount<Stage, options::Radix_2, options::Data_Real> {
    using External = std::integral_constant<int, 1 << Stage>;
    using Internal = std::integral_constant<int, External::value / 2>;
};

template <std::size_t Stage>
struct SampleCount<Stage, options::Radix_Split_2_4, options::Data_Real> {
    using External = std::integral_constant<int, 1 << Stage>;
    using Internal = std::integral_constant<int, External::value / 2>;
};

template <std::size_t Stage>
struct SampleCount<Stage, options::Radix_4, options::Data_Real> {
    using Internal = std::integral_constant<int, 1 << (Stage << 1)>;
    using External = std::integral_constant<int, Internal::value * 2>;
};

} // namespace jb

#endif // JB_SAMPLECOUNT_HPP_