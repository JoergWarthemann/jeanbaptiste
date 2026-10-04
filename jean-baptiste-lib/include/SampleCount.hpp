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