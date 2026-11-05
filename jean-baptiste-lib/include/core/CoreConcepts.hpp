#ifndef JB_CORE_CONCEPTS_HPP_
#define JB_CORE_CONCEPTS_HPP_

#include <type_traits>

namespace jb::core {

template <typename T>
concept IntegralConstant = std::is_integral_v<typename T::value_type>;

template <typename T>
concept ComplexLike = std::is_floating_point_v<typename T::value_type>;

template <typename SampleCnt, typename DirectionFactor, typename Complex>
concept FftParameters = IntegralConstant<SampleCnt> &&
    IntegralConstant<DirectionFactor> &&
    ComplexLike<Complex>;

} // namespace jb::core

#endif // JB_CORE_CONCEPTS_HPP_