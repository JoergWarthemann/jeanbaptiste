#pragma once

#include <boost/math/constants/constants.hpp>
#include <type_traits>

namespace constants = boost::math::constants;

/** std::decay returns the underlying type by removing reference and const.
    Use the shorthand std::decay_t<const X&> instead of std::decay<const X&>::type.
    e.g.: std::decay_t<const X&>    returns X
          std::decay_t<const X>     returns X
          std::decay_t<X&>          returns X
*/

namespace jeanbaptiste::basic
{
    /** Recursively calculates a single term of a Horner schematized power series for sin/cos.
        Causes the recursive calculation of previous terms until (seriesStart >= seriesEnd) -
        thats where the break condition in the constexpr recursion matches.
        \param[in] seriesStart ... Calculate series with terms starting with seriesStart.
        \param[in] seriesEnd ... Calculate series with terms until seriesEnd.
        \param[in] x ... The value whose sin/cos is to be calculated.
    */
    template <typename T = double>
    constexpr std::decay_t<T> sineCosineSeries(const std::size_t seriesStart, const std::size_t seriesEnd, const T x)
    {
        return
            (seriesStart >= seriesEnd)
            ? 1.0
            : 1 - x * x / seriesStart / (seriesStart + 1) * sineCosineSeries(seriesStart + 2, seriesEnd, x);
    }

    /** Calculation of sine by a Horner schematized power series:
        sin x = x (1 - x^2 (1/3! - x^2/5! + x^4/7! - x^6/9! ...))
              = x (1 - x^2 (1/3! - x^2 (1/5! - x^2 (1/7! - x^2 (1/9! ...)))))
       \param[in] x ... The value whose sine is to be calculated.
    */
    template <typename T = double>
    constexpr std::decay_t<T> sine(const T x)
    {
        static_assert(std::is_floating_point<T>::value, "Trying to generate sine using a non floating point type.");

        // Sine starts with the 2nd element of the power series: x^1 = x.
        // The number of terms affects the precision of the calculated sin value.
        // 34 terms is enough for 8 byte datatypes like double.
        // 24 terms is enough for 4 byte datatypes like float.
        return x * sineCosineSeries<T>(2, (sizeof(T) > 4) ? 34 : 24, x);
    }

    /* Calculation of cosine by a Horner schematized power series:
       cos x = 1 - x^2 (1/2! - x^2/4! + x^4/6! - x^6/8! ...)
             = 1 - x^2 (1/2! - x^2 (1/4! - x^2 (1/6! - x^2 (1/8! ...))))
       \param[in] x ... The value whose sine is to be calculated.
    */
    template <typename T = double>
    constexpr std::decay_t<T> cosine(const T x)
    {
        static_assert(std::is_floating_point<T>::value, "Trying to generate cosine using a non floating point type.");

        // Cosine starts with the 1st element of the power series : x^0 = 1.
        // The number of terms affects the precision of the calculated cos value.
        // 33 terms is enough for 8 byte datatypes like double.
        // 23 terms is enough for 4 byte datatypes like float.
        return sineCosineSeries<T>(1, (sizeof(T) > 4) ? 33 : 23, x);
    }
}