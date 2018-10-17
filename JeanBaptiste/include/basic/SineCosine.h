#pragma once

#include <boost/math/constants/constants.hpp>
#include "Round.h"
#include <type_traits>

/** std::decay returns the underlying type by removing reference and const.
    Use the shorthand std::decay_t<const X&> instead of std::decay<const X&>::type.
    e.g.: std::decay_t<const X&>    returns X
          std::decay_t<const X>     returns X
          std::decay_t<X&>          returns X
*/

namespace constants = boost::math::constants;

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
            : 1.0 - x * x / seriesStart / (seriesStart + 1) * sineCosineSeries(seriesStart + 2, seriesEnd, x);
    }

    /** Executes Cody-Waite range reduction on x.
    */
    template <typename T = double>
    constexpr std::decay_t<T> reduceRange(T x, int32_t& quadrant)
    {
        // Find out in which quadrant of the unit circle x is situated.
        quadrant = static_cast<int32_t>(round<T>(x * constants::two_div_pi<T>()));
        auto temp = x - quadrant * constants::half_pi<T>();
        //quadrant = (quadrant + 1) >> 1;
        
        return (((quadrant & 1) == 0) ? T{} : constants::pi<T>()) - temp;
    }

    ///////////////////////////////////////////

    template <typename T = double>
    constexpr std::decay_t<T> testSine(const T x)
    {
        static_assert(std::is_floating_point<T>::value, "Trying to generate sine using a non floating point type.");

        int32_t quadrant{0};
        auto reducedX = reduceRange<T>(x, quadrant);

        // Sine starts with the 2nd element of the power series: x^1 = x.
        // The number of terms affects the precision of the calculated sin value.
        // It differs when dealing with 4 or 8 byte datatypes respecting rounding errors.
        return reducedX * sineCosineSeries<T>(2, (sizeof(T) > 4) ? 34 : 24, reducedX)/* * (((quadrant & 2) == 2) ? -1.0 : 1.0)*/;
    }

    template <typename T = double>
    constexpr std::decay_t<T> testCosine(const T x)
    {
        static_assert(std::is_floating_point<T>::value, "Trying to generate cosine using a non floating point type.");

        int32_t quadrant{0};
        auto reducedX = reduceRange<T>(x, quadrant);

        // Cosine starts with the 1st element of the power series : x^0 = 1.
        // The number of terms affects the precision of the calculated cos value.
        // It differs when dealing with 4 or 8 byte datatypes respecting rounding errors.
        return sineCosineSeries<T>(1, (sizeof(T) > 4) ? 33 : 23, reducedX)/* * (((quadrant & 2) == 2) ? -1.0 : 1.0)*/;
    }

    ///////////////////////////////////////////

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
        // It differs when dealing with 4 or 8 byte datatypes respecting rounding errors.
        return x * sineCosineSeries<T>(2, (sizeof(T) > 4) ? 76 : 66, x);
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
        // It differs when dealing with 4 or 8 byte datatypes respecting rounding errors.
        return sineCosineSeries<T>(1, (sizeof(T) > 4) ? 75 : 65, x);
    }
}