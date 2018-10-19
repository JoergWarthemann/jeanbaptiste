#pragma once

#include <boost/math/constants/constants.hpp>
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
    namespace internal
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

        /** Executes Cody-Waite range reduction on x to keep it in [-pi, pi].
        */
        template <typename T = double>
        constexpr std::decay_t<T> reduceRange(T x)
        {
            // We reduce x into [-pi, pi], so c=2*pi and c=c1+c2.
            constexpr T c1{6.283203125};
            constexpr T c2{-1.7817819752963e-5};

            // How many times is 2*pi in x?
            int32_t nearestInteger = static_cast<int32_t>(x * constants::one_div_two_pi<T>());

            // Return x within the reduced range.
            return (x - nearestInteger * c1) - nearestInteger * c2;
        }
    }

    template <typename T = double>
    constexpr std::decay_t<T> sine(const T x)
    {
        static_assert(std::is_floating_point<T>::value, "Trying to generate sine using a non floating point type.");

        auto reducedX = internal::reduceRange<T>(x);

        // Sine starts with the 2nd element of the power series: x^1 = x.
        // The number of terms affects the precision of the calculated sin value.
        // It differs when dealing with 4 or 8 byte datatypes respecting rounding errors.
        return reducedX * internal::sineCosineSeries<T>(2, (sizeof(T) > 4) ? 34 : 24, reducedX);
    }

    template <typename T = double>
    constexpr std::decay_t<T> cosine(const T x)
    {
        static_assert(std::is_floating_point<T>::value, "Trying to generate sine using a non floating point type.");

        auto reducedX = internal::reduceRange<T>(x);

        // Cosine starts with the 1st element of the power series : x^0 = 1.
        // The number of terms affects the precision of the calculated cos value.
        // It differs when dealing with 4 or 8 byte datatypes respecting rounding errors.
        return internal::sineCosineSeries<T>(1, (sizeof(T) > 4) ? 33 : 23, reducedX);
    }
}