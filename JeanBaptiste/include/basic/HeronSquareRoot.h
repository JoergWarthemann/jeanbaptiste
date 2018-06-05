#pragma once

namespace jeanbaptiste::basic
{
    /** Recursively calculates an estimation of the square root of an integer.
        \param[in] SeriesStart ... The point to start the recursion.
        \param[in] SeriesEnd ... The end point of the recursion.
        \param[in] Radicant ... The radicand that is to be square rooted.
        \param[in] guess ... The current approach to the actual resulting value.
    */
    template <typename T = double>
    constexpr std::decay_t<T> squareRoot(const std::size_t seriesStart, const std::size_t seriesEnd, const std::size_t radicant, const double guess)
    {
        static_assert(std::is_floating_point<T>::value ||std::is_integral<T>::value, "Trying to calculate the square root on a non integral or floating point type.");

        // Special cases that stop the recursion:
        // - square root of 1 is 1
        // - square root of 0 is 0
        // - start and end point of recursion are equal

        return
            (radicant == 0)
            ? 0.0
            :
                (radicant == 1)
                ? 1.0
                :
                    (seriesStart == seriesEnd)
                    ? (guess + radicant / guess) / 2.0
                    : squareRoot(seriesStart, seriesEnd - 1, radicant, (guess + radicant / guess) / 2.0);
    }

    template <typename T = double>
    constexpr std::decay_t<T> squareRoot(const std::size_t seriesStart, const std::size_t seriesEnd, const std::size_t radicant)
    {
        return squareRoot(seriesStart, seriesEnd, radicant, radicant / 2.0);
    }
}