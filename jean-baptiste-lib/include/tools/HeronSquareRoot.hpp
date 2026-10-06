#ifndef JB_TOOLS_HERONSQUAREROOT_HPP_
#define JB_TOOLS_HERONSQUAREROOT_HPP_

#include <type_traits>

namespace jb::tools {

/**
 * Recursively calculates an estimation of the square root of an integer.
 * \tparam T ... The type of the resulting square root (default is double).
 * \param seriesStart ... The point to start the recursion.
 * \param seriesEnd ... The end point of the recursion.
 * \param radicant ... The radicand that is to be square rooted.
 * \param guess ... The current approach to the actual resulting value.
 */
template <typename T = double>
    requires std::is_floating_point_v<T> ||
    std::is_integral_v<T>
constexpr std::decay_t<T> squareRoot(std::size_t seriesStart, std::size_t seriesEnd, std::size_t radicant, double guess)
{
    // Special cases that stop the recursion:
    // - square root of 1 is 1
    // - square root of 0 is 0
    // - start and end point of recursion are equal

    return (radicant == 0)
        ? 0.0
        : (radicant == 1)
        ? 1.0
        : (seriesStart == seriesEnd)
        ? (guess + radicant / guess) / 2.0
        : squareRoot(seriesStart, seriesEnd - 1, radicant, (guess + radicant / guess) / 2.0);
}

template <typename T = double>
    requires std::is_floating_point_v<T> || std::is_integral_v<T>
constexpr std::decay_t<T> squareRoot(std::size_t seriesStart, std::size_t seriesEnd, std::size_t radicant)
{
    return squareRoot(seriesStart, seriesEnd, radicant, radicant / 2.0);
}

} // namespace jb::tools

#endif // JB_TOOLS_HERONSQUAREROOT_HPP_