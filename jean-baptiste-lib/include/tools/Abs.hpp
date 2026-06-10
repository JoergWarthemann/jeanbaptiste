#ifndef JB_TOOLS_ABS_HPP_
#define JB_TOOLS_ABS_HPP_

#include <type_traits>
namespace jb::tools {

template <typename T = double>
constexpr std::decay_t<T> abs(const T& value)
{
    return (T{} < value) ? value : -value;
}

} // namespace jb::tools

#endif // JB_TOOLS_ABS_HPP_