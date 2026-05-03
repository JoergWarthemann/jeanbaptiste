#ifndef JEANBAPTISTE_TOOLS_ABS_HPP_
#define JEANBAPTISTE_TOOLS_ABS_HPP_

#include <type_traits>
namespace jeanbaptiste::tools {

template <typename T = double>
constexpr std::decay_t<T> abs(const T& value)
{
    return (T{} < value) ? value : -value;
}

}

#endif // JEANBAPTISTE_TOOLS_ABS_HPP_