#pragma once

namespace jeanbaptiste::basic
{
    template <typename T = double>
    constexpr std::decay_t<T> abs(const T& value)
    {
        return (T{} < value) ? value : -value;
    }
}