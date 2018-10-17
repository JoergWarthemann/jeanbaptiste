#pragma once

namespace jeanbaptiste::basic
{
    template <typename T = double>
    constexpr std::decay_t<T> ceil(const T& value)
    {
        return (static_cast<float>(static_cast<int32_t>(value)) == value)
            ? value
            : static_cast<std::decay_t<T>>(static_cast<int32_t>(value) + ((value > T{}) ? 1 : T{}));
    }
}