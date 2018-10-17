#pragma once

namespace jeanbaptiste::basic
{
    template <typename T = double>
    constexpr std::decay_t<T> floor(const T& value)
    {
        return (static_cast<float>(static_cast<int32_t>(value)) == value)
            ? value
            : static_cast<std::decay_t<T>>( static_cast<int32_t>(value) - ((T{} < value) ? T{} : 1));
    }
}