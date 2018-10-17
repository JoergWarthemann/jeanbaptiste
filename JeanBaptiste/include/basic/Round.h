#pragma once

#include "Abs.h"
#include "Floor.h"

namespace jeanbaptiste::basic
{
    template <typename T = double>
    constexpr std::decay_t<T> round(const T& value)
    {
        auto temp = floor<T>(abs<T>(value) + 0.5);
        return (T{} < value) ? temp : -temp;
    }
}