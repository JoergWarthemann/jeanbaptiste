#pragma once

namespace jeanbaptiste::windowing
{
    template <typename Complex>
    class WindowExecutionFunctor
    {
        using ValueType = typename Complex::value_type;

    public:
        Complex operator()(const Complex& factor1, const ValueType& factor2)
        {
            return Complex(factor1.real() * factor2, factor1.imag());
        }
    };
}