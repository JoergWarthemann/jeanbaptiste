#ifndef JEANBAPTISTE_WINDOWING_EXECUTEWINDOWONCOMPLEXDATA_HPP_
#define JEANBAPTISTE_WINDOWING_EXECUTEWINDOWONCOMPLEXDATA_HPP_

namespace jeanbaptiste::windowing {

template <typename Complex>
class ExecuteWindowOnComplexData {
public:
    Complex operator()(const Complex& factor1, const Complex::value_type& factor2)
    {
        return Complex(factor1.real() * factor2, factor1.imag());
    }
};

} // namespace jeanbaptiste::windowing

#endif // JEANBAPTISTE_WINDOWING_EXECUTEWINDOWONCOMPLEXDATA_HPP_