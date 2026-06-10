#ifndef JB_WINDOWS_EXECUTEWINDOWONCOMPLEXDATA_HPP_
#define JB_WINDOWS_EXECUTEWINDOWONCOMPLEXDATA_HPP_

namespace jb::windows {

template <typename Complex>
class ExecuteWindowOnComplexData {
public:
    Complex operator()(const Complex& factor1, const Complex::value_type& factor2)
    {
        return Complex(factor1.real() * factor2, factor1.imag());
    }
};

} // namespace jb::windows

#endif // JB_WINDOWS_EXECUTEWINDOWONCOMPLEXDATA_HPP_