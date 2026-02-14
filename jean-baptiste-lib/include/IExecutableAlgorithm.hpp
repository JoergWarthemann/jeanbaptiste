#ifndef JEANBAPTISTE_EXECUTABLE_ALGORITHM_HPP_
#define JEANBAPTISTE_EXECUTABLE_ALGORITHM_HPP_

#include <complex>

namespace jeanbaptiste
{
   /**
     * Defines a unique interface for dynamic algorithms.
     * @param Complex Complex data type.
     */
    template <typename Complex = std::complex<double>>
    class IExecutableAlgorithm
    {
    public:
        virtual ~IExecutableAlgorithm()
        {}

        virtual void operator()(Complex* data) const
        {}
    };
}

#endif // JEANBAPTISTE_EXECUTABLE_ALGORITHM_HPP_