#ifndef JEANBAPTISTE_IEXECUTABLE_ALGORITHM_HPP_
#define JEANBAPTISTE_IEXECUTABLE_ALGORITHM_HPP_

#include <span>

namespace jeanbaptiste
{
   /**
     * Defines a unique interface for dynamic algorithms.
     * @param Complex Complex data type.
     */
    template <typename Complex>
        requires std::is_floating_point_v<typename Complex::value_type>
    class IExecutableAlgorithm
    {
    public:
        virtual ~IExecutableAlgorithm() = default;
        virtual void operator()(std::span<Complex> data) const = 0;
    };
}

#endif // JEANBAPTISTE_IEXECUTABLE_ALGORITHM_HPP_