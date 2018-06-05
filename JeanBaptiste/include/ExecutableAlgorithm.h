#pragma once

namespace jeanbaptiste
{
    /** Defines a unique interface for dynamic algorithms.
        \param Complex ... Complex data type.
    */
    template <typename Complex = std::complex<double>>
    class ExecutableAlgorithm
    {
    public:
        virtual ~ExecutableAlgorithm()
        {}

        virtual void operator()(Complex* data) const
        {}
    };
}