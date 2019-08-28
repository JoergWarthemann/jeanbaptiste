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

        virtual std::size_t numberOfSamples(void) const
        {
            return 0;
        }
        
        virtual std::size_t numberOfFrequencies(void) const
        {
            return 0;
        }
    };
}