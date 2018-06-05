#pragma once

namespace jeanbaptiste
{
    /** Uses CRTP to define a unique interface for compile time sub tasks.
        \param Derived ... The derived class being used in compile time inheritance.
        \param Complex ... Complex data type.
    */
    template <typename Derived,
              typename Complex>
    class SubTask
    {
    public:
        void operator()(Complex* data) const
        {
            static_cast<Derived*>(this)->operator()(data);
        }
    };
}