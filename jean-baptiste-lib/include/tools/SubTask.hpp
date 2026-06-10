#ifndef JB_SUBTASK_HPP_
#define JB_SUBTASK_HPP_

#include <span>

namespace jb::tools {
/**
 * Uses CRTP to define a unique interface for compile time sub tasks.
 * \param Derived ... The derived class being used in compile time inheritance.
 * \param Complex ... Complex data type.
 */
template <typename Derived, typename Complex>
class SubTask {
public:
    // void operator()(Complex* data) const
    //{
    //    static_cast<Derived*>(this)->operator()(data);
    //}
    template <typename SampleCnt>
    void operator()(std::span<Complex, SampleCnt::value> data) const
    {
        static_cast<Derived*>(this)->operator()(data);
    }
};

} // namespace jb::tools

#endif // JB_TOOLS_SUBTASK_HPP_
