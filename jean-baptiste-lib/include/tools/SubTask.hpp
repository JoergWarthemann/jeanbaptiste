#ifndef JB_SUBTASK_HPP_
#define JB_SUBTASK_HPP_

#include <span>

namespace jb::tools {
/**
 * Uses CRTP to define a unique interface for compile time sub tasks.
 * \tparam Derived ... The derived class being used in compile time inheritance.
 * \tparam Complex ... Complex data type.
 */
template <typename Derived, typename Complex>
class SubTask {
public:
    template <typename SampleCnt>
    void operator()(std::span<Complex, SampleCnt::value> data) const
    {
        static_cast<Derived*>(this)->operator()(data);
    }
};

} // namespace jb::tools

#endif // JB_TOOLS_SUBTASK_HPP_
