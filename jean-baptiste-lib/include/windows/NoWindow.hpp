#ifndef JB_WINDOWS_NOWINDOW_HPP_
#define JB_WINDOWS_NOWINDOW_HPP_

#include "tools/SubTask.hpp"

namespace jb::windows {

/**
 * Creates an empty window (rectangular) for a specified sample count.
 * \param SampleCnt ... The count of samples to be processed in this recursion level (stage)
 * \param Complex ... The complex type.
 */
template <typename SampleCnt, typename Complex>
    requires std::is_integral_v<typename SampleCnt::value_type> &&
    std::is_floating_point_v<typename Complex::value_type>
class NoWindow
    : public tools::SubTask<NoWindow<SampleCnt, Complex>, Complex> {
public:
    void operator()(std::span<Complex, SampleCnt::value> data) const
    {}
};

} // namespace jb::windows

#endif // JB_WINDOWS_NOWINDOW_HPP_