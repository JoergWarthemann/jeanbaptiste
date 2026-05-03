#ifndef JEANBAPTISTE_WINDOWING_NOWINDOW_HPP_
#define JEANBAPTISTE_WINDOWING_NOWINDOW_HPP_

#include "tools/SubTask.hpp"

namespace jeanbaptiste::windowing {

/**
 * Creates an empty window (rectangular) for a specified sample count.
 * \param SampleCnt ... The count of samples to be processed in this recursion level (stage)
 * \param Complex ... The complex type.
 */
template <typename SampleCnt, typename Complex>
class NoWindow
    : public SubTask<NoWindow<SampleCnt, Complex>, Complex> {
public:
    void operator()(std::span<Complex, SampleCnt::value> data) const
    {}
};

} // namespace jeanbaptiste::windowing

#endif // JEANBAPTISTE_WINDOWING_NOWINDOW_HPP_