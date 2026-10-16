#ifndef JB_IEXECUTABLE_ALGORITHM_HPP_
#define JB_IEXECUTABLE_ALGORITHM_HPP_

#include <span>

namespace jb {
/**
 * Defines a unique interface for dynamic algorithms.
 * \tparam Complex Complex data type.
 */
template <typename Complex>
    requires std::is_floating_point_v<typename Complex::value_type>
class IExecutableAlgorithm {
public:
    virtual ~IExecutableAlgorithm() = default;
    virtual void operator()(std::span<Complex> data) const = 0;
};
} // namespace jb

#endif // JB_IEXECUTABLE_ALGORITHM_HPP_