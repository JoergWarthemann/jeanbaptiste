#ifndef JB_TESTING_ALGORITHM_FIXTURE_HPP_
#define JB_TESTING_ALGORITHM_FIXTURE_HPP_

#include <complex>
#include <memory>
#include <vector>

#include "AlgorithmResultAnalysis.hpp"
#include "IExecutableAlgorithm.hpp"

namespace jb::testing {

/** Fixture for executable algorithms.
 *   Provides common data and functions for algorithm tests.
 */
class AlgorithmFixture {
protected:
    using TAlgorithmType = std::unique_ptr<IExecutableAlgorithm<std::complex<double>>>;
    using TDataSetType = std::vector<double>;
    using TComplexDataSetType = std::vector<std::complex<double>>;

    TComplexDataSetType mDataSetA;
    TComplexDataSetType mDataSetB;

public:
    AlgorithmFixture(void) = default;
    virtual ~AlgorithmFixture(void) = default;

    void verifyAlgorithm(TAlgorithmType algorithm)
    {
        (*algorithm)(mDataSetA);
        AlgorithmResultAnalysis::compareComplexDataSets(mDataSetA, mDataSetB);
    }

    template <typename... Algorithm>
    void verifyAlgorithms(Algorithm... algorithm)
    {
        (verifyAlgorithm(algorithm), ...);
    }
};

} // namespace jb::testing

#endif // JB_TESTING_ALGORITHM_FIXTURE_HPP_
