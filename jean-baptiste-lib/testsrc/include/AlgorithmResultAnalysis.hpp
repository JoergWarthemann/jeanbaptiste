#ifndef JB_TESTING_ALGORITHMRESULTANALYSIS_HPP_
#define JB_TESTING_ALGORITHMRESULTANALYSIS_HPP_

#include <complex>
#include <format>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <string_view>
#include <vector>

#include <boost/math/special_functions/next.hpp>

#include "TestCaseLoader.hpp"

namespace jb::testing {
class AlgorithmResultAnalysis {
public:
    /**
     * Loads the test case data.
     * \param testDataFile ... The path of a test data file.
     * \param tagA ... The identifier tag of the 1st data set.
     * \param dataSetA ... The 1st data set to be loaded.
     * \param tagB ... The identifier tag of the 2nd data set.
     * \param dataSetB ... The 2nd data set to be loaded.
     */
    template <typename TDataSetType>
    static bool initialize(const std::string_view testDataFile, const std::string_view tagA, TDataSetType& dataSetA,
        const std::string_view tagB, TDataSetType& dataSetB)
    {
        try {
            TestCaseLoader loader(testDataFile);
            return loader.getData(tagA, dataSetA, tagB, dataSetB);
        } catch (const std::exception& ex) {
            std::cerr << "Exception occurred when loading test data. " << ex.what() << std::endl;
        }

        return false;
    }

    /**
     * Compares the content of 2 data sets.
     * \param dataSetA ... The first data set to be compared.
     * \param dataSetB ... The second data set to be compared.
     */
    static void compareComplexDataSets(std::vector<std::complex<double>>& dataSetA,
        std::vector<std::complex<double>>& dataSetB)
    {
        ASSERT_EQ(dataSetA.size(), dataSetB.size()) << "The lengths of both vectors need to be equal.";

        for (auto i = 0; i < dataSetA.size(); ++i) {
            EXPECT_NEAR(dataSetA[i].real(), dataSetB[i].real(), kPrecision)
                << std::format("Real part mismatch at position {}", i);
            EXPECT_NEAR(dataSetA[i].imag(), dataSetB[i].imag(), kPrecision)
                << std::format("Imaginary part mismatch at position {}", i);
        }
    }

private:
    static constexpr double kPrecision{0.0001};
};

} // namespace jb::testing

#endif // JB_TESTING_ALGORITHMRESULTANALYSIS_HPP_