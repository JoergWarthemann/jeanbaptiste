#ifndef JEANBAPTISTE_TESTING_ALGORITHMRESULTANALYSIS_HPP_
#define JEANBAPTISTE_TESTING_ALGORITHMRESULTANALYSIS_HPP_

#include <boost/math/special_functions/next.hpp>
//#include <boost/test/unit_test.hpp>
#include <complex>
#include <format>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <iostream>
#include <string_view>
#include <vector>

#include "TestCaseLoader.hpp"

//namespace tt = boost::test_tools;

namespace jeanbaptiste::testing {

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
    template<typename TDataSetType>
    static bool initialize(const std::string_view testDataFile, const std::string_view tagA, TDataSetType& dataSetA,
        const std::string_view tagB, TDataSetType& dataSetB)
    {
        std::cout << "Initialize data." << std::endl;

        try {
            TestCaseLoader loader(testDataFile);
            return loader.getData(tagA, dataSetA, tagB, dataSetB);
        }
        catch (std::exception ex) {
            std::cout << "Exception occurred when loading test data. " << ex.what() << std::endl;
        }

        return false;
    }

    // /** Loads the test case data.
    //     \param[in] testDataFile ... The path of a test data file.
    //     \param[in] identifier1 ... The 1st data tag identifier.
    //     \param[out] workingSet ... The set of working (input) data.
    //     \param[out] expected1 ... The 1st set of expected (output) data.
    //     \param[in] identifier2 ... The 2nd data tag identifier.
    //     \param[out] expected2 ... The 2nd set of expected (output) data.
    // */
    // bool initialize(const std::string_view testDataFile, const std::string_view tagA, TComplexDataSetType& dataSetA,
    //     const std::string_view tagB, TComplexDataSetType& dataSetB)
    // {
    //     //BOOST_TEST_MESSAGE("Initialize data.");
    //     std::cout << "Initialize data." << std::endl;

    //     try
    //     {
    //         TestCaseLoader<T> loader(testDataFile);
    //         return loader.getData(identifier1, workingSet, expected1, identifier2, expected2);
    //     }
    //     catch (std::exception ex)
    //     {
    //         //BOOST_TEST(false, (boost::format("Exception caught loading test data. %s") % ex.what()).str());
    //         std::cout << "Exception occurred when loading test data. " << ex.what() << std::endl;
    //     }

    //     return false;
    // }

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

            // bool realMatch = std::abs(dataSetA[i].real() - dataSetB[i].real()) <= static_cast<T>(kPrecision);
            // bool imagMatch = std::abs(dataSetA[i].imag() - dataSetB[i].imag()) <= static_cast<T>(kPrecision);
            // EXPECT_TRUE(realMatch && imagMatch)
            //     << std::format("Mismatch at position {}: ({}, {}i) != ({}, {}i)",
            //         i,
            //         dataSetA[i].real(),
            //         dataSetA[i].imag(),
            //         dataSetB[i].real(),
            //         dataSetB[i].imag());
        }
    }

private:
    static constexpr double kPrecision{0.0001};
};

}

#endif // JEANBAPTISTE_TESTING_ALGORITHMRESULTANALYSIS_HPP_