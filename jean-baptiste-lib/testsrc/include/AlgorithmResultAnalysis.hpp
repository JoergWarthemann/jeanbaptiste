#ifndef JEANBAPTISTE_UTILITIES_ALGORITHMRESULTANALYSIS_HPP_
#define JEANBAPTISTE_UTILITIES_ALGORITHMRESULTANALYSIS_HPP_

#include <boost/math/special_functions/next.hpp>
//#include <boost/test/unit_test.hpp>
#include <complex>
#include <format>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <vector>

#include "TestCaseLoader.hpp"

//namespace tt = boost::test_tools;

namespace Utilities
{

template<typename T>
class AlgorithmResultAnalysis
{
    const double kPrecision = 0.0001;

public:
    /** Loads the test case data.
        \param[in] file ... The file path.
        \param[in] identifier1 ... The 1st data tag identifier.
        \param[out] workingSet ... The set of working (input) data.
        \param[in] identifier2 ... The 2nd data tag identifier.
        \param[out] expected ... The set of expected (output) data.
    */
    bool initialize(const std::string& file, std::string identifier1, std::vector<T>& workingSet, std::string identifier2, std::vector<T>& expected)
    {
        std::cout << "Initialize data." << std::endl;
        //BOOST_TEST_MESSAGE("Initialize data.");

        try
        {
            TestCaseLoader<T> loader(file);
            return loader.getData(identifier1, workingSet, identifier2, expected);
        }
        catch (std::exception ex)
        {
            //BOOST_TEST(false, (boost::format("Exception caught loading test data. %s") % ex.what()).str());
            std::cout << "Exception occurred when loading test data. " << ex.what() << std::endl;
        }

        return false;
    }

    /** Loads the test case data.
        \param[in] file ... The file path.
        \param[in] identifier1 ... The 1st data tag identifier.
        \param[out] workingSet ... The set of working (input) data.
        \param[out] expected1 ... The 1st set of expected (output) data.
        \param[in] identifier2 ... The 2nd data tag identifier.
        \param[out] expected2 ... The 2nd set of expected (output) data.
    */
    bool initialize(const std::string& file, std::string identifier1, std::vector<std::complex<T>>& workingSet, std::vector<std::complex<T>>& expected1, 
        std::string identifier2, std::vector<std::complex<T>>& expected2)
    {
        //BOOST_TEST_MESSAGE("Initialize data.");
        std::cout << "Initialize data." << std::endl;

        try
        {
            TestCaseLoader<T> loader(file);
            return loader.getData(identifier1, workingSet, expected1, identifier2, expected2);
        }
        catch (std::exception ex)
        {
            //BOOST_TEST(false, (boost::format("Exception caught loading test data. %s") % ex.what()).str());
            std::cout << "Exception occurred when loading test data. " << ex.what() << std::endl;
        }

        return false;
    }

    /** Checks the algorithm output data against a set of expected output data.
        \param[out] workingSet ... The set of calculated data.
        \param[out] expectedOutput ... The set of expected output data.
    */
    void checkOutput(std::vector<std::complex<T>>& workingSet, std::vector<std::complex<T>>& expectedOutput)
    {
        ASSERT_EQ(workingSet.size(), expectedOutput.size()) << "The lengths of both vectors need to be equal.";

        for (auto i = 0; i < expectedOutput.size(); ++i)
        {
            bool realMatch = std::abs(workingSet[i].real() - expectedOutput[i].real()) <= static_cast<T>(kPrecision);
            bool imagMatch = std::abs(workingSet[i].imag() - expectedOutput[i].imag()) <= static_cast<T>(kPrecision);
            EXPECT_TRUE(realMatch && imagMatch)
                << std::format("Mismatch at position {}: ({}, {}i) != ({}, {}i)",
                    i,
                    workingSet[i].real(),
                    workingSet[i].imag(),
                    expectedOutput[i].real(),
                    expectedOutput[i].imag());
        }
    }
};

}

#endif // JEANBAPTISTE_UTILITIES_ALGORITHMRESULTANALYSIS_HPP_