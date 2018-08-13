#include <boost/format.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/test/unit_test.hpp>
#include <complex>
#include <string>
#include "TestCaseLoader.h"
#include <vector>

namespace tt = boost::test_tools;

namespace Utilities
{

template<typename T>
class AlgorithmResultAnalysis
{
    const double kPrecision = 0.00001;

public:
    /** Loads the test case data.
        \param[in] file ... The file path.
        \param[in] identifier1 ... The 1st data tag identifier.
        \param[out] workingSet ... The set of working (input) data.
        \param[in] identifier2 ... The 2nd data tag identifier.
        \param[out] expected ... The set of expected (output) data.
    */
    void initialize(const std::string& file, std::string identifier1, std::vector<T>& workingSet, std::string identifier2, std::vector<T>& expected)
    {
        BOOST_TEST_MESSAGE("Initialize data.");

        try
        {
            TestCaseLoader<T> loader(file);
            BOOST_TEST(loader.getData(identifier1, workingSet, identifier2, expected), "Error loading test data.");
        }
        catch (std::exception ex)
        {
            BOOST_TEST(false, (boost::format("Exception caught loading test data. %s") % ex.what()).str());
        }
    }

    /** Loads the test case data.
        \param[in] file ... The file path.
        \param[in] identifier1 ... The 1st data tag identifier.
        \param[out] workingSet ... The set of working (input) data.
        \param[out] expected1 ... The 1st set of expected (output) data.
        \param[in] identifier2 ... The 2nd data tag identifier.
        \param[out] expected2 ... The 2nd set of expected (output) data.
    */
    void initialize(const std::string& file, std::string identifier1, std::vector<std::complex<T>>& workingSet, std::vector<std::complex<T>>& expected1, 
        std::string identifier2, std::vector<std::complex<T>>& expected2)
    {
        BOOST_TEST_MESSAGE("Initialize data.");

        try
        {
            TestCaseLoader<T> loader(file);
            BOOST_TEST(loader.getData(identifier1, workingSet, expected1, identifier2, expected2), "Error loading test data.");
        }
        catch (std::exception ex)
        {
            BOOST_TEST(false, (boost::format("Exception caught loading test data. %s") % ex.what()).str());
        }
    }

    /** Checks the algorithm output data against a set of expected output data.
        \param[out] workingSet ... The set of calculated data.
        \param[out] expectedOutput ... The set of expected output data.
    */
    void checkOutput(std::vector<std::complex<T>>& workingSet, std::vector<std::complex<T>>& expectedOutput)
    {
        //BOOST_CHECK_EQUAL_COLLECTIONS(workingSet.begin(), workingSet.end(), expectedOutput.begin(), expectedOutput.end());

        if (workingSet.size() != expectedOutput.size())
           BOOST_TEST(false, "The lengths of both vectors need to be equal.");

        for (auto i = 0; i < expectedOutput.size(); ++i)
           BOOST_TEST(
                  ((std::abs(workingSet[i].real() - expectedOutput[i].real()) <= static_cast<T>(kPrecision))
               &&  (std::abs(workingSet[i].imag() - expectedOutput[i].imag()) <= static_cast<T>(kPrecision))),
               (boost::format("Mismatch at position %s: (%s, %si) != (%s, %si)")
                  % i
                  % workingSet[i].real()
                  % workingSet[i].imag()
                  % expectedOutput[i].real()
                  % expectedOutput[i].imag()).str());
    }
};

}