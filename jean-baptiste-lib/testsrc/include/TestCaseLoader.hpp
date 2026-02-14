#ifndef JEANBAPTISTE_UTILITIES_TESTCASELOADER_HPP_
#define JEANBAPTISTE_UTILITIES_TESTCASELOADER_HPP_

#include <charconv>
#include <complex>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
//#include <boost/algorithm/string.hpp>
//#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace Utilities {

inline std::string_view trim(std::string_view str) {
    auto start = str.find_first_not_of(" \t\n\r\f\v");
    if (start == std::string_view::npos) {
        return std::string_view(); // All whitespace
    }
    
    auto end = str.find_last_not_of(" \t\n\r\f\v");
    
    return str.substr(start, end - start + 1);
}

template<typename T>
class TestCaseLoader
{
    //boost::filesystem::path file_;
    std::filesystem::path mFile{};
    //std::string file_;

    /** Converts a single string with 2 numbers into a complex number.
        \param[in] line ... The string.
    */
    //std::complex<T> extractComplexNumber(std::string& line)
    std::pair<T, T> extractComplexArguments(std::string_view line)
    {
        line = trim(line);

        std::string_view::size_type posEnd = std::string_view::npos;
        std::string_view::size_type posStart = std::string_view::npos;
        if ((posStart = line.find('\t')) != std::string_view::npos) {
            T real, imag;
            std::string_view tmp = trim(line.substr(0, posStart));
            std::from_chars(tmp.data(), tmp.data() + tmp.size(), real);

            // std::istringstream iss(tmp);
            // iss >> real;
            
            // iss.clear();

            // tmp = trim(line.substr(posStart, posEnd - posStart));
            // iss.str(tmp);
            // iss >> imag;

            std::from_chars(tmp.data(), tmp.data() + tmp.size(), imag);

            return {real, imag};
        }

        throw std::invalid_argument("Invalid complex number format");
    }
    
    /** Converts a single string with 1 number into a number.
        \param[in] line ... The string.
    */
    T extractRealNumber(std::string& line)
    {
        T result = T(0);

        line = trim(line);
        std::from_chars(line.data(), line.data() + line.size(), result);

        // boost::algorithm::trim(line);
        // std::istringstream iss(line);
        // iss >> result;

        return result;
    }

public:
    TestCaseLoader(std::string_view file)
        : mFfile(std::filesystem::system_complete(file))
    {}

    ~TestCaseLoader(void)
    {}

    /** Opens the file and reads the test case data in.
        \param[in] identifier1 ... The 1st data tag identifier.
        \param[out] dataIn ... The set of input data.
        \param[out] expected1 ... The set of expected 1st (output) data.
        \param[out] identifier2 ... The 2nd data tag identifier.
        \param[out] expected2 ... The set of expected 2nd (output) data.
    */
   // TODO: turn identifier1 into string_view, keep the vectors
    bool getData(std::string_view identifier1, std::vector< std::complex<T> >& dataIn, std::vector< std::complex<T> >& expected1,
        std::string_view identifier2, std::vector< std::complex<T> >& expected2)
    {
        if (std::filesystem::exists(file_)) {
            // Create an empty property tree object
            boost::property_tree::ptree tree;
            boost::property_tree::read_xml(file_.string(), tree);

            std::string line;

            // Load dataIn and expected2 data.
            std::stringstream linesIn1(tree.get<std::string>(identifier1));
            while (std::getline(linesIn1, line)) {
                boost::algorithm::trim(line);
                if (!line.empty()) {
                    auto [real, imag] = extractComplexArguments(line);
                    dataIn.emplace_back(real, imag);
                    expected1.emplace_back(real, imag);

                    // dataIn.push_back(tmp);
                    // expected1.push_back(tmp);
                }
            }

            // Load expected2 data.
            std::stringstream linesIn2(tree.get<std::string>(identifier2));
            while (std::getline(linesIn2, line)) {
                boost::algorithm::trim(line);
                if (!line.empty()) {
                    expected2.emplace_back(extractComplexArguments(line));
                }
            }

            if ((dataIn.size() > 0)
                && (dataIn.size() == expected1.size())
                && (expected2.size() == expected1.size()))
                return true;
        }

        return false;
    }

    /** Opens the file and reads the test case data in.
        \param[in] identifier1 ... The 1st data tag identifier.
        \param[out] dataIn ... The set of input data.
        \param[in] identifier2 ... The 2nd data tag identifier.
        \param[out] expected ... The set of expected 2nd (output) data.
    */
    bool getData(std::string_view identifier1, std::vector<T>& dataIn, std::string_view identifier2, std::vector<T>& expected)
    {
        if (std::filesystem::exists(mFfile)) {
            // Create an empty property tree object
            boost::property_tree::ptree tree;
            boost::property_tree::read_xml(file_.string(), tree);

            std::string line;

            // Load dataIn.
            std::stringstream linesIn1(tree.get<std::string>(identifier1));
            while (std::getline(linesIn1, line)) {
                boost::algorithm::trim(line);
                if (!line.empty())
                    dataIn.emplace_back(extractRealNumber(line));
            }

            // Load expected data.
            std::stringstream linesIn2(tree.get<std::string>(identifier2));
            while (std::getline(linesIn2, line))
            {
                boost::algorithm::trim(line);
                if (!line.empty()) {
                    expected.emplace_back(extractRealNumber(line));
                }
            }

            if ((dataIn.size() > 0)
                && (dataIn.size() == expected.size())) {
                return true;
            }
        }

        return false;
    }
};

}

#endif // JEANBAPTISTE_UTILITIES_TESTCASELOADER_HPP_