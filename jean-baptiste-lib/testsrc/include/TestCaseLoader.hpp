#ifndef JB_TESTING_TESTCASELOADER_HPP_
#define JB_TESTING_TESTCASELOADER_HPP_

#include <charconv>
#include <complex>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace jb::testing {

inline std::string_view trim(std::string_view str)
{
    auto start = str.find_first_not_of(" \t\n\r\f\v");
    if (start == std::string_view::npos) {
        return std::string_view(); // All whitespace
    }

    auto end = str.find_last_not_of(" \t\n\r\f\v");

    return str.substr(start, end - start + 1);
}

/**
 * Converts a single string with 2 numbers into a complex number.
 * \param line ... The string.
 */
template <typename T>
inline std::pair<T, T> extractComplexArguments(std::string_view line)
{
    line = trim(line);

    std::string_view::size_type posEnd = std::string_view::npos;
    std::string_view::size_type posStart = std::string_view::npos;
    if ((posStart = line.find('\t')) != std::string_view::npos) {
        T real{}, imag{};

        auto left = trim(line.substr(0, posStart));
        auto right = trim(line.substr(posStart + 1));

        real = static_cast<T>(std::strtod(left.data(), nullptr));
        imag = static_cast<T>(std::strtod(right.data(), nullptr));

        return {real, imag};
    }

    throw std::invalid_argument("Invalid complex number format");
}

/**
 * Converts a single string with 1 number into a number.
 * \param line ... The string.
 */
template <typename T>
inline T extractRealNumber(std::string& line)
{
    T result = T(0);

    line = trim(line);
    result = static_cast<T>(std::strtod(line.data(), nullptr));

    return result;
}

class TestCaseLoader {
public:
    TestCaseLoader(std::string_view file)
        : mFile(std::filesystem::absolute(file))
    {}

    ~TestCaseLoader(void)
    {}

    /**
     * Opens the file and reads complex the test case data in.
     * \param tagA ... Identifier of data set A.
     * \param dataSetA ... The data set A.
     * \param tagB ... Identifier of data set B.
     * \param dataSetB ... The data set B.
     */
    bool getData(std::string_view tagA, std::vector<std::complex<double>>& dataSetA, std::string_view tagB,
        std::vector<std::complex<double>>& dataSetB);

    /**
     * Opens the file and reads the floating point test case data in.
     * \param tagA ... Identifier of data set A.
     * \param dataSetA ... The data set A.
     * \param tagB ... Identifier of data set B.
     * \param dataSetB ... The data set B.
     */
    bool getData(std::string_view tagA, std::vector<double>& dataSetA, std::string_view tagB,
        std::vector<double>& dataSetB);

private:
    std::filesystem::path mFile{};
};

} // namespace jb::testing

#endif // JB_TESTING_TESTCASELOADER_HPP_