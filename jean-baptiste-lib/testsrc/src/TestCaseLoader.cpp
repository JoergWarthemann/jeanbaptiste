#include "TestCaseLoader.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace jb::testing {

bool TestCaseLoader::getData(std::string_view tagA, std::vector<std::complex<double>>& dataSetA,
    std::string_view tagB, std::vector<std::complex<double>>& dataSetB)
{
    if (std::filesystem::exists(mFile)) {
        // Create an empty property tree object
        boost::property_tree::ptree tree;
        boost::property_tree::read_xml(mFile.string(), tree);

        std::string line;

        // Read data into dataSetA..
        std::stringstream linesIn1(tree.get<std::string>(std::string(tagA)));
        while (std::getline(linesIn1, line)) {
            boost::algorithm::trim(line);
            if (!line.empty()) {
                auto [real, imag] = extractComplexArguments<double>(line);
                dataSetA.emplace_back(real, imag);
            }
        }

        // Read data into dataSetB..
        std::stringstream linesIn2(tree.get<std::string>(std::string(tagB)));
        while (std::getline(linesIn2, line)) {
            boost::algorithm::trim(line);
            if (!line.empty()) {
                auto [real, imag] = extractComplexArguments<double>(line);
                dataSetB.emplace_back(real, imag);
            }
        }

        return ((dataSetA.size() > 0) && (dataSetA.size() == dataSetB.size()));
    }

    return false;
}

bool TestCaseLoader::getData(std::string_view tagA, std::vector<double>& dataSetA, std::string_view tagB, std::vector<double>& dataSetB)
{
    if (std::filesystem::exists(mFile)) {
        // Create an empty property tree object
        boost::property_tree::ptree tree;
        boost::property_tree::read_xml(mFile.string(), tree);

        std::string line;

        // Read data into dataSetA.
        std::stringstream linesIn1(tree.get<std::string>(std::string(tagA)));
        while (std::getline(linesIn1, line)) {
            boost::algorithm::trim(line);
            if (!line.empty()) {
                dataSetA.emplace_back(extractRealNumber<double>(line));
            }
        }

        // Read data into dataSetB.
        std::stringstream linesIn2(tree.get<std::string>(std::string(tagB)));
        while (std::getline(linesIn2, line)) {
            boost::algorithm::trim(line);
            if (!line.empty()) {
                dataSetB.emplace_back(extractRealNumber<double>(line));
            }
        }

        return ((dataSetA.size() > 0) && (dataSetA.size() == dataSetB.size()));
    }

    return false;
}

} // namespace jb::testing