#pragma once

#include <vector>
#include <complex>
#include <string>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "StringConversion.h"

namespace Utilities {

template<typename T>
class TestCaseLoader
{
    boost::filesystem::path file_;
    //std::string file_;

    /** Converts a single string with 2 numbers into a complex number.
        \param[in] line ... The string.
    */
    std::complex<T> extractComplexNumber(std::string& line)
    {
        boost::algorithm::trim(line);

        std::string::size_type posEnd = std::string::npos;
        std::string::size_type posStart = std::string::npos;
        if    ((posStart = line.find('\t')) != std::string::npos)
        {
            T real, imag;
            std::string tmp = line.substr(0, posStart);
            boost::algorithm::trim(tmp);
            std::istringstream iss(tmp);
            iss >> real;
            
            iss.clear();

            tmp = line.substr(posStart, posEnd - posStart);
            boost::algorithm::trim(tmp);
            iss.str(tmp);
            iss >> imag;

            return std::complex<T>(real, imag);
        }

        return std::complex<T>(0, 0);
    }
    
    /** Converts a single string with 1 number into a number.
        \param[in] line ... The string.
    */
    T extractRealNumber(std::string& line)
    {
        T result = T(0);

        boost::algorithm::trim(line);
        std::istringstream iss(line);
        iss >> result;

        return result;
    }

public:
    TestCaseLoader(const std::string& file)
        : file_(boost::filesystem::system_complete(file))
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
    bool getData(std::string identifier1, std::vector< std::complex<T> >& dataIn, std::vector< std::complex<T> >& expected1,
        std::string identifier2, std::vector< std::complex<T> >& expected2)
    {
        if (boost::filesystem::exists(file_))
        {
            // Create an empty property tree object
            boost::property_tree::ptree tree;
            boost::property_tree::read_xml(file_.string(), tree);

            std::string line;

            // Load dataIn and expected2 data.
            std::stringstream linesIn1(tree.get<std::string>(identifier1));
            while (std::getline(linesIn1, line))
            {
                boost::algorithm::trim(line);
                if (!line.empty())
                {
                    std::complex<T> tmp = extractComplexNumber(line);
                    dataIn.push_back(tmp);
                    expected1.push_back(tmp);
                }
            }

            // Load expected2 data.
            std::stringstream linesIn2(tree.get<std::string>(identifier2));
            while (std::getline(linesIn2, line))
            {
                boost::algorithm::trim(line);
                if (!line.empty())
                    expected2.push_back(extractComplexNumber(line));
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
    bool getData(std::string identifier1, std::vector<T>& dataIn, std::string identifier2, std::vector<T>& expected)
    {
        if (boost::filesystem::exists(file_))
        {
            // Create an empty property tree object
            boost::property_tree::ptree tree;
            boost::property_tree::read_xml(file_.string(), tree);

            std::string line;

            // Load dataIn.
            std::stringstream linesIn1(tree.get<std::string>(identifier1));
            while (std::getline(linesIn1, line))
            {
                boost::algorithm::trim(line);
                if (!line.empty())
                    dataIn.push_back(extractRealNumber(line));
            }

            // Load expected data.
            std::stringstream linesIn2(tree.get<std::string>(identifier2));
            while (std::getline(linesIn2, line))
            {
                boost::algorithm::trim(line);
                if (!line.empty())
                    expected.push_back(extractRealNumber(line));
            }

            if ((dataIn.size() > 0)
                && (dataIn.size() == expected.size()))
                return true;
        }

        return false;
    }
};

}