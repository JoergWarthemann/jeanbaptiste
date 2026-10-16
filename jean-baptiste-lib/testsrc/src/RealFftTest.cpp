#include <complex>
#include <vector>

#include <gtest/gtest.h>

#include "AlgorithmFactory.hpp"
#include "AlgorithmResultAnalysis.hpp"
#include "Options.hpp"

namespace jb::testing {

class RealFftTest : public ::testing::Test {
protected:
    using Complex = std::complex<double>;

    static std::vector<Complex> createRealInput(const std::size_t sampleCount)
    {
        std::vector<Complex> data(sampleCount);

        for (auto idx = 0U; idx < sampleCount; ++idx) {
            data[idx].real(static_cast<double>((idx % 5) - 2) + static_cast<double>(idx % 3) * 0.25);
        }

        return data;
    }

    jb::AlgorithmFactory<
        2, 5,
        jb::options::Radix_2,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Forward,
        jb::options::Window_None,
        jb::options::Normalization_No,
        Complex>
        mComplexFftFactory;

    jb::AlgorithmFactory<
        2, 5,
        jb::options::Radix_2,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Forward,
        jb::options::Window_None,
        jb::options::Normalization_No,
        Complex,
        jb::options::Data_Real>
        mRealFftFactory;

    jb::AlgorithmFactory<
        2, 5,
        jb::options::Radix_2,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Backward,
        jb::options::Window_None,
        jb::options::Normalization_Division_By_Length,
        Complex,
        jb::options::Data_Real>
        mRealIfftFactory;
};

TEST_F(RealFftTest, ForwardRealFftMatchesComplexFftForRealInput)
{
    std::vector<Complex> realData{
        {1.0, 0.0},
        {2.0, 0.0},
        {-1.0, 0.0},
        {0.5, 0.0},
        {3.0, 0.0},
        {-2.0, 0.0},
        {0.0, 0.0},
        {1.5, 0.0},
    };
    auto complexData = realData;

    (*mRealFftFactory.getAlgorithm(3))(realData);
    (*mComplexFftFactory.getAlgorithm(3))(complexData);

    AlgorithmResultAnalysis::compareComplexDataSets(realData, complexData);
}

TEST_F(RealFftTest, RealFftRoundTripReconstructsInput)
{
    std::vector<Complex> data{
        {1.0, 0.0},
        {2.0, 0.0},
        {-1.0, 0.0},
        {0.5, 0.0},
        {3.0, 0.0},
        {-2.0, 0.0},
        {0.0, 0.0},
        {1.5, 0.0},
    };
    auto expected = data;

    (*mRealFftFactory.getAlgorithm(3))(data);
    (*mRealIfftFactory.getAlgorithm(3))(data);

    AlgorithmResultAnalysis::compareComplexDataSets(data, expected);
}

TEST_F(RealFftTest, ForwardSplitRadixRealFftMatchesComplexFftForRealInput)
{
    jb::AlgorithmFactory<
        3, 4,
        jb::options::Radix_Split_2_4,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Forward,
        jb::options::Window_None,
        jb::options::Normalization_No,
        Complex,
        jb::options::Data_Real>
        realFftFactory;

    auto realData = createRealInput(8);
    auto complexData = realData;

    (*realFftFactory.getAlgorithm(3))(realData);
    (*mComplexFftFactory.getAlgorithm(3))(complexData);

    AlgorithmResultAnalysis::compareComplexDataSets(realData, complexData);
}

TEST_F(RealFftTest, ForwardRadix4RealFftMatchesComplexFftForRealInput)
{
    jb::AlgorithmFactory<
        5, 6,
        jb::options::Radix_2,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Forward,
        jb::options::Window_None,
        jb::options::Normalization_No,
        Complex>
        complexFftFactory;
    jb::AlgorithmFactory<
        2, 3,
        jb::options::Radix_4,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Forward,
        jb::options::Window_None,
        jb::options::Normalization_No,
        Complex,
        jb::options::Data_Real>
        realFftFactory;

    auto realData = createRealInput(32);
    auto complexData = realData;

    (*realFftFactory.getAlgorithm(2))(realData);
    (*complexFftFactory.getAlgorithm(5))(complexData);

    AlgorithmResultAnalysis::compareComplexDataSets(realData, complexData);
}

} // namespace jb::testing