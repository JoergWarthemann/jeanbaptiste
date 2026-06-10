#include <complex>
#include <format>
#include <string_view>

#include <gtest/gtest.h>

#include "AlgorithmFactory.hpp"
#include "AlgorithmFixture.hpp"

using namespace testing;

namespace jb::testing {

struct Radix4TestConfig {
    std::string_view testFileName;
    std::size_t stage;
};

class Radix4Test
    : public TestWithParam<Radix4TestConfig>
    , public AlgorithmFixture {
public:
    Radix4Test() = default;

    void SetUp() override
    {
        mDataSetA.clear();
        mDataSetB.clear();
        mInitialized = false;
    }

    ~Radix4Test() override = default;

protected:
    bool mInitialized{};

    // Use factory to define algorithm Radix-4 DIT FFT algorithms for sample counts 4 ... 256.
    jb::AlgorithmFactory<
        1, 4,
        jb::options::Radix_4,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Forward,
        jb::options::Window_None,
        jb::options::Normalization_Square_Root,
        std::complex<double>>
        mFFTFactory;

    // Use factory to define algorithm Radix-4 DIT IFFT algorithms for sample counts 4 ... 256.
    jb::AlgorithmFactory<
        1, 4,
        jb::options::Radix_4,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Backward,
        jb::options::Window_None,
        jb::options::Normalization_Square_Root,
        std::complex<double>>
        mIFFtFactory;
};

INSTANTIATE_TEST_SUITE_P(
    Radix4Tests,
    Radix4Test,
    ::testing::Values(
        Radix4TestConfig{.testFileName = "square pulse (n=64)", .stage = 3}));

TEST_P(Radix4Test, SuccessfullyCalculatesRadix4Dit)
{
    EXPECT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetA,
        "fft.out", mDataSetB));

    // Run algorithm with mDataSetA and check result against mDatasetB
    verifyAlgorithm(mFFTFactory.getAlgorithm(GetParam().stage));
}

TEST_P(Radix4Test, SuccessfullyCalculatesInverseRadix4Dit)
{
    EXPECT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetB,
        "fft.out", mDataSetA));

    // Run algorithm with mDataSetA and check result against mDatasetB
    verifyAlgorithm(mIFFtFactory.getAlgorithm(GetParam().stage));
}

TEST_P(Radix4Test, SuccessfullyCalculatesRadix4Dif)
{
    EXPECT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetA,
        "fft.out", mDataSetB));

    // Run algorithm with mDataSetA and check result against mDatasetB
    verifyAlgorithm(mFFTFactory.getAlgorithm(GetParam().stage));
}

TEST_P(Radix4Test, SuccessfullyCalculatesInverseRadix4Dif)
{
    EXPECT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetB,
        "fft.out", mDataSetA));

    // Run algorithm with mDataSetA and check result against mDatasetB
    verifyAlgorithm(mIFFtFactory.getAlgorithm(GetParam().stage));
}

} // namespace jb::testing