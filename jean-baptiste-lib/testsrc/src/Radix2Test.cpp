#include <complex>
#include <format>
#include <string_view>

#include <gtest/gtest.h>

#include "AlgorithmFactory.hpp"
#include "AlgorithmFixture.hpp"

using namespace testing;

namespace jb::testing {

struct Radix2TestConfig {
    std::string_view testFileName;
    std::size_t stage;
};

class Radix2Test
    : public TestWithParam<Radix2TestConfig>
    , public AlgorithmFixture {
public:
    Radix2Test() = default;

    void SetUp() override
    {
        mDataSetA.clear();
        mDataSetB.clear();
        mInitialized = false;
    }

    ~Radix2Test() override = default;

protected:
    bool mInitialized{};

    // Use factory to define algorithm Radix-2 DIT FFT algorithms for sample counts 64 ... 256.
    jb::AlgorithmFactory<
        6, 8,
        jb::options::Radix_2,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Forward,
        jb::options::Window_None,
        jb::options::Normalization_Square_Root,
        std::complex<double>>
        mFFTFactory;

    // Use factory to define algorithm Radix-2 DIT IFFT algorithms for sample counts 64 ... 256.
    jb::AlgorithmFactory<
        6, 8,
        jb::options::Radix_2,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Backward,
        jb::options::Window_None,
        jb::options::Normalization_Square_Root,
        std::complex<double>>
        mIFFtFactory;
};

INSTANTIATE_TEST_SUITE_P(
    Radix2Tests,
    Radix2Test,
    ::testing::Values(
        Radix2TestConfig{.testFileName = "square pulse (n=64)", .stage = 6},
        Radix2TestConfig{.testFileName = "square pulse (n=128)", .stage = 7},
        Radix2TestConfig{.testFileName = "cosine (n=128)", .stage = 7}));

TEST_P(Radix2Test, SuccessfullyCalculatesRadix2Dit)
{
    EXPECT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetA,
        "fft.out", mDataSetB));

    // Run algorithm with mDataSetA and check result against mDatasetB
    verifyAlgorithm(mFFTFactory.getAlgorithm(GetParam().stage));
}

TEST_P(Radix2Test, SuccessfullyCalculatesInverseRadix2Dit)
{
    EXPECT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetB,
        "fft.out", mDataSetA));

    // Run algorithm with mDataSetA and check result against mDatasetB
    verifyAlgorithm(mIFFtFactory.getAlgorithm(GetParam().stage));
}

TEST_P(Radix2Test, SuccessfullyCalculatesRadix2Dif)
{
    EXPECT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetA,
        "fft.out", mDataSetB));

    // Run algorithm with mDataSetA and check result against mDatasetB
    verifyAlgorithm(mFFTFactory.getAlgorithm(GetParam().stage));
}

TEST_P(Radix2Test, SuccessfullyCalculatesInverseRadix2Dif)
{
    EXPECT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetB,
        "fft.out", mDataSetA));

    // Run algorithm with mDataSetA and check result against mDatasetB
    verifyAlgorithm(mIFFtFactory.getAlgorithm(GetParam().stage));
}

} // namespace jb::testing