#include <complex>
#include <format>
#include <string_view>

#include <gtest/gtest.h>

#include "AlgorithmFactory.hpp"
#include "AlgorithmFixture.hpp"
#include "Options.hpp"

using namespace testing;

namespace jb::testing {

struct RadixSplit24TestConfig {
    std::string_view testFileName;
    std::size_t stage;
};

class RadixSplit24Test
    : public TestWithParam<RadixSplit24TestConfig>
    , public AlgorithmFixture {
public:
    RadixSplit24Test() = default;

    void SetUp() override
    {
        mDataSetA.clear();
        mDataSetB.clear();
        mInitialized = false;
    }

    ~RadixSplit24Test() override = default;

protected:
    bool mInitialized{};

    // // Create Split-Radix-2-4 DIF FFT algorithms for sample counts 2 ... 256.
    //     jb::AlgorithmFactory<1, 8, jbo::Radix_Split_2_4, jbo::Decimation_In_Frequency, jbo::Direction_Forward, jbo::Window_None,
    //         jbo::Normalization_Square_Root, std::complex<double>> fftFactory;

    //     // Create Split-Radix-2-4 DIF IFFT algorithms for sample counts 2 ... 256.
    //     jb::AlgorithmFactory<1, 8, jbo::Radix_Split_2_4, jbo::Decimation_In_Frequency, jbo::Direction_Backward, jbo::Window_None,
    //         jbo::Normalization_Square_Root, std::complex<double>> ifftFactory;

    // Use factory to define algorithm Split-Radix-2-4 DIT FFT algorithms for sample counts 2 ... 256.
    jb::AlgorithmFactory<
        1, 8,
        jb::options::Radix_Split_2_4,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Forward,
        jb::options::Window_None,
        jb::options::Normalization_Square_Root,
        std::complex<double>>
        mFFTDITFactory;

    // Use factory to define algorithm Split-Radix-2-4 DIF FFT algorithms for sample counts 2 ... 256.
    jb::AlgorithmFactory<
        1, 8,
        jb::options::Radix_Split_2_4,
        jb::options::Decimation_In_Frequency,
        jb::options::Direction_Forward,
        jb::options::Window_None,
        jb::options::Normalization_Square_Root,
        std::complex<double>>
        mFFTDIFFactory;

    // Use factory to define algorithm Split-Radix-2-4 DIT IFFT algorithms for sample counts 2 ... 256.
    jb::AlgorithmFactory<
        1, 8,
        jb::options::Radix_Split_2_4,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Backward,
        jb::options::Window_None,
        jb::options::Normalization_Square_Root,
        std::complex<double>>
        mIFFTDITFactory;

    // Use factory to define algorithm Split-Radix-2-4 DIF IFFT algorithms for sample counts 2 ... 256.
    jb::AlgorithmFactory<
        1, 8,
        jb::options::Radix_Split_2_4,
        jb::options::Decimation_In_Frequency,
        jb::options::Direction_Backward,
        jb::options::Window_None,
        jb::options::Normalization_Square_Root,
        std::complex<double>>
        mIFFTDIFFactory;
};

class RadixSplit24OneOffTest : public Test {
protected:
    // Use factory to define algorithm Split-Radix-2-4 DIT FFT algorithms for sample counts 2 ... 256.
    jb::AlgorithmFactory<
        1, 8,
        jb::options::Radix_Split_2_4,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Forward,
        jb::options::Window_None,
        jb::options::Normalization_Square_Root,
        std::complex<double>>
        mFFTDITFactory;
};

INSTANTIATE_TEST_SUITE_P(
    RadixSplit24Tests,
    RadixSplit24Test,
    ::testing::Values(
        RadixSplit24TestConfig{.testFileName = "square pulse (n=64)", .stage = 6}));

TEST_P(RadixSplit24Test, SuccessfullyCalculatesRadixSplit24Dit)
{
    EXPECT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetA,
        "fft.out", mDataSetB));

    // Run algorithm with mDataSetA and check result against mDatasetB
    verifyAlgorithm(mFFTDITFactory.getAlgorithm(GetParam().stage));
}

TEST_P(RadixSplit24Test, SuccessfullyCalculatesInverseRadixSplit24Dit)
{
    EXPECT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetB,
        "fft.out", mDataSetA));

    // Run algorithm with mDataSetA and check result against mDatasetB
    verifyAlgorithm(mIFFTDITFactory.getAlgorithm(GetParam().stage));
}

TEST_P(RadixSplit24Test, SuccessfullyCalculatesRadixSplit24Dif)
{
    EXPECT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetA,
        "fft.out", mDataSetB));

    // Run algorithm with mDataSetA and check result against mDatasetB
    verifyAlgorithm(mFFTDIFFactory.getAlgorithm(GetParam().stage));
}

TEST_P(RadixSplit24Test, SuccessfullyCalculatesInverseRadixSplit24Dif)
{
    EXPECT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetB,
        "fft.out", mDataSetA));

    // Run algorithm with mDataSetA and check result against mDatasetB
    verifyAlgorithm(mIFFTDIFFactory.getAlgorithm(GetParam().stage));
}

TEST_P(RadixSplit24Test, SuccessfullyReconstructsInputWithRadixSplit24DitRoundTrip)
{
    ASSERT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetA,
        "fft.out", mDataSetB));

    auto expected = mDataSetA;
    auto fftDIT = mFFTDITFactory.getAlgorithm(GetParam().stage);
    auto iFFTDIT = mIFFTDITFactory.getAlgorithm(GetParam().stage);

    (*fftDIT)(mDataSetA);
    (*iFFTDIT)(mDataSetA);

    AlgorithmResultAnalysis::compareComplexDataSets(mDataSetA, expected);
}

TEST_P(RadixSplit24Test, SuccessfullyReconstructsInputWithRadixSplit24DifRoundTrip)
{
    ASSERT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetA,
        "fft.out", mDataSetB));

    auto expected = mDataSetA;
    auto fftDIF = mFFTDIFFactory.getAlgorithm(GetParam().stage);
    auto iFFTDIF = mIFFTDIFFactory.getAlgorithm(GetParam().stage);

    (*fftDIF)(mDataSetA);
    (*iFFTDIF)(mDataSetA);

    AlgorithmResultAnalysis::compareComplexDataSets(mDataSetA, expected);
}

TEST_P(RadixSplit24Test, MatchesRadixSplit24DitAndDifForwardTransformation)
{
    ASSERT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetA,
        "fft.out", mDataSetB));

    auto ditData = mDataSetA;
    auto difData = mDataSetA;
    auto fftDIT = mFFTDITFactory.getAlgorithm(GetParam().stage);
    auto fftDIF = mFFTDIFFactory.getAlgorithm(GetParam().stage);

    (*fftDIT)(ditData);
    (*fftDIF)(difData);

    AlgorithmResultAnalysis::compareComplexDataSets(ditData, difData);
}

TEST_P(RadixSplit24Test, MatchesRadixSplit24DitAndDifBackwardTransformation)
{
    ASSERT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetA,
        "fft.out", mDataSetB));

    auto ditData = mDataSetB;
    auto difData = mDataSetB;
    auto iFFTDIT = mIFFTDITFactory.getAlgorithm(GetParam().stage);
    auto iFFTDIF = mIFFTDIFFactory.getAlgorithm(GetParam().stage);

    (*iFFTDIT)(ditData);
    (*iFFTDIF)(difData);

    AlgorithmResultAnalysis::compareComplexDataSets(ditData, difData);
}

TEST_F(RadixSplit24OneOffTest, SuccessfullyCalculatesSmallestSupportedStage)
{
    using Complex = std::complex<double>;

    std::vector<Complex> data{
        {1.0, 0.0},
        {-1.0, 0.0},
    };

    std::vector<Complex> expected{
        {0.0, 0.0},
        {std::sqrt(2.0), 0.0},
    };

    mFFTDITFactory.getAlgorithm(1)->operator()(data);

    AlgorithmResultAnalysis::compareComplexDataSets(data, expected);
}

// Impulse and Constant Signals: These use stage 1, i.e. four samples. With square-root normalization, the scale factor
// is 1/sqrt(4)=0.5
TEST_F(RadixSplit24OneOffTest, ForwardTransformationOfImpulseCalculatesFlatSpectrum)
{
    using Complex = std::complex<double>;

    std::vector<Complex> data{
        {1.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
    };

    std::vector<Complex> expected{
        {0.5, 0.0},
        {0.5, 0.0},
        {0.5, 0.0},
        {0.5, 0.0},
    };

    mFFTDITFactory.getAlgorithm(2)->operator()(data);

    AlgorithmResultAnalysis::compareComplexDataSets(data, expected);
}

TEST_F(RadixSplit24OneOffTest, ForwardTransformationOfConstantSignalCalculatesOnlyDcComponent)
{
    using Complex = std::complex<double>;

    std::vector<Complex> data{
        {1.0, 0.0},
        {1.0, 0.0},
        {1.0, 0.0},
        {1.0, 0.0},
    };

    std::vector<Complex> expected{
        {2.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
        {0.0, 0.0},
    };

    mFFTDITFactory.getAlgorithm(2)->operator()(data);

    AlgorithmResultAnalysis::compareComplexDataSets(data, expected);
}

// Invalid Factory Stage Death Test: This only makes sense when assertions are enabled. In NDEBUG builds, the assert() in getAlgorithm() disappears, so the test should skip.
TEST_F(RadixSplit24OneOffTest, DiesWhenUsedWithUnknownStage)
{
#if defined(NDEBUG)
    GTEST_SKIP() << "AlgorithmFactory unknown-stage behavior is guarded by assert().";
#else
    using Complex = std::complex<double>;

    EXPECT_DEATH(
        static_cast<void>(mFFTDITFactory.getAlgorithm(9)),
        "Trying to find algorithm of unknown stage.");
#endif
}

} // namespace jb::testing