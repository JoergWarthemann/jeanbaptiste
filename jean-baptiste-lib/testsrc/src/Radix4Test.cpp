#include <complex>
#include <format>
#include <string_view>

#include <gtest/gtest.h>

#include "AlgorithmFactory.hpp"
#include "AlgorithmFixture.hpp"
#include "Options.hpp"

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
    }

    ~Radix4Test() override = default;

protected:
    // Use factory to define algorithm Radix-4 DIT FFT algorithms for sample counts 4 ... 256.
    jb::AlgorithmFactory<
        1, 4,
        jb::options::Radix_4,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Forward,
        jb::options::Window_None,
        jb::options::Normalization_Square_Root,
        std::complex<double>>
        mFFTDITFactory;

    // Use factory to define algorithm Radix-4 DIF FFT algorithms for sample counts 4 ... 256.
    jb::AlgorithmFactory<
        1, 4,
        jb::options::Radix_4,
        jb::options::Decimation_In_Frequency,
        jb::options::Direction_Forward,
        jb::options::Window_None,
        jb::options::Normalization_Square_Root,
        std::complex<double>>
        mFFTDIFFactory;

    // Use factory to define algorithm Radix-4 DIT IFFT algorithms for sample counts 4 ... 256.
    jb::AlgorithmFactory<
        1, 4,
        jb::options::Radix_4,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Backward,
        jb::options::Window_None,
        jb::options::Normalization_Square_Root,
        std::complex<double>>
        mIFFTDITFactory;

    // Use factory to define algorithm Radix-4 DIF IFFT algorithms for sample counts 4 ... 256.
    jb::AlgorithmFactory<
        1, 4,
        jb::options::Radix_4,
        jb::options::Decimation_In_Frequency,
        jb::options::Direction_Backward,
        jb::options::Window_None,
        jb::options::Normalization_Square_Root,
        std::complex<double>>
        mIFFTDIFFactory;
};

class Radix4OneOffTest : public Test {
protected:
    // Use factory to define algorithm Radix-4 DIT FFT algorithms for sample counts 4 ... 256.
    jb::AlgorithmFactory<
        1, 4,
        jb::options::Radix_4,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Forward,
        jb::options::Window_None,
        jb::options::Normalization_Square_Root,
        std::complex<double>>
        mFFTDITFactory;
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
    verifyAlgorithm(mFFTDITFactory.getAlgorithm(GetParam().stage));
}

TEST_P(Radix4Test, SuccessfullyCalculatesInverseRadix4Dit)
{
    EXPECT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetB,
        "fft.out", mDataSetA));

    // Run algorithm with mDataSetA and check result against mDatasetB
    verifyAlgorithm(mIFFTDITFactory.getAlgorithm(GetParam().stage));
}

TEST_P(Radix4Test, SuccessfullyCalculatesRadix4Dif)
{
    EXPECT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetA,
        "fft.out", mDataSetB));

    // Run algorithm with mDataSetA and check result against mDatasetB
    verifyAlgorithm(mFFTDIFFactory.getAlgorithm(GetParam().stage));
}

TEST_P(Radix4Test, SuccessfullyCalculatesInverseRadix4Dif)
{
    EXPECT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in", mDataSetB,
        "fft.out", mDataSetA));

    // Run algorithm with mDataSetA and check result against mDatasetB
    verifyAlgorithm(mIFFTDIFFactory.getAlgorithm(GetParam().stage));
}

TEST_P(Radix4Test, SuccessfullyReconstructsInputWithRadix4DitRoundTrip)
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

TEST_P(Radix4Test, SuccessfullyReconstructsInputWithRadix4DifRoundTrip)
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

TEST_P(Radix4Test, MatchesRadix4DitAndDifForwardTransformation)
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

TEST_P(Radix4Test, MatchesRadix4DitAndDifBackwardTransformation)
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

TEST_F(Radix4OneOffTest, SuccessfullyCalculatesSmallestSupportedStage)
{
    using Complex = std::complex<double>;

    std::vector<Complex> data{
        {1.0, 0.0},
        {-1.0, 0.0},
        {1.0, 0.0},
        {-1.0, 0.0},
    };

    std::vector<Complex> expected{
        {0.0, 0.0},
        {0.0, 0.0},
        {2.0, 0.0},
        {0.0, 0.0},
    };

    mFFTDITFactory.getAlgorithm(1)->operator()(data);

    AlgorithmResultAnalysis::compareComplexDataSets(data, expected);
}

// Impulse and Constant Signals: These use stage 1, i.e. four samples. With square-root normalization, the scale factor
// is 1/sqrt(4)=0.5
TEST_F(Radix4OneOffTest, ForwardTransformationOfImpulseCalculatesFlatSpectrum)
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

    mFFTDITFactory.getAlgorithm(1)->operator()(data);

    AlgorithmResultAnalysis::compareComplexDataSets(data, expected);
}

TEST_F(Radix4OneOffTest, ForwardTransformationOfConstantSignalCalculatesOnlyDcComponent)
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

    mFFTDITFactory.getAlgorithm(1)->operator()(data);

    AlgorithmResultAnalysis::compareComplexDataSets(data, expected);
}

// Invalid Factory Stage Death Test: This only makes sense when assertions are enabled. In NDEBUG builds, the assert() in getAlgorithm() disappears, so the test should skip.
TEST_F(Radix4OneOffTest, DiesWhenUsedWithUnknownStage)
{
#if defined(NDEBUG)
    GTEST_SKIP() << "AlgorithmFactory unknown-stage behavior is guarded by assert().";
#else
    using Complex = std::complex<double>;

    EXPECT_DEATH(
        static_cast<void>(mFFTDITFactory.getAlgorithm(5)),
        "Trying to find algorithm of unknown stage.");
#endif
}

} // namespace jb::testing