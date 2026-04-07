#include <algorithm>
#include <array>
#include <complex>
#include <format>
#include <ranges>
#include <span>
#include <string_view>

#include <gtest/gtest.h>

#include "core/Radix2.hpp"
#include "AlgorithmFixture.hpp"
#include "AlgorithmFactory.hpp"

using namespace testing;

namespace jeanbaptiste::testing {

struct Radix2TestConfig {
    std::string_view testFileName;
};

class Radix2Test
    //: public ::testing::Test
    : public TestWithParam<Radix2TestConfig>
    , public AlgorithmFixture {
public:
    Radix2Test()
    {
        // mInputCopy = mInput;

        // constexpr std::string_view testName{"square pulse (n=64)"};
        // mInitialized = mAlgorithmResult.initialize(std::format("{}/{}.xml", TEST_DATA_DIR, testName),
        //     "fft.in",
        //     mInput,
        //     mInputCopy,
        //     "fft.out",
        //     mOutput);
    }

    ~Radix2Test() override = default;

protected:
    bool mInitialized{};
    //std::vector<std::complex<double>> mDataSetACopy;
};

INSTANTIATE_TEST_SUITE_P(
    Radix2Tests,
    Radix2Test,
    ::testing::Values(
        Radix2TestConfig{.testFileName = "square pulse (n=64)"}
    )
);

TEST_P(Radix2Test, SuccessfullyCalculatesRadix2Dit)
{
    mDataSetA.clear();
    //mDataSetACopy.clear();
    mDataSetB.clear();

    EXPECT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in",
        mDataSetA,
        "fft.out",
        mDataSetB));
    
    // Use factory to define algorithm.
    jeanbaptiste::
    // Run algorithm with mDataSetA
    // Check output against mDataSetB
}

TEST_P(Radix2Test, SuccessfullyCalculatesInverseRadix2Dit)
{
    mDataSetA.clear();
    //mDataSetACopy.clear();
    mDataSetB.clear();

    AlgorithmResultAnalysis::initialize(std::format("{}/{}.xml", TEST_DATA_DIR, GetParam().testFileName),
        "fft.in",
        mDataSetA,
        "fft.out",
        mDataSetB);
}

TEST_F(Radix2Test, SuccessfullyCalculatesRadix2Dif)
{
    GTEST_SKIP() << "TODO: Add Radix2 test for DIF FFT";
}

} // namespace jeanbaptiste::testing