#include <cmath>
#include <complex>
#include <format>
#include <span>
#include <string_view>

#include <gtest/gtest.h>

#include "AlgorithmFactory.hpp"
#include "AlgorithmFixture.hpp"
#include "Options.hpp"

using namespace testing;

namespace jb::testing {

template <typename TWindow>
struct WindowTestConfig {
    using WindowType = TWindow;
    static constexpr std::size_t stage{7};
};

template <>
struct WindowTestConfig<jb::options::Window_Bartlett> {
    using WindowOption = jb::options::Window_Bartlett;
    template <typename SampleCnt, typename Complex>
    using WindowType = jb::windows::BartlettWindow<SampleCnt, Complex>;
    static constexpr std::string_view testFileName = "WinBartlettTest";
};

template <>
struct WindowTestConfig<jb::options::Window_Blackman> {
    using WindowOption = jb::options::Window_Blackman;
    template <typename SampleCnt, typename Complex>
    using WindowType = jb::windows::BlackmanWindow<SampleCnt, Complex>;
    static constexpr std::string_view testFileName = "WinBlackmanTest";
};

template <>
struct WindowTestConfig<jb::options::Window_BlackmanHarris> {
    using WindowOption = jb::options::Window_BlackmanHarris;
    template <typename SampleCnt, typename Complex>
    using WindowType = jb::windows::BlackmanHarrisWindow<SampleCnt, Complex>;
    static constexpr std::string_view testFileName = "WinBlackmanHarrisTest";
};

template <>
struct WindowTestConfig<jb::options::Window_Cosine> {
    using WindowOption = jb::options::Window_Cosine;
    template <typename SampleCnt, typename Complex>
    using WindowType = jb::windows::CosineWindow<SampleCnt, Complex>;
    static constexpr std::string_view testFileName = "WinCosineTest";
};

template <>
struct WindowTestConfig<jb::options::Window_FlatTop> {
    using WindowOption = jb::options::Window_FlatTop;
    template <typename SampleCnt, typename Complex>
    using WindowType = jb::windows::FlatTopWindow<SampleCnt, Complex>;
    static constexpr std::string_view testFileName = "WinFlatTopTest";
};

template <>
struct WindowTestConfig<jb::options::Window_Hamming> {
    using WindowOption = jb::options::Window_Hamming;
    template <typename SampleCnt, typename Complex>
    using WindowType = jb::windows::HammingWindow<SampleCnt, Complex>;
    static constexpr std::string_view testFileName = "WinHammingTest";
};

template <>
struct WindowTestConfig<jb::options::Window_vonHann> {
    using WindowOption = jb::options::Window_vonHann;
    template <typename SampleCnt, typename Complex>
    using WindowType = jb::windows::VonHannWindow<SampleCnt, Complex>;
    static constexpr std::string_view testFileName = "WinvonHannTest";
};

template <>
struct WindowTestConfig<jb::options::Window_Welch> {
    using WindowOption = jb::options::Window_Welch;
    template <typename SampleCnt, typename Complex>
    using WindowType = jb::windows::WelchWindow<SampleCnt, Complex>;
    static constexpr std::string_view testFileName = "WinWelchTest";
};

template <typename TConfig>
class WindowTest
    : public Test
    , public AlgorithmFixture {
public:
    WindowTest() = default;

    void SetUp() override
    {
        mDataSetA.clear();
        mDataSetB.clear();
    }

    ~WindowTest() override = default;

protected:
    static constexpr std::size_t kSampleCnt{128};
};

using WindowTypes = ::testing::Types<
    WindowTestConfig<jb::options::Window_Bartlett>,
    WindowTestConfig<jb::options::Window_Blackman>,
    WindowTestConfig<jb::options::Window_BlackmanHarris>,
    WindowTestConfig<jb::options::Window_Cosine>,
    WindowTestConfig<jb::options::Window_FlatTop>,
    WindowTestConfig<jb::options::Window_Hamming>,
    WindowTestConfig<jb::options::Window_vonHann>,
    WindowTestConfig<jb::options::Window_Welch>>;

TYPED_TEST_SUITE(WindowTest, WindowTypes);

TYPED_TEST(WindowTest, SuccessfullyCalculatesWindowSamples)
{
    ASSERT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, TypeParam::testFileName),
        "win.in", this->mDataSetA,
        "win.out", this->mDataSetB));

    using SampleCnt = std::integral_constant<int, static_cast<int>(TestFixture::kSampleCnt)>;
    using Complex = std::complex<double>;
    using Window = typename TypeParam::template WindowType<SampleCnt, Complex>;

    ASSERT_EQ(this->mDataSetA.size(), TestFixture::kSampleCnt);
    ASSERT_EQ(this->mDataSetB.size(), TestFixture::kSampleCnt);

    Window win;
    win(std::span<Complex, TestFixture::kSampleCnt>(TestFixture::mDataSetA.data(), TestFixture::kSampleCnt));

    AlgorithmResultAnalysis::compareComplexDataSets(TestFixture::mDataSetA, TestFixture::mDataSetB);
}

TYPED_TEST(WindowTest, SuccessfullyCalculateFFTOnWindowedData)
{
    constexpr auto kStage{7};

    std::vector<std::complex<double>> dataSetA;
    std::vector<std::complex<double>> dataSetB;

    jb::AlgorithmFactory<
        kStage, kStage + 1,
        jb::options::Radix_2,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Forward,
        typename TypeParam::WindowOption,
        jb::options::Normalization_Square_Root,
        std::complex<double>>
        fftFactory;
    jb::AlgorithmFactory<
        kStage, kStage + 1,
        jb::options::Radix_2,
        jb::options::Decimation_In_Time,
        jb::options::Direction_Backward,
        jb::options::Window_None,
        jb::options::Normalization_Square_Root,
        std::complex<double>>
        ifftFactory;

    EXPECT_TRUE(AlgorithmResultAnalysis::initialize(
        std::format("{}/{}.xml", TEST_DATA_DIR, TypeParam::testFileName),
        "win.in", dataSetA,
        "win.out", dataSetB));

    auto fft = fftFactory.getAlgorithm(kStage);
    auto ifft = ifftFactory.getAlgorithm(kStage);

    (*fft)(dataSetA);
    (*ifft)(dataSetA);
    AlgorithmResultAnalysis::compareComplexDataSets(dataSetA, dataSetB);
}

} // namespace jb::testing