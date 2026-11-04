#include <cmath>
#include <complex>
#include <numbers>
#include <span>
#include <vector>

#include <gtest/gtest.h>

#include "AlgorithmFactory.hpp"
#include "Options.hpp"
#include "tools/OverlapAdd.hpp"

namespace jb::testing {

namespace {

using Complex = std::complex<double>;

constexpr double kTolerance = 1e-9;

std::vector<Complex> createSignal(const std::size_t size)
{
    std::vector<Complex> data(size);
    for (auto idx = 0U; idx < size; ++idx) {
        data[idx] = Complex(std::sin(0.37 * static_cast<double>(idx)) + 0.25 * static_cast<double>(idx % 3), 0.5);
    }

    return data;
}

// Slices the signal into frames, overlap-adds them and returns the concatenated finished samples.
template <std::size_t FrameSize, std::size_t HopSize>
std::vector<double> overlapAddFrames(const std::vector<Complex>& frames, const std::size_t numFrames,
    tools::OverlapAdd<std::integral_constant<int, FrameSize>, std::integral_constant<int, HopSize>, Complex>& overlapAdd)
{
    std::vector<double> result;
    for (std::size_t frame = 0; frame < numFrames; ++frame) {
        std::array<double, HopSize> out{};
        overlapAdd(std::span<const Complex, FrameSize>(frames.data() + frame * FrameSize, FrameSize), out);
        result.insert(result.end(), out.begin(), out.end());
    }

    return result;
}

} // namespace

TEST(OverlapAddTest, ReturnsOnlyTheFirstFrameBeforeAnyOverlapExists)
{
    constexpr auto kSignalSize{8};
    constexpr auto kHopSize{4};
    tools::OverlapAdd<
        std::integral_constant<int, kSignalSize>,
        std::integral_constant<int, kHopSize>,
        Complex>
        overlapAdd;
    const auto frame = createSignal(kSignalSize);
    std::array<double, kHopSize> out{};

    overlapAdd(std::span<const Complex, kSignalSize>(frame), out);

    for (auto idx = 0U; idx < out.size(); ++idx) {
        EXPECT_DOUBLE_EQ(out[idx], frame[idx].real());
    }
}

TEST(OverlapAddTest, AddsTheTailOfThePreviousFrameWith50PercentOverlap)
{
    constexpr auto kSignalSize{8};
    constexpr auto kHopSize{4};
    tools::OverlapAdd<
        std::integral_constant<int, kSignalSize>,
        std::integral_constant<int, kHopSize>,
        Complex>
        overlapAdd;
    const auto first = createSignal(kSignalSize);

    std::vector<Complex> second(kSignalSize);
    for (auto idx = 0U; idx < second.size(); ++idx) {
        second[idx] = Complex(10.0 + static_cast<double>(idx), 7.0);
    }

    std::array<double, kHopSize> out{};
    overlapAdd(std::span<const Complex, kSignalSize>(first), out);
    overlapAdd(std::span<const Complex, kSignalSize>(second), out);

    for (auto idx = 0U; idx < out.size(); ++idx) {
        EXPECT_NEAR(out[idx], first[idx + 4].real() + second[idx].real(), kTolerance) << "sample " << idx;
    }
}

TEST(OverlapAddTest, IgnoresTheImaginaryPart)
{
    constexpr auto kSignalSize{8};
    constexpr auto kHopSize{4};
    tools::OverlapAdd<
        std::integral_constant<int, kSignalSize>,
        std::integral_constant<int, kHopSize>,
        Complex>
        overlapAdd;
    const std::vector<Complex> frame(kSignalSize, Complex(1.0, 1000.0));
    std::array<double, kHopSize> out{};

    overlapAdd(std::span<const Complex, kSignalSize>(frame), out);
    overlapAdd(std::span<const Complex, kSignalSize>(frame), out);

    for (const auto sample : out) {
        EXPECT_DOUBLE_EQ(sample, 2.0);
    }
}

TEST(OverlapAddTest, SumsAllOverlappingFramesWith75PercentOverlap)
{
    constexpr std::size_t kFrameSize{8};
    constexpr std::size_t kHopSize{2};
    tools::OverlapAdd<
        std::integral_constant<int, kFrameSize>,
        std::integral_constant<int, kHopSize>,
        Complex>
        overlapAdd;

    constexpr std::size_t kNumFrames{8};
    // Constant frames of one, so every output sample equals the number of frames which cover it.
    const std::vector<Complex> frames(kFrameSize * kNumFrames, Complex(1.0, 0.0));
    const auto result = overlapAddFrames<kFrameSize, kHopSize>(frames, kNumFrames, overlapAdd);

    ASSERT_EQ(result.size(), kNumFrames * kHopSize);
    EXPECT_DOUBLE_EQ(result[0], 1.0);
    EXPECT_DOUBLE_EQ(result[2], 2.0);
    EXPECT_DOUBLE_EQ(result[4], 3.0);
    for (auto idx = 6U; idx < result.size(); ++idx) {
        EXPECT_DOUBLE_EQ(result[idx], 4.0) << "sample " << idx;
    }
}

TEST(OverlapAddTest, SimplyConcatenatesFramesWithoutOverlap)
{
    constexpr std::size_t kFrameSize{8};
    constexpr std::size_t kHopSize{kFrameSize};
    tools::OverlapAdd<
        std::integral_constant<int, kFrameSize>,
        std::integral_constant<int, kHopSize>,
        Complex>
        overlapAdd;
    constexpr std::size_t kNumFrames{3};
    const auto frames = createSignal(kFrameSize * kNumFrames);

    const auto result = overlapAddFrames<kFrameSize, kHopSize>(frames, kNumFrames, overlapAdd);

    ASSERT_EQ(result.size(), frames.size());
    for (auto idx = 0U; idx < result.size(); ++idx) {
        EXPECT_DOUBLE_EQ(result[idx], frames[idx].real());
    }
}

TEST(OverlapAddTest, ResetDiscardsTheUnfinishedTail)
{
    constexpr auto kSignalSize{8};
    constexpr auto kHopSize{4};
    tools::OverlapAdd<
        std::integral_constant<int, kSignalSize>,
        std::integral_constant<int, kHopSize>,
        Complex>
        overlapAdd;
    const std::vector<Complex> frame(kSignalSize, Complex(1.0, 0.0));
    std::array<double, kHopSize> out{};

    overlapAdd(std::span<const Complex, 8>(frame), out);
    overlapAdd.reset();
    overlapAdd(std::span<const Complex, 8>(frame), out);

    for (const auto sample : out) {
        EXPECT_DOUBLE_EQ(sample, 1.0);
    }
}

// Sine analysis and synthesis windows with 50 % overlap reconstruct the signal one hop later.
TEST(OverlapAddTest, ReconstructsSignalFromWindowedFftFrames)
{
    constexpr std::size_t kStage{6};
    constexpr std::size_t kFrameSize{1 << kStage};
    constexpr std::size_t kHopSize{kFrameSize / 2};
    constexpr std::size_t kNumFrames{12};

    auto analysis =
        AlgorithmFactory<
            kStage,
            kStage + 1,
            options::Radix_2,
            options::Decimation_In_Time,
            options::Direction_Forward,
            options::Window_Cosine,
            options::Normalization_No,
            Complex,
            options::Data_Real,
            options::Window_Before>{}
            .getAlgorithm(kStage);
    auto synthesis =
        AlgorithmFactory<
            kStage,
            kStage + 1,
            options::Radix_2,
            options::Decimation_In_Time,
            options::Direction_Backward,
            options::Window_Cosine,
            options::Normalization_Division_By_Length,
            Complex,
            options::Data_Real,
            options::Window_After>{}
            .getAlgorithm(kStage);
    tools::OverlapAdd<
        std::integral_constant<int, kFrameSize>,
        std::integral_constant<int, kHopSize>,
        Complex>
        overlapAdd;

    const auto signal = createSignal(kHopSize * (kNumFrames + 1));

    std::vector<double> output;
    for (std::size_t frame = 0; frame < kNumFrames; ++frame) {
        std::vector<Complex> data(
            signal.begin() + static_cast<std::ptrdiff_t>(frame * kHopSize),
            signal.begin() + static_cast<std::ptrdiff_t>(frame * kHopSize + kFrameSize));

        for (auto& element : data) {
            element.imag(0.0);
        }

        (*analysis)(data);
        (*synthesis)(data);

        std::array<double, kHopSize> out{};
        overlapAdd(std::span<const Complex, kFrameSize>(data), out);
        output.insert(output.end(), out.begin(), out.end());
    }

    // The first hop is covered by one frame only, from there on output[n] equals the input.
    for (std::size_t idx = kHopSize; idx < output.size(); ++idx) {
        EXPECT_NEAR(output[idx], signal[idx].real(), kTolerance) << "sample " << idx;
    }
}

} // namespace jb::testing
