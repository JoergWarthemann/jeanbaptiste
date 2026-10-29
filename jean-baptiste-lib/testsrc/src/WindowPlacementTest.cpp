#include <cmath>
#include <complex>
#include <numbers>
#include <vector>

#include <gtest/gtest.h>

#include "AlgorithmFactory.hpp"
#include "Options.hpp"
#include "windows/ExecuteWindowOnComplexData.hpp"

namespace jb::testing {

namespace {

using Complex = std::complex<double>;

constexpr double kTolerance = 1e-9;

// The sine window, i.e. jb::windows::CosineWindow, evaluated independently of the library.
double sineWindow(const std::size_t index, const std::size_t size)
{
    return std::sin(std::numbers::pi * static_cast<double>(index) / static_cast<double>(size));
}

std::vector<Complex> createComplexSignal(const std::size_t size, const unsigned seed = 0)
{
    std::vector<Complex> data(size);
    for (auto idx = 0U; idx < size; ++idx) {
        data[idx] = Complex(
            std::sin(0.37 * static_cast<double>(idx + seed)) + 0.25 * static_cast<double>((idx + seed) % 3),
            std::cos(0.91 * static_cast<double>(idx + seed)) - 0.1 * static_cast<double>((idx + seed) % 5));
    }

    return data;
}

std::vector<Complex> createRealSignal(const std::size_t size, const unsigned seed = 0)
{
    auto data = createComplexSignal(size, seed);
    for (auto& element : data) {
        element.imag(0.0);
    }

    return data;
}

void expectNear(const std::vector<Complex>& actual, const std::vector<Complex>& expected)
{
    ASSERT_EQ(actual.size(), expected.size());
    for (auto idx = 0U; idx < actual.size(); ++idx) {
        EXPECT_NEAR(actual[idx].real(), expected[idx].real(), kTolerance) << "real part of sample " << idx;
        EXPECT_NEAR(actual[idx].imag(), expected[idx].imag(), kTolerance) << "imaginary part of sample " << idx;
    }
}

template <typename Radix, std::size_t Stage, typename Data, typename Direction, typename Window, typename Normalization, typename WindowPlacement>
auto createAlgorithm()
{
    return AlgorithmFactory<Stage, Stage + 1, Radix, options::Decimation_In_Time, Direction, Window, Normalization, Complex, Data, WindowPlacement>{}
        .getAlgorithm(Stage);
}

// Window_After on a backward transform must equal the plain backward transform multiplied by the window.
template <typename Radix, std::size_t Stage, typename Data>
void expectWindowAfterMatchesWindowedResult(const std::vector<Complex>& spectrum)
{
    const auto size = spectrum.size();

    auto withWindow = spectrum;
    createAlgorithm<Radix, Stage, Data, options::Direction_Backward, options::Window_Cosine, options::Normalization_Division_By_Length, options::Window_After>()
        ->operator()(withWindow);

    auto expected = spectrum;
    createAlgorithm<Radix, Stage, Data, options::Direction_Backward, options::Window_None, options::Normalization_Division_By_Length, options::Window_Before>()
        ->operator()(expected);
    for (auto idx = 0U; idx < size; ++idx) {
        expected[idx] *= sineWindow(idx, size);
    }

    expectNear(withWindow, expected);
}

// Window_Before on a forward transform must equal the plain forward transform of the windowed input.
template <typename Radix, std::size_t Stage, typename Data, typename WindowPlacement>
void expectWindowBeforeMatchesWindowedInput(const std::vector<Complex>& signal)
{
    const auto size = signal.size();

    auto withWindow = signal;
    createAlgorithm<Radix, Stage, Data, options::Direction_Forward, options::Window_Cosine, options::Normalization_No, WindowPlacement>()
        ->operator()(withWindow);

    auto expected = signal;
    for (auto idx = 0U; idx < size; ++idx) {
        expected[idx] *= sineWindow(idx, size);
    }
    createAlgorithm<Radix, Stage, Data, options::Direction_Forward, options::Window_None, options::Normalization_No, options::Window_Before>()
        ->operator()(expected);

    expectNear(withWindow, expected);
}

// Sine analysis window + sine synthesis window with 50 % overlap reconstruct the signal where two frames overlap.
template <typename Radix, std::size_t Stage, typename Data>
void expectOverlapAddReconstructsSignal(const std::vector<Complex>& signal, const std::size_t size)
{
    const auto hop = size / 2;
    auto analysis = createAlgorithm<Radix, Stage, Data, options::Direction_Forward, options::Window_Cosine, options::Normalization_No, options::Window_Before>();
    auto synthesis = createAlgorithm<Radix, Stage, Data, options::Direction_Backward, options::Window_Cosine, options::Normalization_Division_By_Length, options::Window_After>();

    std::vector<Complex> result(signal.size());
    for (std::size_t start = 0; start + size <= signal.size(); start += hop) {
        std::vector<Complex> frame(signal.begin() + static_cast<std::ptrdiff_t>(start), signal.begin() + static_cast<std::ptrdiff_t>(start + size));
        (*analysis)(frame);
        (*synthesis)(frame);
        for (std::size_t idx = 0; idx < size; ++idx) {
            result[start + idx] += frame[idx];
        }
    }

    // The first and last hop samples are covered by only one frame.
    const std::vector<Complex> interiorResult(result.begin() + static_cast<std::ptrdiff_t>(hop), result.end() - static_cast<std::ptrdiff_t>(hop));
    const std::vector<Complex> interiorSignal(signal.begin() + static_cast<std::ptrdiff_t>(hop), signal.end() - static_cast<std::ptrdiff_t>(hop));
    expectNear(interiorResult, interiorSignal);
}

} // namespace

TEST(WindowPlacementTest, WindowFunctorScalesRealAndImaginaryPart)
{
    const Complex scaled = windows::ExecuteWindowOnComplexData<Complex>{}(Complex(2.0, -3.0), 0.5);

    EXPECT_DOUBLE_EQ(scaled.real(), 1.0);
    EXPECT_DOUBLE_EQ(scaled.imag(), -1.5);
}

TEST(WindowPlacementTest, ComplexRadix2BackwardWindowAfterMultipliesResultByWindow)
{
    expectWindowAfterMatchesWindowedResult<options::Radix_2, 6, options::Data_Complex>(createComplexSignal(64));
}

TEST(WindowPlacementTest, ComplexRadix4BackwardWindowAfterMultipliesResultByWindow)
{
    expectWindowAfterMatchesWindowedResult<options::Radix_4, 3, options::Data_Complex>(createComplexSignal(64));
}

TEST(WindowPlacementTest, ComplexRadixSplit24BackwardWindowAfterMultipliesResultByWindow)
{
    expectWindowAfterMatchesWindowedResult<options::Radix_Split_2_4, 6, options::Data_Complex>(createComplexSignal(64));
}

TEST(WindowPlacementTest, RealBackwardWindowAfterMultipliesResultByWindow)
{
    // A valid half spectrum is produced by a forward transform of a real signal.
    auto spectrum = createRealSignal(64);
    createAlgorithm<options::Radix_2, 6, options::Data_Real, options::Direction_Forward, options::Window_None, options::Normalization_No, options::Window_Before>()
        ->operator()(spectrum);

    expectWindowAfterMatchesWindowedResult<options::Radix_2, 6, options::Data_Real>(spectrum);
}

TEST(WindowPlacementTest, ComplexForwardWindowBeforeMultipliesInputByWindow)
{
    expectWindowBeforeMatchesWindowedInput<options::Radix_2, 6, options::Data_Complex, options::Window_Before>(createComplexSignal(64));
}

TEST(WindowPlacementTest, RealForwardWindowBeforeMultipliesInputByWindow)
{
    expectWindowBeforeMatchesWindowedInput<options::Radix_2, 6, options::Data_Real, options::Window_Before>(createRealSignal(64));
}

TEST(WindowPlacementTest, WindowBeforeIsTheDefaultPlacement)
{
    const auto signal = createComplexSignal(64);

    auto byDefault = signal;
    AlgorithmFactory<6, 7, options::Radix_2, options::Decimation_In_Time, options::Direction_Forward, options::Window_vonHann, options::Normalization_No, Complex>{}
        .getAlgorithm(6)
        ->operator()(byDefault);

    auto explicitBefore = signal;
    createAlgorithm<options::Radix_2, 6, options::Data_Complex, options::Direction_Forward, options::Window_vonHann, options::Normalization_No, options::Window_Before>()
        ->operator()(explicitBefore);

    expectNear(byDefault, explicitBefore);
}

TEST(WindowPlacementTest, WindowAfterOnForwardTransformMultipliesSpectrumByWindow)
{
    const auto signal = createComplexSignal(64);

    auto withWindow = signal;
    createAlgorithm<options::Radix_2, 6, options::Data_Complex, options::Direction_Forward, options::Window_Cosine, options::Normalization_No, options::Window_After>()
        ->operator()(withWindow);

    auto expected = signal;
    createAlgorithm<options::Radix_2, 6, options::Data_Complex, options::Direction_Forward, options::Window_None, options::Normalization_No, options::Window_Before>()
        ->operator()(expected);
    for (auto idx = 0U; idx < expected.size(); ++idx) {
        expected[idx] *= sineWindow(idx, expected.size());
    }

    expectNear(withWindow, expected);
}

TEST(WindowPlacementTest, ComplexOverlapAddReconstructsComplexSignal)
{
    expectOverlapAddReconstructsSignal<options::Radix_2, 6, options::Data_Complex>(createComplexSignal(64 * 8), 64);
}

TEST(WindowPlacementTest, ComplexRadix4OverlapAddReconstructsComplexSignal)
{
    expectOverlapAddReconstructsSignal<options::Radix_4, 3, options::Data_Complex>(createComplexSignal(64 * 8), 64);
}

TEST(WindowPlacementTest, RealOverlapAddReconstructsRealSignal)
{
    expectOverlapAddReconstructsSignal<options::Radix_2, 6, options::Data_Real>(createRealSignal(64 * 8), 64);
}

TEST(WindowPlacementTest, RealOverlapAddReconstructsRealSignalWithLargeFrames)
{
    // Frame sizes above 256 samples exercise the bit reversal for larger stages.
    expectOverlapAddReconstructsSignal<options::Radix_2, 11, options::Data_Real>(createRealSignal(2048 * 4), 2048);
}

} // namespace jb::testing
