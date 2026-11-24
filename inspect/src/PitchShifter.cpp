#include "PitchShifter.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numbers>

namespace inspect {

namespace {

constexpr double kTwoPi = 2.0 * std::numbers::pi;
constexpr std::size_t kPeakRadius = 2; // a peak is the strongest bin within this distance

// Wraps a phase to [-pi, pi].
double wrapPhase(double phase) noexcept
{
    return phase - kTwoPi * std::round(phase / kTwoPi);
}

} // namespace

void PitchShifter::prepare(std::size_t fftSize, std::size_t hopSize)
{
    mNumBins = fftSize / 2 + 1;
    // The forward transform of jean-baptiste-lib uses the conjugate of the usual convention, so the phase of a
    // stationary sinusoid decreases from frame to frame.
    mExpectedPhaseChangePerHop = -kTwoPi * static_cast<double>(hopSize) / static_cast<double>(fftSize);
    mPreviousPhases.assign(mNumBins, 0.0);
    mSynthesisPhase.assign(mNumBins, 0.0);
    mCurrentMagnitudes.assign(mNumBins, 0.0);
    mCurrentPhases.assign(mNumBins, 0.0);
    mTrueFrequency.assign(mNumBins, 0.0);
    mOutMagnitude.assign(mNumBins, 0.0);
    mOutPhase.assign(mNumBins, 0.0);
    mPeaks.assign(mNumBins, 0);
    mHasHistory = false;
}

void PitchShifter::reset() noexcept
{
    std::ranges::fill(mPreviousPhases, 0.0);
    std::ranges::fill(mSynthesisPhase, 0.0);
    mHasHistory = false;
}

void PitchShifter::process(std::span<Complex> bins, double ratio) noexcept
{
    assert(bins.size() == mNumBins);

    // Go through all bins and extract magnitude and phase.
    for (std::size_t bin = 0; bin < mNumBins; ++bin) {
        mCurrentMagnitudes[bin] = std::abs(bins[bin]);
        mCurrentPhases[bin] = std::arg(bins[bin]);

        const double binPosition = static_cast<double>(bin);
        const double expectedPhaseAdvance = mExpectedPhaseChangePerHop * binPosition;
        // The deviation from the phase advance of the bin's centre frequency gives the true frequency.
        // It is current phase - previous phase - expected phase change per hop * the bin position.
        const double deviation =
            mHasHistory ? wrapPhase(mCurrentPhases[bin] - mPreviousPhases[bin] - expectedPhaseAdvance) : 0.0;
        mTrueFrequency[bin] = binPosition + deviation / mExpectedPhaseChangePerHop;
        mPreviousPhases[bin] = mCurrentPhases[bin];
    }
    mHasHistory = true;

    // The peaks are the local maxima of the magnitude.
    std::size_t numPeaks = 0;
    for (std::size_t bin = 0; bin < mNumBins; ++bin) {
        const std::size_t first = bin >= kPeakRadius ? bin - kPeakRadius : 0;
        const std::size_t last = std::min(bin + kPeakRadius, mNumBins - 1);
        bool isPeak = mCurrentMagnitudes[bin] > 0.0;
        for (std::size_t other = first; isPeak && other <= last; ++other) {
            isPeak = other == bin || mCurrentMagnitudes[bin] > mCurrentMagnitudes[other];
        }
        if (isPeak) {
            mPeaks[numPeaks++] = bin;
        }
    }

    // Bins without content keep advancing at their own frequency.
    std::ranges::fill(mOutMagnitude, 0.0);
    for (std::size_t bin = 0; bin < mNumBins; ++bin) {
        mOutPhase[bin] = wrapPhase(mSynthesisPhase[bin] + mExpectedPhaseChangePerHop * static_cast<double>(bin));
    }

    // Every peak moves together with the bins up to the middle to its neighbouring peaks, which keeps the phase
    // relations within a sinusoid's main lobe intact. Only the peak's phase is accumulated, the other bins follow it.
    for (std::size_t idx = 0; idx < numPeaks; ++idx) {
        const std::size_t peak = mPeaks[idx];
        const double target = std::round(static_cast<double>(peak) * ratio);
        if (target >= static_cast<double>(mNumBins)) {
            break;
        }

        const auto targetPeak = static_cast<std::size_t>(target);
        const std::size_t regionBegin = idx == 0 ? 0 : (mPeaks[idx - 1] + peak) / 2 + 1;
        const std::size_t regionEnd = idx + 1 == numPeaks ? mNumBins - 1 : (peak + mPeaks[idx + 1]) / 2;

        const double peakPhase =
            wrapPhase(mSynthesisPhase[targetPeak] + mExpectedPhaseChangePerHop * mTrueFrequency[peak] * ratio);

        for (std::size_t bin = regionBegin; bin <= regionEnd; ++bin) {
            // Moves the bin by the same integer offset as the peak. The target is within range if the offset is.
            const std::ptrdiff_t shifted =
                static_cast<std::ptrdiff_t>(bin) + static_cast<std::ptrdiff_t>(targetPeak) - static_cast<std::ptrdiff_t>(peak);
            if (shifted < 0 || shifted >= static_cast<std::ptrdiff_t>(mNumBins)) {
                continue;
            }

            const auto out = static_cast<std::size_t>(shifted);
            if (bin == peak || mCurrentMagnitudes[bin] > mOutMagnitude[out]) {
                mOutMagnitude[out] = mCurrentMagnitudes[bin];
                mOutPhase[out] = wrapPhase(peakPhase + mCurrentPhases[bin] - mCurrentPhases[peak]);
            }
        }
    }

    for (std::size_t bin = 0; bin < mNumBins; ++bin) {
        mSynthesisPhase[bin] = mOutPhase[bin];
        bins[bin] = std::polar(mOutMagnitude[bin], mOutPhase[bin]);
    }

    // The inverse transform ignores the imaginary parts of DC and Nyquist.
    bins.front() = Complex(bins.front().real(), 0.0);
    bins.back() = Complex(bins.back().real(), 0.0);
}

void PitchShifter::passThrough(std::span<const Complex> bins) noexcept
{
    assert(bins.size() == mNumBins);

    for (std::size_t bin = 0; bin < mNumBins; ++bin) {
        mPreviousPhases[bin] = std::arg(bins[bin]);
        mSynthesisPhase[bin] = mPreviousPhases[bin];
    }
    mHasHistory = true;
}

} // namespace inspect
