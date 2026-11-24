#ifndef INSPECT_PITCHSHIFTER_HPP_
#define INSPECT_PITCHSHIFTER_HPP_

#include <complex>
#include <cstddef>
#include <span>
#include <vector>

namespace inspect {

/**
 * Phase vocoder pitch shifter working on the spectra of consecutive STFT frames of one channel.
 *
 * Each frame, the true frequency of every bin is estimated from the phase advance since the previous frame. The
 * spectral peaks are moved to the bins at their shifted frequencies, each together with the bins around it, and
 * the phase of every peak is rebuilt by accumulating its shifted frequency while the surrounding bins keep their
 * phase relation to it (identity phase locking after Laroche and Dolson). This keeps the phases of consecutive
 * frames consistent, so the overlap-add of the synthesis doesn't smear.
 *
 * The instance holds the phase state of exactly one channel and has to see every frame of that channel in order.
 * prepare() allocates, all other methods never allocate or lock.
 */
class PitchShifter final {
public:
    using Complex = std::complex<double>;

    /**
     * Sizes the state for the given frame geometry and clears it. Not for the audio thread.
     *
     * \param fftSize The length of the frames. The spectra have fftSize / 2 + 1 bins.
     * \param hopSize The number of samples the frames advance by.
     */
    void prepare(std::size_t fftSize, std::size_t hopSize);

    /**
     * Forgets the phase history, so the next frame starts a new stream.
     */
    void reset() noexcept;

    /**
     * Shifts all frequencies of a frame by the given ratio, in place.
     * Content shifted above Nyquist is discarded, bins which receive nothing become zero.
     *
     * \param bins The non-redundant half of the spectrum (DC ... Nyquist).
     * \param ratio The frequency ratio, above 1 shifts up and below 1 shifts down.
     */
    void process(std::span<Complex> bins, double ratio) noexcept;

    /**
     * Leaves the frame untouched, but keeps the phase state in sync with it.
     * Call this for frames which aren't shifted, so that shifting can resume without a glitch.
     *
     * \param bins The non-redundant half of the spectrum (DC ... Nyquist).
     */
    void passThrough(std::span<const Complex> bins) noexcept;

private:
    std::size_t mNumBins{0};
    double mExpectedPhaseChangePerHop{0.0}; // The expected phase change per hop of a sinusoid exactly at bin centre.
    bool mHasHistory{false};
    std::vector<double> mPreviousPhases; // the phase of every bin in the previous frame
    std::vector<double> mSynthesisPhase; // the accumulated output phase of every bin

    // Scratch buffers of one frame.
    std::vector<double> mCurrentMagnitudes;
    std::vector<double> mCurrentPhases;
    std::vector<double> mTrueFrequency; // in bins
    std::vector<double> mOutMagnitude;
    std::vector<double> mOutPhase;
    std::vector<std::size_t> mPeaks;
};

} // namespace inspect

#endif // INSPECT_PITCHSHIFTER_HPP_
