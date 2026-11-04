#ifndef INSPECT_SPECTRUMPROCESSOR_HPP_
#define INSPECT_SPECTRUMPROCESSOR_HPP_

#include <array>
#include <complex>
#include <memory>
#include <span>
#include <type_traits>
#include <vector>

#include <juce_audio_utils/juce_audio_utils.h>

#include "IExecutableAlgorithm.hpp"
#include "SpectrumQueue.hpp"
#include "tools/OverlapAdd.hpp"
#include "tools/RingBuffer.hpp"

namespace inspect {

/**
 * Spectral processing node (STFT with overlap-add). The channel mix is analysed every spectrum::kHopSize samples
 * with a real FFT (jean-baptiste-lib) and a sine window, the spectrum is handed to manipulateSpectrum(), transformed
 * back with an inverse FFT and overlap-added (jb::tools::OverlapAdd). The analysis window is applied by the forward transform, the synthesis
 * window by the inverse transform (jb::options::Window_Before / Window_After). The result replaces the audio in the
 * buffer (identical on all channels). Sine analysis and synthesis windows with 50 % overlap reconstruct the input
 * perfectly as long as manipulateSpectrum() leaves the spectrum untouched.
 *
 * The output is delayed by spectrum::kFftSize samples, which is reported as latency.
 * processBlock never allocates or locks. The manipulated spectrum is pushed to a lock-free queue for the GUI.
 * If the GUI doesn't keep up and the queue is full, frames are dropped.
 */
class SpectrumProcessor : public juce::AudioProcessor {
public:
    SpectrumProcessor();
    ~SpectrumProcessor() override;

    /**
     * Gets the queue which receives the magnitudes of the manipulated spectrum.
     * The queue is single producer (this processor) / single consumer, so only one reader may pop from it.
     * Shared ownership keeps it valid for a display that outlives the node.
     */
    std::shared_ptr<spectrum::Queue> getQueue() const noexcept { return mQueue; }

    void prepareToPlay(double sampleRate, int maximumExpectedSamplesPerBlock) override;
    void releaseResources() override {}
    using AudioProcessor::processBlock;
    void processBlock(juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    /**
     * Checks if the given bus layout is supported by the processor.
     *
     * \param layouts The bus layout to check.
     * \return True for equal mono or stereo input and output, false otherwise.
     */
    bool isBusesLayoutSupported(const BusesLayout& layouts) const override;

    const juce::String getName() const override { return "Spectrum Processor"; }
    double getTailLengthSeconds() const override { return 0.0; }
    bool acceptsMidi() const override { return false; }
    bool producesMidi() const override { return false; }

    juce::AudioProcessorEditor* createEditor() override { return nullptr; }
    bool hasEditor() const override { return false; }

    int getNumPrograms() override { return 1; }
    int getCurrentProgram() override { return 0; }
    void setCurrentProgram(int) override {}
    const juce::String getProgramName(int) override { return {}; }
    void changeProgramName(int, const juce::String&) override {}

    void getStateInformation(juce::MemoryBlock&) override {}
    void setStateInformation(const void*, int) override {}

protected:
    using Complex = std::complex<double>;

    /**
     * Hook for changing the spectrum of one frame. Runs on the audio thread, so it must not allocate or lock.
     * The default implementation leaves the spectrum untouched.
     *
     * Only the non-redundant half is passed: bin k is at k * sampleRate / spectrum::kFftSize Hz and the scaling is
     * that of an unnormalized FFT of the windowed frame. The negative frequencies are implied by conjugate symmetry,
     * and the imaginary parts of DC (bin 0) and Nyquist (last bin) are ignored by the inverse transform.
     *
     * \param bins The spectrum bins 0 ... spectrum::kNumBins - 1.
     */
    virtual void manipulateSpectrum(std::span<Complex, spectrum::kNumBins> bins) noexcept;

private:
    static_assert(2 * spectrum::kHopSize == spectrum::kFftSize, "Sine window overlap-add requires 50 % overlap.");

    /**
     * Takes the next output sample and feeds one input sample.
     * \return The output sample, which is the input from kFftSize samples ago after spectral processing.
     */
    float processSample(float input) noexcept;
    void processLatestFrame() noexcept;
    void publishSpectrum() noexcept;

    std::shared_ptr<spectrum::Queue> mQueue = std::make_shared<spectrum::Queue>();
    std::unique_ptr<jb::IExecutableAlgorithm<Complex>> mForward;
    std::unique_ptr<jb::IExecutableAlgorithm<Complex>> mInverse;

    jb::tools::RingBuffer<Complex> mRing; // the most recent kFftSize samples, sized in prepareToPlay
    std::vector<Complex> mFftData; // in-place FFT buffer, sized in prepareToPlay
    jb::tools::OverlapAdd<
        std::integral_constant<int, spectrum::kFftSize>,
        std::integral_constant<int, spectrum::kHopSize>,
        Complex>
        mOverlapAdd;
    std::array<double, spectrum::kHopSize> mPending{}; // the finished output samples of the last frame
    double mMagnitudeScale{1.0}; // makes a full-scale sine read 1.0 at its bin
    spectrum::Frame mFrame;
    std::size_t mPendingPosition{0}; // the next sample of mPending to output
    std::size_t mSamplesSinceLastFrame{0};

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SpectrumProcessor)
};

} // namespace inspect

#endif // INSPECT_SPECTRUMPROCESSOR_HPP_
