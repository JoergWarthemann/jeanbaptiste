#ifndef INSPECT_SPECTRUMPROCESSOR_HPP_
#define INSPECT_SPECTRUMPROCESSOR_HPP_

#include <algorithm>
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
 * Spectral processing node (STFT with overlap-add). Each channel is analysed every spectrum::kHopSize samples
 * with a real FFT (jean-baptiste-lib) and a sine window, the spectrum is handed to manipulateSpectrum(), transformed
 * back with an inverse FFT and overlap-added (jb::tools::OverlapAdd). The analysis window is applied by the forward transform, the synthesis
 * window by the inverse transform (jb::options::Window_Before / Window_After). The result replaces the audio in the
 * buffer independently per channel. Sine analysis and synthesis windows with 50 % overlap reconstruct the input
 * perfectly as long as manipulateSpectrum() leaves the spectrum untouched.
 *
 * The output is delayed by spectrum::kFftSize samples, which is reported as latency.
 * processBlock never allocates or locks. The mean channel magnitudes are pushed to a lock-free queue for the GUI.
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
     * Hook for changing the spectrum of one channel's frame. Runs once per channel on the audio thread,
     * so it must not allocate or lock.
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

    struct ChannelState {
        std::unique_ptr<jb::IExecutableAlgorithm<Complex>> forward;
        std::unique_ptr<jb::IExecutableAlgorithm<Complex>> inverse;
        jb::tools::RingBuffer<Complex> ring;
        std::vector<Complex> fftData;
        jb::tools::OverlapAdd<
            std::integral_constant<int, spectrum::kFftSize>,
            std::integral_constant<int, spectrum::kHopSize>,
            Complex>
            overlapAdd;
        std::array<double, spectrum::kHopSize> pending{};
    };

    void exchangeSamples(std::span<float> samples, std::span<const double> output, ChannelState& channel) noexcept;

    void processSampleFrameAndSpectrum(int numChannels) noexcept;
    void publishSpectrum(int numChannels) noexcept;

    std::shared_ptr<spectrum::Queue> mQueue = std::make_shared<spectrum::Queue>();
    std::array<ChannelState, 2> mChannels;
    double mMagnitudeScale{1.0}; // makes a full-scale sine read 1.0 at its bin
    spectrum::Frame mFrame;
    std::size_t mPendingPosition{0}; // the next pending sample to output for each channel
    std::size_t mSamplesSinceLastFrame{0};

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SpectrumProcessor)
};

} // namespace inspect

#endif // INSPECT_SPECTRUMPROCESSOR_HPP_
