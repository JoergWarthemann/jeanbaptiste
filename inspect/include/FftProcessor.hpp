#ifndef INSPECT_FFTPROCESSOR_HPP_
#define INSPECT_FFTPROCESSOR_HPP_

#include <array>
#include <complex>
#include <memory>
#include <type_traits>
#include <vector>

#include <juce_audio_utils/juce_audio_utils.h>

#include "IExecutableAlgorithm.hpp"
#include "SpectrumQueue.hpp"

namespace inspect {

/**
 * Pass-through analyser node. The audio is left untouched, while a windowed real FFT (jean-baptiste-lib)
 * of the channel mix is computed every spectrum::kHopSize samples and pushed to a lock-free queue for the GUI.
 * processBlock never allocates or locks. If the GUI doesn't keep up and the queue is full, frames are dropped.
 */
class FftProcessor final : public juce::AudioProcessor {
public:
    FftProcessor();
    ~FftProcessor() override;

    /**
     * Gets the queue which receives the spectrum frames.
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

    const juce::String getName() const override { return "FFT Analyser"; }
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

private:
    using Complex = std::complex<double>;

    void pushSample(float sample) noexcept;
    void analyseLatestFrame() noexcept;

    std::shared_ptr<spectrum::Queue> mQueue = std::make_shared<spectrum::Queue>();
    std::unique_ptr<jb::IExecutableAlgorithm<Complex>> mAlgorithm;

    std::array<float, spectrum::kFftSize> mRing{}; // the most recent kFftSize samples
    std::vector<Complex> mFftData; // in-place FFT buffer, sized in prepareToPlay
    spectrum::Frame mFrame;
    std::size_t mWritePosition = 0;
    std::size_t mSamplesSinceLastFrame = 0;
    std::size_t mNumValidSamples = 0;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(FftProcessor)
};

} // namespace inspect

#endif // INSPECT_FFTPROCESSOR_HPP_
