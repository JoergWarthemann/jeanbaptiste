#ifndef INSPECT_FILEPLAYERPROCESSOR_HPP_
#define INSPECT_FILEPLAYERPROCESSOR_HPP_

#include <juce_audio_utils/juce_audio_utils.h>

namespace inspect {

/**
 * Audio processor that plays back an audio file with read-ahead buffering and resampling.
 * Reading happens on a TimeSliceThread (via AudioTransportSource read-ahead), so processBlock
 * never touches the disk. The file is resampled to the device sample rate by the transport source.
 */
class FilePlayerProcessor final : public juce::AudioProcessor {
public:
    FilePlayerProcessor();
    ~FilePlayerProcessor() override;

    /**
     * Loads an audio file for playback.
     *
     * \param file The audio file to load.
     * \return False (and keeps the previous file) if the file can't be read.
     */
    bool loadFile(const juce::File& file);
    bool hasFile() const noexcept { return mReaderSource != nullptr; }

    void start() { mTransport.start(); }
    void stop() { mTransport.stop(); }
    bool isPlaying() const { return mTransport.isPlaying(); }
    void setLooping(bool shouldLoop);
    void setPositionInSeconds(double seconds) { mTransport.setPosition(seconds); }
    double getPositionInSeconds() const { return mTransport.getCurrentPosition(); }
    double getLengthInSeconds() const { return mTransport.getLengthInSeconds(); }

    /**
     * Prepares the processor for playback, allocating any necessary resources.
     *
     * \param sampleRate The sample rate of the audio device.
     * \param maximumExpectedSamplesPerBlock The maximum number of samples per block that will be passed to processBlock.
     */
    void prepareToPlay(double sampleRate, int maximumExpectedSamplesPerBlock) override;

    void releaseResources() override;

    /**
     * Processes a block of audio data.
     *
     * \param buffer The audio buffer to fill with audio data.
     * \param midiMessages The MIDI buffer containing any incoming MIDI messages.
     */
    void processBlock(juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    // Keep the base class's double-precision overload visible.
    using juce::AudioProcessor::processBlock;

    /**
     * Checks if the given bus layout is supported by the processor.
     *
     * \param layouts The bus layout to check.
     * \return True if the layout is supported, false otherwise.
     */
    bool isBusesLayoutSupported(const BusesLayout& layouts) const override;

    const juce::String getName() const override { return "File Player"; }
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
    // Members are destroyed in reverse order: the transport must go before its
    // source and before the read-ahead thread.
    juce::TimeSliceThread mReadAheadThread { "FilePlayer read-ahead" };
    juce::AudioFormatManager mFormatManager;
    std::unique_ptr<juce::AudioFormatReaderSource> mReaderSource;
    juce::AudioTransportSource mTransport;
    bool mLooping = false;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(FilePlayerProcessor)
};

} // namespace inspect

#endif // INSPECT_FILEPLAYERPROCESSOR_HPP_
