#include "FilePlayerProcessor.hpp"

namespace inspect {

FilePlayerProcessor::FilePlayerProcessor()
    : AudioProcessor(BusesProperties().withOutput("Output", juce::AudioChannelSet::stereo(), true))
{
    mFormatManager.registerBasicFormats();
    mReadAheadThread.startThread();
}

FilePlayerProcessor::~FilePlayerProcessor()
{
    mTransport.stop();
    mTransport.setSource(nullptr);
    mReadAheadThread.stopThread(2000);
}

bool FilePlayerProcessor::loadFile(const juce::File& file)
{
    std::unique_ptr<juce::AudioFormatReader> reader(mFormatManager.createReaderFor(file));
    if (reader == nullptr) {
        return false;
    }

    const double fileSampleRate = reader->sampleRate;
    auto newSource = std::make_unique<juce::AudioFormatReaderSource>(reader.release(), true);
    newSource->setLooping(mLooping);

    // Detach the old source from the transport before it is destroyed.
    mTransport.stop();
    mTransport.setSource(nullptr);
    mReaderSource = std::move(newSource);
    mTransport.setSource(mReaderSource.get(), 32768, &mReadAheadThread, fileSampleRate);

    return true;
}

void FilePlayerProcessor::setLooping(bool shouldLoop)
{
    mLooping = shouldLoop;
    if (mReaderSource != nullptr) {
        mReaderSource->setLooping(shouldLoop);
    }
}

void FilePlayerProcessor::prepareToPlay(double sampleRate, int maximumExpectedSamplesPerBlock)
{
    mTransport.prepareToPlay(maximumExpectedSamplesPerBlock, sampleRate);
}

void FilePlayerProcessor::releaseResources()
{
    mTransport.releaseResources();
}

void FilePlayerProcessor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer&)
{
    juce::ScopedNoDenormals noDenormals;

    // Clears the buffer itself when no source is set or the transport isn't playing.
    juce::AudioSourceChannelInfo info(&buffer, 0, buffer.getNumSamples());
    mTransport.getNextAudioBlock(info);
}

bool FilePlayerProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
    const auto& out = layouts.getMainOutputChannelSet();
    return layouts.getMainInputChannelSet().isDisabled()
        && (out == juce::AudioChannelSet::mono() || out == juce::AudioChannelSet::stereo());
}

} // namespace inspect
