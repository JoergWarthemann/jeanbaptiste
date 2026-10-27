#include "MainComponent.hpp"

namespace inspect {

MainComponent::MainComponent(const juce::File& audioFile)
{
    auto player = std::make_unique<FilePlayerProcessor>();
    mpFilePlayer = player.get();
    mFilePlayerNode = mAudioOutput.addProcessor(std::move(player));
    auto fftNode = mAudioOutput.addProcessor(std::make_unique<SpectrumProcessor>());
    mAudioOutput.connectChain({mFilePlayerNode, fftNode});

    if (audioFile == juce::File()) {
        mStatus = "No audio file given.";
    }
    else if (mpFilePlayer->loadFile(audioFile)) {
        mpFilePlayer->start();
        mpFilePlayer->setLooping(true);
        mStatus = "Playing " + audioFile.getFileName();
    }
    else {
        mStatus = "Could not read " + audioFile.getFullPathName();
    }

    if (mAudioOutput.getDeviceError().isNotEmpty())
        mStatus << " (audio device error: " << mAudioOutput.getDeviceError() << ")";

    setSize(600, 400);
}

//==============================================================================
void MainComponent::paint(juce::Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll(getLookAndFeel().findColour(juce::ResizableWindow::backgroundColourId));

    g.setFont(juce::FontOptions(16.0f));
    g.setColour(juce::Colours::white);
    g.drawText(mStatus, getLocalBounds(), juce::Justification::centred, true);
}

void MainComponent::resized()
{
    // This is called when the MainComponent is resized.
    // If you add any child components, this is where you should
    // update their positions.
}

} // namespace inspect
