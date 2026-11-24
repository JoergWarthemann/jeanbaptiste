#include "MainComponent.hpp"

namespace inspect {

MainComponent::MainComponent(const juce::File& audioFile)
{
    auto player = std::make_unique<FilePlayerProcessor>();
    mpFilePlayer = player.get();
    mFilePlayerNode = mAudioOutput.addProcessor(std::move(player));
    auto spectrumProcessor = std::make_unique<SpectrumProcessor>();
    mpSpectrumProcessor = spectrumProcessor.get();
    auto fftNode = mAudioOutput.addProcessor(std::move(spectrumProcessor));
    mAudioOutput.connectChain({mFilePlayerNode, fftNode});

    mPitchKnob.setSliderStyle(juce::Slider::RotaryHorizontalVerticalDrag);
    mPitchKnob.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 80, 20);
    mPitchKnob.setRange(-12.0, 12.0, 0.01);
    mPitchKnob.setDoubleClickReturnValue(true, 0.0);
    mPitchKnob.setTextValueSuffix(" st");
    mPitchKnob.setValue(0.0, juce::dontSendNotification);
    mPitchKnob.onValueChange = [this] {
        mpSpectrumProcessor->setPitchShift(static_cast<float>(mPitchKnob.getValue()));
    };
    addAndMakeVisible(mPitchKnob);

    mPitchLabel.setText("Pitch", juce::dontSendNotification);
    mPitchLabel.setJustificationType(juce::Justification::centred);
    mPitchLabel.attachToComponent(&mPitchKnob, false);
    addAndMakeVisible(mPitchLabel);

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
    g.drawText(mStatus, getLocalBounds().removeFromTop(30), juce::Justification::centred, true);
}

void MainComponent::resized()
{
    mPitchKnob.setBounds(getLocalBounds().withSizeKeepingCentre(120, 140));
}

} // namespace inspect
