#ifndef INSPECT_MAINCOMPONENT_HPP_
#define INSPECT_MAINCOMPONENT_HPP_

// CMake builds don't use an AppConfig.h, so it's safe to include juce module headers
// directly. If you need to remain compatible with Projucer-generated builds, and
// have called `juce_generate_juce_header(<thisTarget>)` in your CMakeLists.txt,
// you could `#include <JuceHeader.h>` here instead, to make all your module headers visible.
#include <juce_gui_extra/juce_gui_extra.h>

#include "AudioOutput.hpp"
#include "FilePlayerProcessor.hpp"
#include "SpectrumProcessor.hpp"

namespace inspect {

/**
 * Main component that manages the audio output and file player processor.
 */
class MainComponent final : public juce::Component {
public:
    explicit MainComponent(const juce::File& audioFile = {});

    void paint(juce::Graphics&) override;
    void resized() override;

private:
    AudioOutput mAudioOutput;
    AudioOutput::Node::Ptr mFilePlayerNode; // owns the processor below
    FilePlayerProcessor* mpFilePlayer = nullptr;
    juce::String mStatus;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MainComponent)
};

} // namespace inspect

#endif // INSPECT_MAINCOMPONENT_HPP_
