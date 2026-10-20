#include "AudioOutput.hpp"

namespace inspect {

AudioOutput::AudioOutput()
{
    // Stereo default so the graph topology is valid even without a working device;
    // the player reconfigures the graph for the real device on start.
    mGraph.setPlayConfigDetails(0, 2, 44100.0, 512);

    using IO = juce::AudioProcessorGraph::AudioGraphIOProcessor;
    mOutputNode = mGraph.addNode(std::make_unique<IO>(IO::audioOutputNode));

    mDeviceError = mDeviceManager.initialiseWithDefaultDevices(0, 2);
    if (mDeviceError.isNotEmpty()) {
        juce::Logger::writeToLog("AudioOutput: could not open audio device: " + mDeviceError);
    }

    mPlayer.setProcessor(&mGraph);
    mDeviceManager.addAudioCallback(&mPlayer);
}

AudioOutput::~AudioOutput()
{
    mDeviceManager.removeAudioCallback(&mPlayer);
    mPlayer.setProcessor(nullptr);
}

AudioOutput::Node::Ptr AudioOutput::addProcessor(std::unique_ptr<juce::AudioProcessor> processor)
{
    return mGraph.addNode(std::move(processor));
}

void AudioOutput::removeProcessor(const Node::Ptr& node)
{
    if (node != nullptr && node != mOutputNode) {
        mGraph.removeNode(node->nodeID);
    }
}

bool AudioOutput::connectChain(const std::vector<Node::Ptr>& chain)
{
    for (const auto& node : chain) {
        mGraph.disconnectNode(node->nodeID);
    }

    mGraph.disconnectNode(mOutputNode->nodeID);

    bool allConnected = true;
    for (std::size_t i = 0; i < chain.size(); ++i) {
        allConnected &= connectNodes(chain[i], i + 1 < chain.size() ? chain[i + 1] : mOutputNode);
    }

    return allConnected;
}

bool AudioOutput::connectNodes(const Node::Ptr& source, const Node::Ptr& destination)
{
    const int numChannels = std::min(source->getProcessor()->getTotalNumOutputChannels(),
        destination->getProcessor()->getTotalNumInputChannels());

    bool allConnected = numChannels > 0;
    for (int channel = 0; channel < numChannels; ++channel) {
        allConnected &= mGraph.addConnection({{source->nodeID, channel}, {destination->nodeID, channel}});
    }

    return allConnected;
}

} // namespace inspect
