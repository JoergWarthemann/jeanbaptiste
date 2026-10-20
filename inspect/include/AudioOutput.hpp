#ifndef INSPECT_AUDIOOUTPUT_HPP_
#define INSPECT_AUDIOOUTPUT_HPP_

#include <juce_audio_utils/juce_audio_utils.h>

namespace inspect {

/**
 * Owns the audio device and an AudioProcessorGraph, and manages playback through an AudioProcessorPlayer:
 *     AudioDeviceManager --> AudioProcessorPlayer --> AudioProcessorGraph
 * Processors are added as graph nodes and routed in series with connectChain(), ending in the graph's audio
 * output node. All methods are for the message thread.
 */
class AudioOutput final {
public:
    using Node = juce::AudioProcessorGraph::Node;

    AudioOutput();
    ~AudioOutput();

    /**
     * Returns an empty string if the default audio device was opened successfully,
     * otherwise contains an error message.
     */
    const juce::String& getDeviceError() const noexcept { return mDeviceError; }

    /**
     * Returns a reference to the AudioDeviceManager used by this AudioOutput instance.
     *
     * \return Reference to the AudioDeviceManager.
     */
    juce::AudioDeviceManager& getDeviceManager() noexcept { return mDeviceManager; }

    /**
     * Returns a reference to the AudioProcessorGraph used by this AudioOutput instance.
     *
     * \return Reference to the AudioProcessorGraph.
     */
    juce::AudioProcessorGraph& getGraph() noexcept { return mGraph; }

    /**
     * Adds a processor to the graph as a new node.
     * The graph takes ownership of the processor. Keep the returned node (or a raw pointer to the
     * processor) to reach it later; the processor lives as long as the node.
     *
     * \param processor The processor to add.
     * \return A pointer to the newly added node.
     */
    Node::Ptr addProcessor(std::unique_ptr<juce::AudioProcessor> processor);
    void removeProcessor(const Node::Ptr& node);

    /**
     * Replaces all connections of the given nodes by: chain[0] -> chain[1] -> ... -> output.
     *
     * \param chain The sequence of nodes to connect in series.
     * \return False if any connection was rejected by the graph.
     */
    bool connectChain(const std::vector<Node::Ptr>& chain);

private:
    /**
     * Connects the source node to the destination node in the graph.
     *
     * \param source The source node.
     * \param destination The destination node.
     * \return True if the connection was successful, false otherwise.
     */
    bool connectNodes(const Node::Ptr& source, const Node::Ptr& destination);

    // Members are destroyed in reverse order: the player has to be detached from the
    // device manager and the graph before either of them goes away.
    juce::AudioDeviceManager mDeviceManager;
    juce::AudioProcessorGraph mGraph;
    juce::AudioProcessorPlayer mPlayer;
    Node::Ptr mOutputNode;
    juce::String mDeviceError;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AudioOutput)
};

} // namespace inspect

#endif // INSPECT_AUDIOOUTPUT_HPP_
