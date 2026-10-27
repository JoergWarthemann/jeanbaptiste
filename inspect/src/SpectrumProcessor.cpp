#include "SpectrumProcessor.hpp"

#include "AlgorithmFactory.hpp"
#include "Options.hpp"

namespace inspect {

namespace {

// Real FFT, radix-2 DIT, forward, von Hann window, no normalization (scaled below).
using FftFactory = jb::AlgorithmFactory<
    spectrum::kFftStage, spectrum::kFftStage + 1,
    jb::options::Radix_2,
    jb::options::Decimation_In_Time,
    jb::options::Direction_Forward,
    jb::options::Window_vonHann,
    jb::options::Normalization_No,
    std::complex<double>,
    jb::options::Data_Real>;

constexpr double kWindowCoherentGain = 0.5; // mean of the von Hann window

// A sine of amplitude A at bin k gives |X[k]| = A * N / 2 * kWindowCoherentGain.
constexpr double kBinScale = 2.0 / (spectrum::kFftSize * kWindowCoherentGain);

} // namespace

SpectrumProcessor::SpectrumProcessor()
    : AudioProcessor(BusesProperties()
                         .withInput("Input", juce::AudioChannelSet::stereo(), true)
                         .withOutput("Output", juce::AudioChannelSet::stereo(), true))
{
}

SpectrumProcessor::~SpectrumProcessor() = default;

//==============================================================================
void SpectrumProcessor::prepareToPlay(double, int)
{
    mAlgorithm = FftFactory{}.getAlgorithm(spectrum::kFftStage);
    mFftData.assign(spectrum::kFftSize, Complex{});

    mRing.reset(spectrum::kFftSize);
    mSamplesSinceLastFrame = 0;
}

void SpectrumProcessor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer&)
{
    juce::ScopedNoDenormals noDenormals;

    const int numChannels = std::min(buffer.getNumChannels(), getTotalNumInputChannels());
    if (numChannels == 0 || mAlgorithm == nullptr) {
        return;
    }

    // Mix down all input channels to a single mono signal.
    const float channelWeight = 1.0f / static_cast<float>(numChannels);
    for (int i = 0; i < buffer.getNumSamples(); ++i) {
        float mix = 0.0f;
        for (int channel = 0; channel < numChannels; ++channel) {
            mix += buffer.getSample(channel, i);
        }

        // Push the mixed and scaled sample into the ring buffer for FFT analysis.
        pushSample(mix * channelWeight);
    }
}

bool SpectrumProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
    const auto& in = layouts.getMainInputChannelSet();
    return in == layouts.getMainOutputChannelSet() && (in == juce::AudioChannelSet::mono() || in == juce::AudioChannelSet::stereo());
}

//==============================================================================
void SpectrumProcessor::pushSample(float sample) noexcept
{
    mRing.insert(Complex(sample, 0.0));

    if (++mSamplesSinceLastFrame < spectrum::kHopSize) {
        return;
    }

    mSamplesSinceLastFrame = 0;
    if (mRing.currentSize() == spectrum::kFftSize) {
        analyseLatestFrame();
    }
}

void SpectrumProcessor::analyseLatestFrame() noexcept
{
    // Copies the last kFftSize samples, oldest first.
    mRing.copy(mFftData.data(), spectrum::kFftSize);

    (*mAlgorithm)(mFftData);

    for (std::size_t bin = 0; bin < spectrum::kNumBins; ++bin) {
        mFrame.magnitude[bin] = static_cast<float>(std::abs(mFftData[bin]) * kBinScale);
    }

    // DC and Nyquist have no mirrored counterpart, so they must not be doubled.
    mFrame.magnitude.front() *= 0.5f;
    mFrame.magnitude.back() *= 0.5f;

    // A full queue means the GUI is behind: drop the frame instead of blocking the audio thread.
    mQueue->push(mFrame);
}

} // namespace inspect
