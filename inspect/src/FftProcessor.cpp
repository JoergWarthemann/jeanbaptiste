#include "FftProcessor.hpp"

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

FftProcessor::FftProcessor()
    : AudioProcessor(BusesProperties()
                         .withInput("Input", juce::AudioChannelSet::stereo(), true)
                         .withOutput("Output", juce::AudioChannelSet::stereo(), true))
{
}

FftProcessor::~FftProcessor() = default;

//==============================================================================
void FftProcessor::prepareToPlay(double, int)
{
    mAlgorithm = FftFactory{}.getAlgorithm(spectrum::kFftStage);
    mFftData.assign(spectrum::kFftSize, Complex{});

    mRing.fill(0.0f);
    mWritePosition = 0;
    mSamplesSinceLastFrame = 0;
    mNumValidSamples = 0;
}

void FftProcessor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer&)
{
    juce::ScopedNoDenormals noDenormals;

    const int numChannels = std::min(buffer.getNumChannels(), getTotalNumInputChannels());
    if (numChannels == 0 || mAlgorithm == nullptr) {
        return;
    }

    const float channelWeight = 1.0f / static_cast<float>(numChannels);
    for (int i = 0; i < buffer.getNumSamples(); ++i) {
        float mix = 0.0f;
        for (int channel = 0; channel < numChannels; ++channel) {
            mix += buffer.getSample(channel, i);
        }

        pushSample(mix * channelWeight);
    }
}

bool FftProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
    const auto& in = layouts.getMainInputChannelSet();
    return in == layouts.getMainOutputChannelSet() && (in == juce::AudioChannelSet::mono() || in == juce::AudioChannelSet::stereo());
}

//==============================================================================
void FftProcessor::pushSample(float sample) noexcept
{
    mRing[mWritePosition] = sample;
    mWritePosition = (mWritePosition + 1) & (spectrum::kFftSize - 1);

    mNumValidSamples = std::min(mNumValidSamples + 1, spectrum::kFftSize);
    if (++mSamplesSinceLastFrame < spectrum::kHopSize) {
        return;
    }

    mSamplesSinceLastFrame = 0;
    if (mNumValidSamples == spectrum::kFftSize) {
        analyseLatestFrame();
    }
}

void FftProcessor::analyseLatestFrame() noexcept
{
    // mWritePosition is the oldest sample of the ring.
    for (std::size_t i = 0; i < spectrum::kFftSize; ++i) {
        mFftData[i] = Complex(mRing[(mWritePosition + i) & (spectrum::kFftSize - 1)], 0.0);
    }

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
