#include "SpectrumProcessor.hpp"

#include <cmath>
#include <numbers>

#include "AlgorithmFactory.hpp"
#include "Options.hpp"

namespace inspect {

namespace {

using Complex = std::complex<double>;

// Real FFT, radix-2 DIT, forward, sine window applied to the input (analysis), no normalization.
using ForwardFactory = jb::AlgorithmFactory<
    spectrum::kFftStage, spectrum::kFftStage + 1,
    jb::options::Radix_2,
    jb::options::Decimation_In_Time,
    jb::options::Direction_Forward,
    jb::options::Window_Cosine,
    jb::options::Normalization_No,
    Complex,
    jb::options::Data_Real,
    jb::options::Window_Before>;

// The matching inverse with the same sine window applied to the result (synthesis).
// Division by the length makes forward + inverse an identity, and the sine windows reconstruct
// the signal when overlap-added with 50 % overlap.
using InverseFactory = jb::AlgorithmFactory<
    spectrum::kFftStage, spectrum::kFftStage + 1,
    jb::options::Radix_2,
    jb::options::Decimation_In_Time,
    jb::options::Direction_Backward,
    jb::options::Window_Cosine,
    jb::options::Normalization_Division_By_Length,
    Complex,
    jb::options::Data_Real,
    jb::options::Window_After>;

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
    mForward = ForwardFactory{}.getAlgorithm(spectrum::kFftStage);
    mInverse = InverseFactory{}.getAlgorithm(spectrum::kFftStage);
    mFftData.assign(spectrum::kFftSize, Complex{});

    // A sine of amplitude A at bin k gives |X[k]| = A / 2 * sum(w), and the sum of the sine window
    // w(n) = sin(pi * n / N) is cot(pi / (2 * N)).
    mMagnitudeScale = 2.0 * std::tan(std::numbers::pi / (2.0 * spectrum::kFftSize));

    mOverlapAdd.assign(spectrum::kFftSize, 0.0);
    mOverlapPosition = 0;

    mRing.reset(spectrum::kFftSize);
    mSamplesSinceLastFrame = 0;

    setLatencySamples(static_cast<int>(spectrum::kFftSize));
}

void SpectrumProcessor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer&)
{
    juce::ScopedNoDenormals noDenormals;

    const int numChannels = std::min(buffer.getNumChannels(), getTotalNumInputChannels());
    if (numChannels == 0 || mForward == nullptr) {
        return;
    }

    // Mix down all input channels to a single mono signal, process it and write the result to all channels.
    const float channelWeight = 1.0f / static_cast<float>(numChannels);
    for (int i = 0; i < buffer.getNumSamples(); ++i) {
        float mix = 0.0f;
        for (int channel = 0; channel < numChannels; ++channel) {
            mix += buffer.getSample(channel, i);
        }

        const float processed = processSample(mix * channelWeight);
        for (int channel = 0; channel < buffer.getNumChannels(); ++channel) {
            buffer.setSample(channel, i, processed);
        }
    }
}

bool SpectrumProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
    const auto& in = layouts.getMainInputChannelSet();
    return in == layouts.getMainOutputChannelSet() && (in == juce::AudioChannelSet::mono() || in == juce::AudioChannelSet::stereo());
}

void SpectrumProcessor::manipulateSpectrum(std::span<Complex, spectrum::kNumBins>) noexcept
{
}

//==============================================================================
float SpectrumProcessor::processSample(float input) noexcept
{
    // The output sample is finished and read before the input is added, which makes the latency exactly kFftSize.
    const float output = static_cast<float>(mOverlapAdd[mOverlapPosition]);
    mOverlapAdd[mOverlapPosition] = 0.0;
    mOverlapPosition = (mOverlapPosition + 1) % spectrum::kFftSize;

    mRing.insert(Complex(input, 0.0));

    if (++mSamplesSinceLastFrame < spectrum::kHopSize) {
        return output;
    }

    mSamplesSinceLastFrame = 0;
    if (mRing.currentSize() == spectrum::kFftSize) {
        processLatestFrame();
    }

    return output;
}

void SpectrumProcessor::processLatestFrame() noexcept
{
    // Copies the last kFftSize samples, oldest first.
    mRing.copy(mFftData.data(), spectrum::kFftSize);

    (*mForward)(mFftData);
    manipulateSpectrum(std::span<Complex, spectrum::kNumBins>(mFftData.data(), spectrum::kNumBins));
    publishSpectrum();
    (*mInverse)(mFftData);

    // The frame starts at the oldest sample which has not been output yet, which is mOverlapPosition.
    for (std::size_t i = 0; i < spectrum::kFftSize; ++i) {
        mOverlapAdd[(mOverlapPosition + i) % spectrum::kFftSize] += mFftData[i].real();
    }
}

void SpectrumProcessor::publishSpectrum() noexcept
{
    for (std::size_t bin = 0; bin < spectrum::kNumBins; ++bin) {
        mFrame.magnitude[bin] = static_cast<float>(std::abs(mFftData[bin]) * mMagnitudeScale);
    }

    // DC and Nyquist have no mirrored counterpart, so they must not be doubled.
    mFrame.magnitude.front() *= 0.5f;
    mFrame.magnitude.back() *= 0.5f;

    // A full queue means the GUI is behind: drop the frame instead of blocking the audio thread.
    mQueue->push(mFrame);
}

} // namespace inspect
