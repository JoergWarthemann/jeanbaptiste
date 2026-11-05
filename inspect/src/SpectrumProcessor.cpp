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
    for (auto& channel : mChannels) {
        channel.forward = ForwardFactory{}.getAlgorithm(spectrum::kFftStage);
        channel.inverse = InverseFactory{}.getAlgorithm(spectrum::kFftStage);
        channel.fftData.assign(spectrum::kFftSize, Complex{});
        channel.overlapAdd.reset();
        channel.pending.fill(0.0);
        channel.ring.reset(spectrum::kFftSize);
    }

    // A sine of amplitude A at bin k gives |X[k]| = A / 2 * sum(w), and the sum of the sine window
    // w(n) = sin(pi * n / N) is cot(pi / (2 * N)).
    mMagnitudeScale = 2.0 * std::tan(std::numbers::pi / (2.0 * spectrum::kFftSize));

    mPendingPosition = 0;
    mSamplesSinceLastFrame = 0;

    setLatencySamples(static_cast<int>(spectrum::kFftSize));
}

void SpectrumProcessor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer&)
{
    juce::ScopedNoDenormals noDenormals;

    const int numChannels = std::min(buffer.getNumChannels(), getTotalNumInputChannels());
    if (numChannels == 0 || mChannels.front().forward == nullptr) {
        return;
    }

    for (int sample = 0; sample < buffer.getNumSamples(); ++sample) {
        // Process each sample for all channels, exchanging it with the output from the pending buffer.
        for (int channel = 0; channel < numChannels; ++channel) {
            buffer.setSample(
                channel,
                sample,
                exchangeSample(
                    buffer.getSample(channel, sample),
                    mChannels[static_cast<std::size_t>(channel)]));
        }

        // Advance the pending position and check if a new frame should be processed.
        mPendingPosition = (mPendingPosition + 1) % spectrum::kHopSize;
        // The pending position wraps around the hop size, ensuring overlap-add works correctly.
        if (++mSamplesSinceLastFrame == spectrum::kHopSize) {
            mSamplesSinceLastFrame = 0;
            if (mChannels.front().ring.currentSize() == spectrum::kFftSize) {
                processSampleFrameAndSpectrum(numChannels);
            }
        }
    }

    for (int channel = numChannels; channel < buffer.getNumChannels(); ++channel) {
        buffer.clear(channel, 0, buffer.getNumSamples());
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
float SpectrumProcessor::exchangeSample(float input, ChannelState& channel) noexcept
{
    // The finished samples of the last frame are output one by one. The next frame completes right after the last of
    // them has been output, which makes the latency exactly kFftSize.
    const float output = static_cast<float>(channel.pending[mPendingPosition]);
    channel.ring.insert(Complex(input, 0.0));
    return output;
}

void SpectrumProcessor::processSampleFrameAndSpectrum(int numChannels) noexcept
{
    for (int channelIndex = 0; channelIndex < numChannels; ++channelIndex) {
        auto& channel = mChannels[static_cast<std::size_t>(channelIndex)];
        channel.ring.copy(channel.fftData.data(), spectrum::kFftSize);
        (*channel.forward)(channel.fftData);
        manipulateSpectrum(std::span<Complex, spectrum::kNumBins>(channel.fftData.data(), spectrum::kNumBins));
    }

    publishSpectrum(numChannels);

    for (int channelIndex = 0; channelIndex < numChannels; ++channelIndex) {
        auto& channel = mChannels[static_cast<std::size_t>(channelIndex)];
        (*channel.inverse)(channel.fftData);
        channel.overlapAdd(
            std::span<const Complex, spectrum::kFftSize>(channel.fftData.data(), spectrum::kFftSize),
            channel.pending);
    }

    mPendingPosition = 0;
}

void SpectrumProcessor::publishSpectrum(int numChannels) noexcept
{
    for (std::size_t bin = 0; bin < spectrum::kNumBins; ++bin) {
        double magnitude = 0.0;
        for (int channel = 0; channel < numChannels; ++channel) {
            magnitude += std::abs(mChannels[static_cast<std::size_t>(channel)].fftData[bin]);
        }
        mFrame.magnitude[bin] = static_cast<float>(magnitude * mMagnitudeScale / numChannels);
    }

    // DC and Nyquist have no mirrored counterpart, so they must not be doubled.
    mFrame.magnitude.front() *= 0.5f;
    mFrame.magnitude.back() *= 0.5f;

    // A full queue means the GUI is behind: drop the frame instead of blocking the audio thread.
    mQueue->push(mFrame);
}

} // namespace inspect
