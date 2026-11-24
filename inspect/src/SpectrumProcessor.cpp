#include "SpectrumProcessor.hpp"

#include <algorithm>
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
// the signal when overlap-added with 75 % overlap.
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
        channel.pitchShifter.prepare(spectrum::kFftSize, spectrum::kHopSize);
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

    const auto numSamples = static_cast<std::size_t>(buffer.getNumSamples());
    // Split the buffer into chunks of hop size and process each chunk.
    // Incomplete blocks are getting completed in the next block.
    //   [ b0 | b1 | b2 | b3 | b4 | b5 | b6 | b7 | b8 | b9 ]
    //     ╰─── chunk 1 ────╯ ╰─── chunk 2 ────╯ ╰─ chunk 3 ─╯
    for (std::size_t offset = 0; offset < numSamples;) {
        // Determine how many samples to process in this chunk: min of remaining space in current hop and
        // remaining samples in the buffer.
        const auto count = std::min(numSamples - offset, spectrum::kHopSize - mSamplesSinceLastFrame);

        for (int channelIndex = 0; channelIndex < numChannels; ++channelIndex) {
            auto& channel = mChannels[static_cast<std::size_t>(channelIndex)];
            auto samples = std::span<float>{buffer.getWritePointer(channelIndex), numSamples}.subspan(offset, count);
            auto output = std::span<const double>{channel.pending}.subspan(mPendingPosition, count);
            // Exchange the current chunk of samples with the pending output for this channel.
            exchangeSamples(samples, output, channel);
        }

        offset += count;
        mPendingPosition = (mPendingPosition + count) % spectrum::kHopSize;
        mSamplesSinceLastFrame += count;
        // If we've accumulated enough samples for a hop, process the frame and spectrum.
        if (mSamplesSinceLastFrame == spectrum::kHopSize) {
            mSamplesSinceLastFrame = 0;
            // Check if the ring buffer has enough samples for an FFT frame before processing.
            if (mChannels.front().ring.currentSize() == spectrum::kFftSize) {
                processSampleFrameAndSpectrum(numChannels);
            }
        }
    }

    // Silence buffer channels that weren't processed, preventing stale or undefined samples from reaching the output.
    for (int channel = numChannels; channel < buffer.getNumChannels(); ++channel) {
        buffer.clear(channel, 0, buffer.getNumSamples());
    }
}

bool SpectrumProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
    const auto& in = layouts.getMainInputChannelSet();
    return in == layouts.getMainOutputChannelSet() && (in == juce::AudioChannelSet::mono() || in == juce::AudioChannelSet::stereo());
}

void SpectrumProcessor::manipulateSpectrum(std::size_t channel, std::span<Complex, spectrum::kNumBins> bins) noexcept
{
    static constexpr float kPitchShiftThreshold{1e-6f};
    auto& shifter = mChannels[channel].pitchShifter;
    const float semitones = mPitchSemitones.load(std::memory_order_relaxed);
    if (std::abs(semitones) < kPitchShiftThreshold) {
        shifter.passThrough(bins);
        return;
    }

    shifter.process(bins, std::exp2(static_cast<double>(semitones) / 12.0));
}

void SpectrumProcessor::exchangeSamples(std::span<float> samples, std::span<const double> output, ChannelState& channel) noexcept
{
    for (float input : samples) {
        channel.ring.insert(Complex{input, 0.0});
    }

    std::ranges::transform(output, samples.begin(), [](double sample) {
        return static_cast<float>(sample);
    });
}

void SpectrumProcessor::processSampleFrameAndSpectrum(int numChannels) noexcept
{
    for (int channelIndex = 0; channelIndex < numChannels; ++channelIndex) {
        auto& channel = mChannels[static_cast<std::size_t>(channelIndex)];
        channel.ring.copy(channel.fftData.data(), spectrum::kFftSize);
        (*channel.forward)(channel.fftData);
        manipulateSpectrum(
            static_cast<std::size_t>(channelIndex),
            std::span<Complex, spectrum::kNumBins>(channel.fftData.data(), spectrum::kNumBins));
    }

    publishSpectrum(numChannels);

    for (int channelIndex = 0; channelIndex < numChannels; ++channelIndex) {
        auto& channel = mChannels[static_cast<std::size_t>(channelIndex)];
        (*channel.inverse)(channel.fftData);
        // Create pending samples as the input samples for the next FFT by using the inverse FFT output to create
        // overlapping data.
        // input hops:     [ hop 0 ][ hop 1 ][ hop 2 ][ hop 3 ][ hop 4 ][ hop 5 ]
        // ring N=16:      [  ..... frame k, covers hops 0..3 ......  ]
        //                          [  ..... frame k+1, hops 1..4 ..... ]
        //                                   (75 % overlap, one frame per hop)
        // frame k ready:                                   ▲ at the end of hop 3
        // pending (k):                                     [ q(k) played during hop 4 ]
        // output hops:    [ silence ][ silence ][ silence ][ silence ][ q(k) ][ q(k+1) ]
        //
        // The input of hop 0 comes out of hop 4, so the delay is N = kFftSize samples (setLatencySamples).
        channel.overlapAdd(
            std::span<const Complex, spectrum::kFftSize>(channel.fftData.data(), spectrum::kFftSize),
            channel.pending);
        // Analysis and synthesis both apply the window, so the overlapping frames add up to sum(w^2) = N / (2 * H) = 2
        // instead of 1. kOverlapGain removes that factor.
        std::ranges::transform(channel.pending, channel.pending.begin(), [](double sample) {
            return sample * spectrum::kOverlapGain;
        });
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
