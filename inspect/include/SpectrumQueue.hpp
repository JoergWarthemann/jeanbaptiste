#ifndef INSPECT_SPECTRUMQUEUE_HPP_
#define INSPECT_SPECTRUMQUEUE_HPP_

#include <array>
#include <cstddef>

#include <boost/lockfree/spsc_queue.hpp>

namespace inspect::spectrum {

inline constexpr std::size_t kFftStage = 11;
inline constexpr std::size_t kFftSize = std::size_t{1} << kFftStage;
inline constexpr std::size_t kHopSize = kFftSize / 2; // 50 % overlap
inline constexpr std::size_t kNumBins = kFftSize / 2 + 1; // DC ... Nyquist
inline constexpr std::size_t kQueueCapacity = 16;

/**
 * One spectrum snapshot shared between the audio thread (producer) and the GUI thread (consumer).
 * Holds the spectrum after manipulation, i.e. what is heard.
 * Linear magnitudes where a full-scale sine reads 1.0 at its bin. Bin k is at k * sampleRate / kFftSize Hz.
 */
struct Frame {
    std::array<float, kNumBins> magnitude{};
};

/** Fixed capacity, so pushing and popping never allocate. */
using Queue = boost::lockfree::spsc_queue<Frame, boost::lockfree::capacity<kQueueCapacity>>;

} // namespace inspect::spectrum

#endif // INSPECT_SPECTRUMQUEUE_HPP_
