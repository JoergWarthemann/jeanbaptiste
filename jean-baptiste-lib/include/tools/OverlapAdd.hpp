#ifndef JB_TOOLS_OVERLAPADD_HPP_
#define JB_TOOLS_OVERLAPADD_HPP_

#include <array>
#include <cstddef>
#include <span>
#include <type_traits>

namespace jb::tools {

/**
 * Example: FrameSize = 8, HopSize = 2, so kOverlap = 6
 *          The tail holds t0..t5 and the new frame holds f0..f7 (real parts).
 *
 * Step 1: output the first HopSize samples
 *     mTail (old):   [ t0 | t1 | t2 | t3 | t4 | t5 ]
 *                       +    +
 *     frame:         [ f0 | f1 | f2 | f3 | f4 | f5 | f6 | f7 ]
 *                       =    =
 *     out:           [ t0+f0 | t1+f1 ]        -> finished samples
 *
 * Step 2: shift the tail left by HopSize and add the rest of the frame
 *     old tail:      [ t0 | t1 | t2 | t3 | t4 | t5 ]
 *                      ╰────╯    ▲    ▲    ▲    ▲
 *                     consumed   │    │    │    │   shift left by 2
 *                                │    │    │    │   (idx + Hop < kOverlap)
 *     shifted old:             [ t2 | t3 | t4 | t5 | 0  | 0  ]
 *                                 +    +    +    +    +    +
 *     frame[Hop..]:            [ f2 | f3 | f4 | f5 | f6 | f7 ]
 *                                 =    =    =    =    =    =
 *     new tail:                [ t2+f2 | t3+f3 | t4+f4 | t5+f5 |  f6   |  f7   ]
 *                                idx 0   idx 1   idx 2   idx 3   idx 4   idx 5
 */

/**
 * Overlap-add of consecutive frames which advance by HopSize samples, e.g. the synthesis stage of an STFT.
 *
 * Each call adds the real part of one frame to the unfinished tail of the previous frames and returns the HopSize
 * samples which are complete afterwards. Only the tail of FrameSize - HopSize samples is kept as state, so nothing
 * allocates. Frames overlap by FrameSize - HopSize samples, and more than two frames may overlap if HopSize is
 * smaller than half of FrameSize.
 *
 * The type only sums. If the frames' windows don't add up to one for the chosen hop size, the caller has to scale.
 * The first output of a stream contains only the first frame, so it is not complete (it fades in).
 *
 * \tparam FrameSize ... The number of samples of a frame, as std::integral_constant.
 * \tparam HopSize ... The number of samples which the frames advance by, as std::integral_constant.
 * \tparam Complex ... The complex data type of the frames, of which the real part is used.
 */
template <typename FrameSize, typename HopSize, typename Complex>
    requires std::is_integral_v<typename FrameSize::value_type> &&
    std::is_integral_v<typename HopSize::value_type> &&
    std::is_floating_point_v<typename Complex::value_type> &&
    (HopSize::value > 0) && (HopSize::value <= FrameSize::value)
class OverlapAdd {
public:
    using Value = typename Complex::value_type;

    /**
     * Adds a frame and returns the samples which are complete now.
     * \param[in] frame ... The next frame of FrameSize elements.
     * \param[out] out ... Receives the next HopSize finished samples.
     */
    void operator()(std::span<const Complex, FrameSize::value> frame, std::span<Value, HopSize::value> out) noexcept
    {
        // Output the first HopSize frame samples, adding the previous tail where it overlaps.
        for (std::size_t idx = 0; idx < HopSize::value; ++idx) {
            out[idx] = (idx < kOverlap ? mTail[idx] : Value{}) + frame[idx].real();
        }

        // Shift the tail by one hop and add the rest of the frame. Ascending order is safe for the in place shift.
        for (std::size_t idx = 0; idx < kOverlap; ++idx) {
            mTail[idx] =
                (idx + HopSize::value < kOverlap ? mTail[idx + HopSize::value] : Value{}) +
                frame[idx + HopSize::value].real();
        }
    }

    /**
     * Discards the unfinished tail, so the next frame starts a new stream.
     */
    void reset() noexcept
    {
        mTail.fill(Value{});
    }

private:
    static constexpr std::size_t kOverlap = FrameSize::value - HopSize::value;

    std::array<Value, kOverlap> mTail{};
};

} // namespace jb::tools

#endif // JB_TOOLS_OVERLAPADD_HPP_
