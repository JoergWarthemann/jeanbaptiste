#pragma once

#include "../basic/SineCosine.h"
#include <complex>
#include "../SubTask.h"

namespace jeanbaptiste::core
{
    /** Performs a split radix decimation in time FFT using template metaprogramming.
        \param SampleCnt ... The count of samples to be processed in this recursion level (stage)
        \param DirectionFactor ... Specifies the direction of the DFT (forward: 1, backward: -1)
        \param Complex ... The complex type.
    */
    template<typename SampleCnt,
             typename DirectionFactor,
             typename Complex>
    class SplitRadixDIT
        : public SubTask<SplitRadixDIT<SampleCnt, DirectionFactor, Complex>,
                         Complex>
    {
        SplitRadixDIT<std::integral_constant<unsigned, SampleCnt::value / 2>, DirectionFactor, Complex> recursionLevel_;

        void executeSimpleRadix2Butterflies(Complex* data) const
        {
            // Radix-2 butterflies without twiddle factor multiplication.
            for (unsigned int groupNodeIdx = 0, lShapedDualNodeDistance = 4; groupNodeIdx < SampleCnt::value; lShapedDualNodeDistance *= 4)
            {
                for (unsigned int segmentIdx = groupNodeIdx; segmentIdx < SampleCnt::value; segmentIdx += lShapedDualNodeDistance)
                {
                    Complex temp1 = data[segmentIdx] - data[segmentIdx + 1];
                    data[segmentIdx] += data[segmentIdx + 1];
                    data[segmentIdx + 1] = temp1;
                }

                groupNodeIdx = 2 * (lShapedDualNodeDistance - 1);
            }
        }

    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned groupNodeIdx = 0) const
        {
            executeSimpleRadix2Butterflies(data);
            executeRecursion(data, groupNodeIdx, SampleCnt::value);
        }

        void executeRecursion(Complex* data, unsigned groupNodeIdx, unsigned totalSampleCnt) const
        {
            // Recursion goes down. Calculation starts in the last recursion stage with n nodes and goes down: ..., 8, 4.
            recursionLevel_.executeRecursion(data, groupNodeIdx, totalSampleCnt);

            // dualNodeDistance is the distance between elements (successive nodes) of a 
            // dual tuple in a stage, e.g. 4, 8, 16, ... .
            // quaternaryNodeDistance is the distance between elements (successive nodes) of a 
            // quaternary tuple in a stage, e.g. 1, 2, 4, ... .
            // lShapedNodeDistance is the distance between elements (successive nodes) of a
            // 'L' shaped tuple in a stage, e.g. 8, 16, 16, 32, 32, 32, 32, ... .
            // segmentIdx is the index of segments being computed in a stage (considering skips based
            // on the 'L' shape of the tuples), e.g. 0, 12, 0, 1, 0, 1, 2, 3, ... .
            // tupleIdx is the index of 'L' shaped 4 point tuples being computed inside a segment,
            // e.g. 0, 8, 12, 0, 1, 0, 1, 2, 3, ... .
            unsigned dualNodeDistance = SampleCnt::value;
            unsigned quaternaryNodeDistance = SampleCnt::value >> 2;

            Complex twiddleMultiplier(
				static_cast<typename Complex::value_type>(-2.0 * basic::sine<typename Complex::value_type>(1, SampleCnt::value) * basic::sine<typename Complex::value_type>(1, SampleCnt::value)),
                static_cast<typename Complex::value_type>(DirectionFactor::value * basic::sine<typename Complex::value_type>(2, SampleCnt::value)));
            // Create transform factor.
            Complex twiddleFactor(1.0, 0.0);

            for (unsigned currentGroupIdx = groupNodeIdx, groupIdxEnd = (groupNodeIdx + quaternaryNodeDistance); currentGroupIdx < groupIdxEnd; ++currentGroupIdx)
            {
                // Create twiddle factors for Radix-4 butterflies x(4n + 1) and x(4n + 3).
				typename Complex::value_type temp = 1.5 - 0.5 * (twiddleFactor.real() * twiddleFactor.real() + twiddleFactor.imag() * twiddleFactor.imag());
                Complex wn4(twiddleFactor.real() * temp, twiddleFactor.imag() * temp);
                Complex w3n4 = wn4 * wn4 * wn4;

                auto segmentIdx = currentGroupIdx;
                auto lShapedNodeDistance = 2 * dualNodeDistance;

                while (segmentIdx < totalSampleCnt)
                {
                    for (auto tupleIdx = segmentIdx; tupleIdx < totalSampleCnt; tupleIdx += lShapedNodeDistance)
                    {
                        // 'L' shaped DIT split radix butterfly:
                        // X[r]            = Y[r]         +  (Z[r] * W^r + H[r] * W^3r)
                        // X[r + N/4]    = Y[r + N/4] - j(Z[r] * W^r - H[r] * W^3r)
                        // X[r + N/2]    = Y[r]         -  (Z[r] * W^r + H[r] * W^3r)
                        // X[r + 3N/4]    = Y[r + N/4] + j(Z[r] * W^r - H[r] * W^3r)
                        //
                        // with the following correspondences:
                        // data[idxNode0] -> X[r]
                        // data[idxNode1] -> X[r + N/4]
                        // data[idxNode2] -> X[r + N/2]
                        // data[idxNode3] -> X[r + 3N/4]

                        // Create tuple node indexes.
                        auto idxNode0    = tupleIdx;
                        auto idxNode1    = idxNode0        + quaternaryNodeDistance;
                        auto idxNode2    = idxNode1        + quaternaryNodeDistance;
                        auto idxNode3    = idxNode2        + quaternaryNodeDistance;

                        // Save temporary results.
                        auto temp1 = data[idxNode2] * wn4;
                        auto temp2 = data[idxNode3] * w3n4;
                        // Attention: in order to make bit reversal more comfortable slightly customize the formula.
                        // original: temp3 = temp1 - temp2;
                        auto temp3 = temp2 - temp1;
                    
                        temp1 += temp2;
                        temp2 = temp3 * Complex(0, DirectionFactor::value);

                        // Conduct 'L' shaped DIT butterfly.
                        data[idxNode3] = data[idxNode1] + temp2;
                        data[idxNode2] = data[idxNode0] - temp1;
                        data[idxNode1] -= temp2;
                        data[idxNode0] += temp1;
                    }

                    segmentIdx = 2 * lShapedNodeDistance - dualNodeDistance + currentGroupIdx;
                    lShapedNodeDistance <<= 2;
                }

                if ((currentGroupIdx + 1) < groupIdxEnd)
                    twiddleFactor += twiddleMultiplier * twiddleFactor;
            }
        }
    };

    /** Specialization for case SampleCnt=4, direction=1 (forward).
        \param Complex ... The complex type.
    */
    template<typename Complex>
    class SplitRadixDIT<std::integral_constant<unsigned, 4>, std::integral_constant<unsigned, 1>, Complex>
        : public SubTask<SplitRadixDIT<std::integral_constant<unsigned, 4>, std::integral_constant<int, 1>, Complex>,
                         Complex>
    {
    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned groupNodeIdx = 0) const
        {
            executeRecursion(data, groupNodeIdx, 4);
        }

        void executeRecursion(Complex* data, unsigned int groupNodeIdx, unsigned totalSampleCnt) const
        {
            auto segmentIdx = groupNodeIdx;
            auto lShapedNodeDistance = 8;

            while (segmentIdx < totalSampleCnt)
            {
                for (auto tupleIdx = segmentIdx; tupleIdx < totalSampleCnt; tupleIdx += lShapedNodeDistance)
                {
                    // Create tuple node indexes.
                    auto idxNode0    = tupleIdx;
                    auto idxNode1    = idxNode0        + 1;
                    auto idxNode2    = idxNode1        + 1;
                    auto idxNode3    = idxNode2        + 1;

                    // Temporary results.
                    typename Complex::value_type tmpReal0 = data[idxNode3].real() + data[idxNode2].real();
                    typename Complex::value_type tmpImag1 = data[idxNode3].imag() + data[idxNode2].imag();
                    typename Complex::value_type tmpReal2 = data[idxNode3].imag() - data[idxNode2].imag();
                    typename Complex::value_type tmpImag3 = data[idxNode3].real() - data[idxNode2].real();

                    auto temp0 = data[idxNode0];
                    auto temp1 = data[idxNode1];

                    // Conduct 'L' shaped DIT butterfly.
                    data[idxNode0].real(temp0.real() + tmpReal0);
                    data[idxNode0].imag(temp0.imag() + tmpImag1);

                    data[idxNode1].real(temp1.real() + tmpReal2);
                    data[idxNode1].imag(temp1.imag() - tmpImag3);

                    data[idxNode2].real(temp0.real() - tmpReal0);
                    data[idxNode2].imag(temp0.imag() - tmpImag1);

                    data[idxNode3].real(temp1.real() - tmpReal2);
                    data[idxNode3].imag(temp1.imag() + tmpImag3);

                }

                segmentIdx = 2 * lShapedNodeDistance - 4 + groupNodeIdx;
                lShapedNodeDistance <<= 2;
            }
        }
    };

    /** Specialization for case SampleCnt=4, direction=-1 (backward).
        \param Complex ... The complex type.
    */
    template<typename  Complex>
    class SplitRadixDIT<std::integral_constant<unsigned, 4>, std::integral_constant<int, -1>, Complex>
        : public SubTask<SplitRadixDIT<std::integral_constant<unsigned, 4>, std::integral_constant<int, -1>, Complex>,
                         Complex>
    {
    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned groupNodeIdx = 0) const
        {
            executeRecursion(data, groupNodeIdx, 4);
        }

        void executeRecursion(Complex* data, unsigned int groupNodeIdx, unsigned totalSampleCnt) const
        {
            auto segmentIdx = groupNodeIdx;
            auto lShapedNodeDistance = 8;

            while (segmentIdx < totalSampleCnt)
            {
                for (auto tupleIdx = segmentIdx; tupleIdx < totalSampleCnt; tupleIdx += lShapedNodeDistance)
                {
                    // Create tuple node indexes.
                    auto idxNode0    = tupleIdx;
                    auto idxNode1    = idxNode0        + 1;
                    auto idxNode2    = idxNode1        + 1;
                    auto idxNode3    = idxNode2        + 1;

                    // Temporary results.
                    typename Complex::value_type tmpReal0 = data[idxNode3].real() + data[idxNode2].real();
                    typename Complex::value_type tmpImag1 = data[idxNode3].imag() + data[idxNode2].imag();
                    typename Complex::value_type tmpReal2 = data[idxNode3].imag() - data[idxNode2].imag();
                    typename Complex::value_type tmpImag3 = data[idxNode3].real() - data[idxNode2].real();

                    auto temp0 = data[idxNode0];
                    auto temp1 = data[idxNode1];

                    // Conduct 'L' shaped DIT butterfly.
                    data[idxNode0].real(temp0.real() + tmpReal0);
                    data[idxNode0].imag(temp0.imag() + tmpImag1);

                    data[idxNode1].real(temp1.real() - tmpReal2);
                    data[idxNode1].imag(temp1.imag() + tmpImag3);

                    data[idxNode2].real(temp0.real() - tmpReal0);
                    data[idxNode2].imag(temp0.imag() - tmpImag1);

                    data[idxNode3].real(temp1.real() + tmpReal2);
                    data[idxNode3].imag(temp1.imag() - tmpImag3);
                }

                segmentIdx = 2 * lShapedNodeDistance - 4 + groupNodeIdx;
                lShapedNodeDistance <<= 2;
            }
        }
    };

    /** Specialization for case SampleCnt=2.
        \param DirectionFactor ... Specifies the direction of the DFT (forward: 1, backward: -1)
        \param Complex ... The complex type.
    */
    template<typename DirectionFactor,
             typename Complex>
    class SplitRadixDIT<std::integral_constant<unsigned, 2>, DirectionFactor, Complex>
        : public SubTask<SplitRadixDIT<std::integral_constant<unsigned, 2>, DirectionFactor, Complex>,
                         Complex>
    {
    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned groupNodeIdx = 0)const
        {
            executeRecursion(data, groupNodeIdx, 2);
        }

        void executeRecursion(Complex* data, unsigned int groupNodeIdx, unsigned totalSampleCnt) const
        {}
    };

    /** Performs a split radix decimation in frequency FFT using template metaprogramming.
        \param SampleCnt ... The count of samples to be processed in this recursion level (stage)
        \param DirectionFactor ... Specifies the direction of the DFT (forward: 1, backward: -1)
        \param Complex ... The complex type.
    */
    template<typename SampleCnt,
             typename DirectionFactor,
             typename Complex>
    class SplitRadixDIF
        : public SubTask<SplitRadixDIF<SampleCnt, DirectionFactor, Complex>,
                         Complex>
    {
        SplitRadixDIF<std::integral_constant<unsigned, SampleCnt::value / 2>, DirectionFactor, Complex> recursionLevel_;

        void executeSimpleRadix2Butterflies(Complex* data) const
        {
            // Radix-2 butterflies without twiddle factor multiplication.
            for (auto groupNodeIdx = 0, lShapedDualNodeDistance = 4; groupNodeIdx < SampleCnt::value; lShapedDualNodeDistance *= 4)
            {
                for (auto segmentIdx = groupNodeIdx; segmentIdx < SampleCnt::value; segmentIdx += lShapedDualNodeDistance)
                {
                    auto temp1 = data[segmentIdx] - data[segmentIdx + 1];
                    data[segmentIdx] += data[segmentIdx + 1];
                    data[segmentIdx + 1] = temp1;
                }

                groupNodeIdx = 2 * (lShapedDualNodeDistance - 1);
            }
        }

    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned groupNodeIdx = 0) const
        {
            executeRecursion(data, groupNodeIdx, SampleCnt::value);
            executeSimpleRadix2Butterflies(data);
        }

        void executeRecursion(Complex* data, unsigned groupNodeIdx, unsigned totalSampleCnt) const
        {
            // dualNodeDistance is the distance between elements (successive nodes) of a 
            // dual tuple in a stage, e.g. ..., 16, 8, 4.
            // quaternaryNodeDistance is the distance between elements (successive nodes) of a
            // quaternary tuple in a stage, e.g. ..., 4, 2, 1.
            // lShapedNodeDistance is the distance between elements (successive nodes) of a
            // 'L' shaped tuple in a stage, e.g. ..., 32, 32, 32, 32, 16, 16, 8.
            // segmentIdx is the index of segments being computed in a stage (considering skips based
            // on the 'L' shape of the tuples), e.g. ..., 0, 1, 2, 3, 0, 1, 0, 12.
            // tupleIdx is the index of 'L' shaped 4 point tuples being computed inside a segment,
            // e.g. ..., 0, 1, 2, 3, 0, 1, 0, 8, 12.
            auto dualNodeDistance = SampleCnt::value;
            auto quaternaryNodeDistance = SampleCnt::value >> 2;

            Complex twiddleMultiplier(
				static_cast<typename Complex::value_type>(-2.0 * basic::sine<typename Complex::value_type>(1, SampleCnt::value) * basic::sine<typename Complex::value_type>(1, SampleCnt::value)),
                static_cast<typename Complex::value_type>(DirectionFactor::value * basic::sine<typename Complex::value_type>(2, SampleCnt::value)));
            // Create transform factor.
            Complex twiddleFactor(1.0, 0.0);

            for (unsigned currentGroupIdx = groupNodeIdx, groupIdxEnd = (groupNodeIdx + quaternaryNodeDistance); currentGroupIdx < groupIdxEnd; ++currentGroupIdx)
            {
                // Create twiddle factors for Radix-4 butterflies x(4n + 1) and x(4n + 3).
                typename Complex::value_type temp = 1.5 - 0.5 * (twiddleFactor.real() * twiddleFactor.real() + twiddleFactor.imag() * twiddleFactor.imag());
                Complex wn4(twiddleFactor.real() * temp, twiddleFactor.imag() * temp);
                Complex w3n4 = wn4 * wn4 * wn4;

                auto segmentIdx = currentGroupIdx;
                auto lShapedNodeDistance = 2 * dualNodeDistance;

                while (segmentIdx < totalSampleCnt)
                {
                    for (auto tupleIdx = segmentIdx; tupleIdx < totalSampleCnt; tupleIdx += lShapedNodeDistance)
                    {
                        // 'L' shaped DIF split radix butterfly:
                        // data[idxNode0] -> y[l]        = x[l]                   +   x[l + N/2]
                        // data[idxNode1] -> y[l + N/4] = x[l + N/4]           +   x[l + 3N/4]
                        // data[idxNode2] -> z[l]        = ((x[l] - x[l + N/2]) - j(x[l + N/4] - x[l + 3N/4])) * W^l
                        // data[idxNode3] -> h[l]        = ((x[l] - x[l + N/2]) + j(x[l + N/4] - x[l + 3N/4])) * W^3l
                        //
                        // with the following correspondences:
                        // data[idxNode0] -> y[l]
                        // data[idxNode1] -> y[l + N/4]
                        // data[idxNode2] -> z[l]
                        // data[idxNode3] -> h[l]

                        // Create tuple node indexes.
                        auto idxNode0    = tupleIdx;
                        auto idxNode1    = idxNode0        + quaternaryNodeDistance;
                        auto idxNode2    = idxNode1        + quaternaryNodeDistance;
                        auto idxNode3    = idxNode2        + quaternaryNodeDistance;

                        // Save temporary results for Radix-4.
                        auto temp1 = data[idxNode0] - data[idxNode2];
                        auto temp2 = Complex(0, DirectionFactor::value) * (data[idxNode1] - data[idxNode3]);

                        // Conduct 'L' shaped DIF butterfly.
                        data[idxNode0] += data[idxNode2];
                        data[idxNode1] += data[idxNode3];
                        data[idxNode2] = wn4  * (temp1 + temp2);
                        data[idxNode3] = w3n4 * (temp1 - temp2);

                        // Attention: in order to make bit reversal more comfortable slightly customize the formula.
                        // original: data[idxNode2] = wn4  * (temp1 - temp2);
                        //             data[idxNode3] = w3n4 * (temp1 + temp2);
                    }

                    segmentIdx = 2 * lShapedNodeDistance - dualNodeDistance + currentGroupIdx;
                    lShapedNodeDistance <<= 2;
                }

                if ((currentGroupIdx + 1) < groupIdxEnd)
                    twiddleFactor += twiddleMultiplier * twiddleFactor;
            }

            // Recursion goes down. Calculation starts in the last recursion stage with n nodes and goes down: ..., 8, 4.
            recursionLevel_.executeRecursion(data, groupNodeIdx, totalSampleCnt);
        }
    };

    /** Specialization for case SampleCnt=4, direction=1 (forward).
        \param Complex ... The complex type.
    */
    template<typename Complex>
    class SplitRadixDIF<std::integral_constant<unsigned, 4>, std::integral_constant<int, 1>, Complex>
        : public SubTask<SplitRadixDIF<std::integral_constant<unsigned, 4>, std::integral_constant<int, 1>, Complex>,
                         Complex>
    {
    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned groupNodeIdx = 0) const
        {
            executeRecursion(data, groupNodeIdx, 4);
        }

        void executeRecursion(Complex* data, unsigned int groupNodeIdx, unsigned totalSampleCnt) const
        {
            auto segmentIdx = groupNodeIdx;
            auto lShapedNodeDistance = 8;

            while (segmentIdx < totalSampleCnt)
            {
                for (auto tupleIdx = segmentIdx; tupleIdx < totalSampleCnt; tupleIdx += lShapedNodeDistance)
                {
                    // Create tuple node indexes.
                    auto idxNode0    = tupleIdx;
                    auto idxNode1    = idxNode0        + 1;
                    auto idxNode2    = idxNode1        + 1;
                    auto idxNode3    = idxNode2        + 1;

                    // Temporary results.
                    typename Complex::value_type tmpReal3 = data[idxNode0].real() - data[idxNode2].real();
                    typename Complex::value_type tmpReal4 = data[idxNode1].imag() - data[idxNode3].imag();
                    typename Complex::value_type tmpImag3 = data[idxNode0].imag() - data[idxNode2].imag();
                    typename Complex::value_type tmpImag4 = data[idxNode1].real() - data[idxNode3].real();

                    // Conduct 'L' shaped DIF butterfly.
                    data[idxNode0] += data[idxNode2];
                    data[idxNode1] += data[idxNode3];
                
                    data[idxNode2].real(tmpReal3 - tmpReal4);
                    data[idxNode2].imag(tmpImag3 + tmpImag4);

                    data[idxNode3].real(tmpReal3 + tmpReal4);
                    data[idxNode3].imag(tmpImag3 - tmpImag4);

                    // Attention: in order to make bit reversal more comfortable slightly customize the formula.
                    // original: data[idxNode2] = wn4  * (temp1 - temp2);
                    //             data[idxNode3] = w3n4 * (temp1 + temp2);
                }

                segmentIdx = 2 * lShapedNodeDistance - 4 + groupNodeIdx;
                lShapedNodeDistance <<= 2;
            }
        }
    };

    /** Specialization for case SampleCnt=4, direction=-1 (backward).
        \param Complex ... The complex type.
    */
    template<typename Complex>
    class SplitRadixDIF<std::integral_constant<unsigned, 4>, std::integral_constant<int, -1>, Complex>
        : public SubTask<SplitRadixDIF<std::integral_constant<unsigned, 4>, std::integral_constant<int, -1>, Complex>,
                         Complex>
    {
    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned groupNodeIdx = 0) const
        {
            executeRecursion(data, groupNodeIdx, 4);
        }

        void executeRecursion(Complex* data, unsigned int groupNodeIdx, unsigned totalSampleCnt) const
        {
            auto segmentIdx = groupNodeIdx;
            auto lShapedNodeDistance = 8;

            while (segmentIdx < totalSampleCnt)
            {
                for (auto tupleIdx = segmentIdx; tupleIdx < totalSampleCnt; tupleIdx += lShapedNodeDistance)
                {
                    // Create tuple node indexes.
                    auto idxNode0    = tupleIdx;
                    auto idxNode1    = idxNode0        + 1;
                    auto idxNode2    = idxNode1        + 1;
                    auto idxNode3    = idxNode2        + 1;

                    // Temporary results.
                    typename Complex::value_type tmpReal3 = data[idxNode0].real() - data[idxNode2].real();
                    typename Complex::value_type tmpReal4 = data[idxNode1].imag() - data[idxNode3].imag();
                    typename Complex::value_type tmpImag3 = data[idxNode0].imag() - data[idxNode2].imag();
                    typename Complex::value_type tmpImag4 = data[idxNode1].real() - data[idxNode3].real();

                    // Conduct 'L' shaped DIF butterfly.
                    data[idxNode0] += data[idxNode2];
                    data[idxNode1] += data[idxNode3];
                
                    data[idxNode2].real(tmpReal3 + tmpReal4);
                    data[idxNode2].imag(tmpImag3 - tmpImag4);

                    data[idxNode3].real(tmpReal3 - tmpReal4);
                    data[idxNode3].imag(tmpImag3 + tmpImag4);

                    // Attention: in order to make bit reversal more comfortable slightly customize the formula.
                    // original: data[idxNode2] = wn4  * (temp1 - temp2);
                    //             data[idxNode3] = w3n4 * (temp1 + temp2);
                }

                segmentIdx = 2 * lShapedNodeDistance - 4 + groupNodeIdx;
                lShapedNodeDistance <<= 2;
            }
        }
    };

    /** Specialization for case SampleCnt=2.
        \param DirectionFactor ... Specifies the direction of the DFT (forward: 1, backward: -1)
        \param Complex ... The complex type.
    */
    template<typename DirectionFactor,
             typename Complex>
    class SplitRadixDIF<std::integral_constant<unsigned, 2>, DirectionFactor, Complex>
        : public SubTask<SplitRadixDIF<std::integral_constant<unsigned, 2>, DirectionFactor, Complex>,
                         Complex>
    {
    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned groupNodeIdx = 0) const
        {
            executeRecursion(data, groupNodeIdx, 2);
        }

        void executeRecursion(Complex* data, unsigned int groupNodeIdx, unsigned totalSampleCnt) const
        { }
    };
}