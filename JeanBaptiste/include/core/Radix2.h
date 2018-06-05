#pragma once

#include "../basic/SineCosine.h"
#include <complex>
#include "../SubTask.h"

namespace jeanbaptiste::core
{
    /** Performs a radix 2 decimation in time FFT using template metaprogramming.
        \param SampleCnt ... The count of samples to be processed in this recursion level (stage)
        \param DirectionFactor ... Specifies the direction of the DFT (forward: 1, backward: -1)
        \param Complex ... The complex type.
    */
    template<typename SampleCnt,
             typename DirectionFactor,
             typename Complex>
    class Radix2DIT
        : public SubTask<Radix2DIT<SampleCnt, DirectionFactor, Complex>,
                         Complex>
    {
        Radix2DIT<std::integral_constant<unsigned, SampleCnt::value / 2>, DirectionFactor, Complex> recursionLevel_;

    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned groupNodeIdx = 0) const
        {
            // dualNodeDistance is the distance between elements (successive nodes) of a 
            // dual tuple, e.g. ..., 8, 4, 2, 1.
            unsigned dualNodeDistance = SampleCnt::value >> 1;

            // Recursion goes down. Calculation starts in the last recursion stage with 2 nodes and goes up: 4, 8, ...
            recursionLevel_.apply(data, groupNodeIdx);
            recursionLevel_.apply(data, groupNodeIdx + dualNodeDistance);

            // Create twiddle factor multiplier for trigonometric recurrence.
            Complex twiddleMultiplier(
                static_cast<typename Complex::value_type>(-2.0 * basic::sine<typename Complex::value_type>(1, SampleCnt::value) * basic::sine<typename Complex::value_type>(1, SampleCnt::value)),
                static_cast<typename Complex::value_type>(DirectionFactor::value * basic::sine<typename Complex::value_type>(2, SampleCnt::value)));
            // Create transform factor.
            Complex twiddleFactor(1.0, 0.0);

            // Run through dual nodes within the current group.
            for (unsigned idxNode0 = groupNodeIdx, idxEnd = (groupNodeIdx + dualNodeDistance); idxNode0 < idxEnd; ++idxNode0)
            {
                unsigned idxNode1 = idxNode0 + dualNodeDistance;
                // DIT radix-2 butterfly:
                // X[r]            = G[r] + H[r] * W^r
                // X[r + N/2]    = G[r] - H[r] * W^r
                //
                // with the following correspondences:
                // data[nodeOne] -> X[r]
                // data[nodeTwo] -> X[r + N/4]
                //
                // node1: add prod of node2 and twiddle factor.
                // node2: diff of node1 - prod of node2 and twiddle factor.
                Complex product(twiddleFactor * data[idxNode1]);
                data[idxNode1]  = data[idxNode0] - product;
                data[idxNode0] += product;

                // Calculate the next transform factor via trigonometric recurrence.
                if ((idxNode0 + 1) < (groupNodeIdx + dualNodeDistance))
                    twiddleFactor += twiddleMultiplier * twiddleFactor;
            }
        }
    };

    /** Specialization for case SampleCnt=4, direction=1 (forward).
        \param Complex ... The complex type.
    */
    template<typename Complex>
    class Radix2DIT<std::integral_constant<unsigned, 4>, std::integral_constant<int, 1>, Complex>
        : public SubTask<Radix2DIT<std::integral_constant<unsigned, 4>, std::integral_constant<int, 1>, Complex>,
                         Complex>
    {
    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned int groupNodeIdx = 0) const
        {
            // 1st stage butterfly between sequent nodes (distance: 1) - no need for twiddle factor multiplies, since 
            // twiddle factor is 1.
            // Nodes 0 and 1.
            Complex temp = data[groupNodeIdx + 1];
            data[groupNodeIdx + 1] = data[groupNodeIdx] - temp;
            data[groupNodeIdx] += temp;

            // 1st stage butterfly between sequent nodes (distance: 1) - no need for twiddle factor multiplies, since
            // twiddle factor is 1.
            // Nodes 2 and 3.
            temp = data[groupNodeIdx + 3];
            data[groupNodeIdx + 3] = data[groupNodeIdx + 2] - temp;
            data[groupNodeIdx + 2] += temp;

            // 2nd stage butterfly between sequent nodes (distance: 2) - no need for twiddle factor multiplies, since
            // twiddle factor is 1.
            // Nodes 0 and 2.
            temp = data[groupNodeIdx + 2];
            data[groupNodeIdx + 2] = data[groupNodeIdx] - temp;
            data[groupNodeIdx] += temp;

            // 2nd stage butterfly between sequent nodes (distance: 2) - multiplication with twiddle factor of -j 
            // results in a -90° rotation of the complex vector in the complex plane.
            // Nodes 1 and 3.
            temp.real(-data[groupNodeIdx + 3].imag());
            temp.imag(data[groupNodeIdx + 3].real());
            data[groupNodeIdx + 3] = data[groupNodeIdx + 1] - temp;
            data[groupNodeIdx + 1] += temp;
        }
    };

    /** Specialization for case SampleCnt=4, direction=-1 (backward).
        \param Complex ... The complex type.
    */
    template<typename Complex>
    class Radix2DIT<std::integral_constant<unsigned, 4>, std::integral_constant<int, -1>, Complex>
        : public SubTask<Radix2DIT<std::integral_constant<unsigned, 4>, std::integral_constant<int, -1>, Complex>,
                         Complex>
    {
    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned int groupNodeIdx = 0) const
        {
            // 1st stage butterfly between sequent nodes (distance: 1) - no need for twiddle factor multiplies, since 
            // twiddle factor is 1.
            // Nodes 0 and 1.
            Complex temp = data[groupNodeIdx + 1];
            data[groupNodeIdx + 1] = data[groupNodeIdx] - temp;
            data[groupNodeIdx] += temp;

            // 1st stage butterfly between sequent nodes (distance: 1) - no need for twiddle factor multiplies, since
            // twiddle factor is 1.
            // Nodes 2 and 3.
            temp = data[groupNodeIdx + 3];
            data[groupNodeIdx + 3] = data[groupNodeIdx + 2] - temp;
            data[groupNodeIdx + 2] += temp;

            // 2nd stage butterfly between sequent nodes (distance: 2) - no need for twiddle factor multiplies, since
            // twiddle factor is 1.
            // Nodes 0 and 2.
            temp = data[groupNodeIdx + 2];
            data[groupNodeIdx + 2] = data[groupNodeIdx] - temp;
            data[groupNodeIdx] += temp;

            // 2nd stage butterfly between sequent nodes (distance: 2) - multiplication with twiddle factor of -j 
            // results in a -90° rotation of the complex vector in the complex plane.
            // Nodes 1 and 3.
            temp.real(data[groupNodeIdx + 3].imag());
            temp.imag(-data[groupNodeIdx + 3].real());
            data[groupNodeIdx + 3] = data[groupNodeIdx + 1] - temp;
            data[groupNodeIdx + 1] += temp;

        }
    };

    /** Specialization for case SampleCnt=2.
        \param DirectionFactor ... Specifies the direction of the DFT (forward: 1, backward: -1)
        \param Complex ... The complex type.
    */
    template<typename DirectionFactor,
             typename Complex>
    class Radix2DIT<std::integral_constant<unsigned, 2>, DirectionFactor, Complex>
        : public SubTask<Radix2DIT<std::integral_constant<unsigned, 2>, DirectionFactor, Complex>,
                         Complex>
    {
    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned int groupNodeIdx = 0) const
        {
            // 1st stage butterfly between sequent nodes - no need for twiddle factor  multiplies, since twiddle 
            // factor is 1.
            Complex temp = data[groupNodeIdx + 1];
            data[groupNodeIdx + 1] = data[groupNodeIdx] - temp;
            data[groupNodeIdx] += temp;
        }
    };

    /** Specialization for case SampleCnt=1.
        \param DirectionFactor ... Specifies the direction of the DFT (forward: 1, backward: -1)
        \param Complex ... The complex type.
    */
    template<typename DirectionFactor,
             typename Complex>
    class Radix2DIT<std::integral_constant<unsigned, 1>, DirectionFactor, Complex>
        : public SubTask<Radix2DIT<std::integral_constant<unsigned, 1>, DirectionFactor, Complex>,
                         Complex>
    {
    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex*, unsigned int) const
        {}
    };


    /** Performs a radix 2 decimation in frequency FFT using template metaprogramming.
        \param SampleCnt ... The count of samples to be processed in this recursion level (stage)
        \param DirectionFactor ... Specifies the direction of the DFT (forward: 1, backward: -1)
        \param Complex ... The complex type.
    */
    template<typename SampleCnt,
             typename DirectionFactor,
             typename Complex>
    class Radix2DIF
        : public SubTask<Radix2DIF<SampleCnt, DirectionFactor, Complex>,
                         Complex>
    {
        Radix2DIF<std::integral_constant<unsigned, SampleCnt::value / 2>, DirectionFactor, Complex> recursionLevel_;

    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned groupNodeIdx = 0) const
        {
            // dualNodeDistance is the distance between elements (successive nodes) of a 
            // dual tuple, e.g. ..., 8, 4, 2, 1.
            unsigned dualNodeDistance = SampleCnt::value >> 1;

            // Create twiddle factor multiplier for trigonometric recurrence.
            Complex twiddleMultiplier(
                static_cast<typename Complex::value_type>(-2.0 * basic::sine<typename Complex::value_type>(1, SampleCnt::value) * basic::sine<typename Complex::value_type>(1, SampleCnt::value)),
                static_cast<typename Complex::value_type>(DirectionFactor::value * basic::sine<typename Complex::value_type>(2, SampleCnt::value)));
            // Create transform factor.
            Complex twiddleFactor(1.0, 0.0);

            // Run through dual nodes within the current group.
            for (unsigned idxNode0 = groupNodeIdx, idxEnd = (groupNodeIdx + dualNodeDistance); idxNode0 < idxEnd; ++idxNode0)
            {
                unsigned idxNode1 = idxNode0 + dualNodeDistance;
                // DIF radix-2 butterfly:
                // g[l] = x[l] + x[l + N/2]
                // h[l] = (x[l] - x[l + N/2]) * W^l
                //
                // with the following correspondences:
                // data[nodeOne] -> g[l]
                // data[nodeTwo] -> h[l]
                //
                // node1: sum of node1 and node2.
                // node2: (diff off node1 - node2) * twiddle factor.
                Complex sum(data[idxNode0] + data[idxNode1]);
                data[idxNode1] = (data[idxNode0] - data[idxNode1]) * twiddleFactor;
                data[idxNode0] = sum;

                // Calculate the next transform factor via trigonometric recurrence.
                if ((idxNode0 + 1) < (groupNodeIdx + dualNodeDistance))
                    twiddleFactor += twiddleMultiplier * twiddleFactor;
            }

            // Recursion goes down. Calculation starts in the last recursion stage with 2 nodes and goes up: 4, 8, ...
            recursionLevel_.apply(data, groupNodeIdx);
            recursionLevel_.apply(data, groupNodeIdx + dualNodeDistance);
        }
    };

    /** Specialization for case SampleCnt=4, direction=1 (forward).
        \param Complex ... The complex type.
    */
    template<typename Complex>
    class Radix2DIF<std::integral_constant<unsigned, 4>, std::integral_constant<int, 1>, Complex>
        : public SubTask<Radix2DIF<std::integral_constant<unsigned, 4>, std::integral_constant<int, 1>, Complex>,
                         Complex>
    {
    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned int groupNodeIdx = 0) const
        {
            // 1st stage butterfly between sequent nodes (distance: 2) - no need for twiddle factor multiplies, since
            // twiddle factor is 1.
            // Nodes 0 and 2.
            Complex temp = data[groupNodeIdx + 2];
            data[groupNodeIdx + 2] = data[groupNodeIdx] - temp;
            data[groupNodeIdx] += temp;

            // 1st stage butterfly between sequent nodes (distance: 2) - multiplication with twiddle factor of -j
            // results in a -90° rotation of the complex vector in the complex plane.
            // Nodes 1 and 3.
            temp = data[groupNodeIdx + 3];
            data[groupNodeIdx + 3].real(temp.imag() - data[groupNodeIdx + 1].imag());
            data[groupNodeIdx + 3].imag(data[groupNodeIdx + 1].real() - temp.real());
            data[groupNodeIdx + 1] += temp;

            // 2nd stage butterfly between sequent nodes (distance: 1) - no need for twiddle factor multiplies, since
            // twiddle factor is 1.
            // Nodes 0 and 1.
            temp = data[groupNodeIdx + 1];
            data[groupNodeIdx + 1] = data[groupNodeIdx] - temp;
            data[groupNodeIdx] += temp;

            // 2nd stage butterfly between sequent nodes (distance: 1) - no need for twiddle factor multiplies, since
            // twiddle factor is 1.
            // Nodes 2 and 3.
            temp = data[groupNodeIdx + 3];
            data[groupNodeIdx + 3] = data[groupNodeIdx + 2] - temp;
            data[groupNodeIdx + 2] += temp;
        }
    };

    /** Specialization for case SampleCnt=4, direction=-1 (backward).
        \param Complex ... The complex type.
    */
    template<typename Complex>
    class Radix2DIF<std::integral_constant<unsigned, 4>, std::integral_constant<int, -1>, Complex>
        : public SubTask<Radix2DIF<std::integral_constant<unsigned, 4>, std::integral_constant<int, -1>, Complex>,
                         Complex>
    {
    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned int groupNodeIdx = 0) const
        {
            // 1st stage butterfly between sequent nodes (distance: 2) - no need for twiddle factor multiplies, since
            // twiddle factor is 1.
            // Nodes 0 and 2.
            Complex temp = data[groupNodeIdx + 2];
            data[groupNodeIdx + 2] = data[groupNodeIdx] - temp;
            data[groupNodeIdx] += temp;

            // 1st stage butterfly between sequent nodes (distance: 2) - multiplication with twiddle factor of -j
            // results in a -90° rotation of the complex vector in the complex plane.
            // Nodes 1 and 3.
            temp = data[groupNodeIdx + 3];
            data[groupNodeIdx + 3].real(data[groupNodeIdx + 1].imag() - temp.imag());
            data[groupNodeIdx + 3].imag(temp.real() - data[groupNodeIdx + 1].real());
            data[groupNodeIdx + 1] += temp;

            // 2nd stage butterfly between sequent nodes (distance: 1) - no need for twiddle factor multiplies, since
            // twiddle factor is 1.
            // Nodes 0 and 1.
            temp = data[groupNodeIdx + 1];
            data[groupNodeIdx + 1] = data[groupNodeIdx] - temp;
            data[groupNodeIdx] += temp;

            // 2nd stage butterfly between sequent nodes (distance: 1) - no need for twiddle factor multiplies, since
            // twiddle factor is 1.
            // Nodes 2 and 3.
            temp = data[groupNodeIdx + 3];
            data[groupNodeIdx + 3] = data[groupNodeIdx + 2] - temp;
            data[groupNodeIdx + 2] += temp;
        }
    };

    /** Specialization for case SampleCnt=2.
        \param DirectionFactor ... Specifies the direction of the DFT (forward: 1, backward: -1)
        \param Complex ... The complex type.
    */
    template<typename DirectionFactor,
             typename Complex>
    class Radix2DIF<std::integral_constant<unsigned, 2>, DirectionFactor, Complex>
        : public SubTask<Radix2DIF<std::integral_constant<unsigned, 2>, DirectionFactor, Complex>,
                         Complex>
    {
    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned int groupNodeIdx = 0) const
        {
            // 1st stage butterfly between sequent nodes - no need for twiddle factor  multiplies, since twiddle 
            // factor is 1.
            Complex temp = data[groupNodeIdx + 1];
            data[groupNodeIdx + 1] = data[groupNodeIdx] - temp;
            data[groupNodeIdx] += temp;
        }
    };

    /** Specialization for case SampleCnt=1.
        \param DirectionFactor ... Specifies the direction of the DFT (forward: 1, backward: -1)
        \param Complex ... The complex type.
    */
    template<typename DirectionFactor,
             typename Complex>
    class Radix2DIF<std::integral_constant<unsigned, 1>, DirectionFactor, Complex>
        : public SubTask<Radix2DIF<std::integral_constant<unsigned, 1>, DirectionFactor, Complex>,
                         Complex>
    {
    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex*, unsigned int) const
        {}
    };
}