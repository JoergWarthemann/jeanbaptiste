#pragma once

#include "../basic/SineCosine.h"
#include <boost/math/constants/constants.hpp>
#include <complex>
#include "../SubTask.h"

namespace jeanbaptiste::core
{
    /** Performs a radix 4 decimation in time FFT using template metaprogramming.
        \param SampleCnt ... The count of samples to be processed in this recursion level (stage)
        \param DirectionFactor ... Specifies the direction of the DFT (forward: 1, backward: -1)
        \param Complex ... The complex type.
    */
    template<typename SampleCnt,
             typename DirectionFactor,
             typename Complex>
    class Radix4DIT
        : public SubTask<Radix4DIT<SampleCnt, DirectionFactor, Complex>,
                         Complex>
    {
        Radix4DIT<std::integral_constant<unsigned, SampleCnt::value / 4>, DirectionFactor, Complex> recursionLevel_;

    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned groupNodeIdx = 0) const
        {
            using ValueType = typename Complex::value_type;

            // quaternaryNodeDistance is the distance between elements (successive nodes) of a
            // quaternary tuple, e.g. ... 16, 4, 1.
            auto quaternaryNodeDistance = SampleCnt::value >> 2;
                
            // Recursion goes down. Calculation starts in the last recursion stage with 4 nodes and goes up: 16, 64, ... .
            recursionLevel_.apply(data, groupNodeIdx);
            recursionLevel_.apply(data, groupNodeIdx + quaternaryNodeDistance);
            recursionLevel_.apply(data, groupNodeIdx + quaternaryNodeDistance * 2);
            recursionLevel_.apply(data, groupNodeIdx + quaternaryNodeDistance * 3);

            // Create twiddle factor multiplier for trigonometric recurrence.
            constexpr Complex twiddleMultiplier(
                static_cast<ValueType>(-2.0 * basic::sine<ValueType>(1.0 / SampleCnt::value * constants::pi<ValueType>()) * basic::sine<ValueType>(1.0 / SampleCnt::value * constants::pi<ValueType>())),
                static_cast<ValueType>(DirectionFactor::value * basic::sine<ValueType>(2.0 / SampleCnt::value * constants::pi<ValueType>())));
            // Create transform factor.
            Complex twiddleFactor(1.0, 0.0);

            // Run through quaternary tuples of the current group. idxNode0 flags the tuple start.
            for (auto idxNode0 = groupNodeIdx, idxEnd = (groupNodeIdx + quaternaryNodeDistance); idxNode0 < idxEnd; ++idxNode0)
            {
                // Create twiddle factors.
                ValueType temp = 1.5 - 0.5 * (twiddleFactor.real() * twiddleFactor.real() + twiddleFactor.imag() * twiddleFactor.imag());
                Complex wn4(twiddleFactor.real() * temp, twiddleFactor.imag() * temp);
                Complex wn2 = wn4 * wn4;
                Complex w3n4 = wn2 * wn4;

                // Create tuple node indexes.
                auto idxNode1 = idxNode0 + quaternaryNodeDistance;
                auto idxNode2 = idxNode1 + quaternaryNodeDistance;
                auto idxNode3 = idxNode2 + quaternaryNodeDistance;

                // DIT radix-4 butterfly:
                // X[r]            = (Y[r] + G[r] * W^2r) +  (Z[r] * W^r + H[r] * W^3r)
                // X[r + N/4]    = (Y[r] - G[r] * W^2r) - j(Z[r] * W^r - H[r] * W^3r)
                // X[r + N/2]    = (Y[r] + G[r] * W^2r) -  (Z[r] * W^r + H[r] * W^3r)
                // X[r + 3N/4]    = (Y[r] - G[r] * W^2r) + j(Z[r] * W^r - H[r] * W^3r)
                //
                // with the following correspondences:
                // data[idxNode0] -> X[r]
                // data[idxNode1] -> X[r + N/4]
                // data[idxNode2] -> X[r + N/2]
                // data[idxNode3] -> X[r + 3N/4]
                //
                // Attention: In order to use bit reversal instead of digit reversal this original formula has been 
                //            customized slightly (Sidney Burrus).

                // Save temporary results.
                Complex temp1 = data[idxNode0];
                Complex temp2 = data[idxNode2] * wn4;
                Complex temp3 = data[idxNode1] * wn2;
                Complex temp4 = data[idxNode3] * w3n4;

                // Conduct radix-4 DIT butterfly.
                data[idxNode0] = (temp1 + temp3) + (temp2 + temp4);
                data[idxNode1] = (temp1 - temp3) + Complex(0, DirectionFactor::value) * (temp2 - temp4);
                data[idxNode2] = (temp1 + temp3) - (temp2 + temp4);
                data[idxNode3] = (temp1 - temp3) - Complex(0, DirectionFactor::value) * (temp2 - temp4);

                //data[idxNode0] = (temp1 + temp3) + (temp2 + temp4);
                //data[idxNode1] = (temp1 - temp3) - Complex(0, DirectionFactor::value) * (temp2 - temp4);
                //data[idxNode2] = (temp1 + temp3) - (temp2 + temp4);
                //data[idxNode3] = (temp1 - temp3) + Complex(0, DirectionFactor::value) * (temp2 - temp4);

                if ((idxNode0 + 1) < idxEnd)
                    twiddleFactor += twiddleMultiplier * twiddleFactor;
            }
        }
    };

    /** Specialization for case SampleCnt=4, direction=1 (forward).
        \param Complex ... The complex type.
    */
    template<typename Complex>
    class Radix4DIT<std::integral_constant<unsigned, 4>, std::integral_constant<int, 1>, Complex>
        : public SubTask<Radix4DIT<std::integral_constant<unsigned, 4>, std::integral_constant<int, 1>, Complex>,
                         Complex>
    {
    public:
        void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned int groupNodeIdx = 0) const
        {
            // 1st stage butterfly between sequent nodes - no need for twiddle factor multiplies, since twiddle factors 
            // are 1 when N = 4.
            // Multiplication with j results in a 90° rotation of the complex vector in the complex plane.
            // So we just need need 2 * (N = 4) = 8 complex additions.

            using ValueType = typename Complex::value_type;

            // Create tuple node indexes.
            auto idxNode0 = groupNodeIdx;
            auto idxNode1 = idxNode0 + 1;
            auto idxNode2 = idxNode1 + 1;
            auto idxNode3 = idxNode2 + 1;

            // Temporary results.
            ValueType tmpReal1 = data[idxNode0].real() + data[idxNode1].real();
            ValueType tmpReal2 = data[idxNode2].real() + data[idxNode3].real();
            ValueType tmpImag1 = data[idxNode0].imag() + data[idxNode1].imag();
            ValueType tmpImag2 = data[idxNode2].imag() + data[idxNode3].imag();

            ValueType tmpReal3 = data[idxNode0].real() - data[idxNode1].real();
            ValueType tmpReal4 = data[idxNode2].imag() - data[idxNode3].imag();
            ValueType tmpImag3 = data[idxNode0].imag() - data[idxNode1].imag();
            ValueType tmpImag4 = data[idxNode2].real() - data[idxNode3].real();

            data[idxNode0].real(tmpReal1 + tmpReal2);
            data[idxNode0].imag(tmpImag1 + tmpImag2);

            data[idxNode1].real(tmpReal3 - tmpReal4);
            data[idxNode1].imag(tmpImag3 + tmpImag4);

            data[idxNode2].real(tmpReal1 - tmpReal2);
            data[idxNode2].imag(tmpImag1 - tmpImag2);

            data[idxNode3].real(tmpReal3 + tmpReal4);
            data[idxNode3].imag(tmpImag3 - tmpImag4);
        }
    };

    /** Specialization for case SampleCnt=4, direction=-1 (backward).
        \param Complex ... The complex type.
    */
    template<typename  Complex>
    class Radix4DIT<std::integral_constant<unsigned, 4>, std::integral_constant<int, -1>, Complex>
		: public SubTask<Radix4DIT<std::integral_constant<unsigned, 4>, std::integral_constant<int, -1>, Complex>,
                         Complex>
    {
    public:
		void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned int groupNodeIdx = 0) const
        {
            // 1st stage butterfly between sequent nodes - no need for twiddle factor multiplies, since twiddle factors 
            // are 1 when N = 4.
            // Multiplication with j results in a 90� rotation of the complex vector in the complex plane.
            // So we just need need 2 * (N = 4) = 8 complex additions.

            using ValueType = typename Complex::value_type;

            // Create tuple node indexes.
            auto idxNode0 = groupNodeIdx;
            auto idxNode1 = idxNode0 + 1;
            auto idxNode2 = idxNode1 + 1;
            auto idxNode3 = idxNode2 + 1;

            // Temporary results.
            ValueType tmpReal1 = data[idxNode0].real() + data[idxNode1].real();
            ValueType tmpReal2 = data[idxNode2].real() + data[idxNode3].real();
            ValueType tmpImag1 = data[idxNode0].imag() + data[idxNode1].imag();
            ValueType tmpImag2 = data[idxNode2].imag() + data[idxNode3].imag();

            ValueType tmpReal3 = data[idxNode0].real() - data[idxNode1].real();
            ValueType tmpReal4 = data[idxNode2].imag() - data[idxNode3].imag();
            ValueType tmpImag3 = data[idxNode0].imag() - data[idxNode1].imag();
            ValueType tmpImag4 = data[idxNode2].real() - data[idxNode3].real();

            data[idxNode0].real(tmpReal1 + tmpReal2);
            data[idxNode0].imag(tmpImag1 + tmpImag2);

            data[idxNode1].real(tmpReal3 + tmpReal4);
            data[idxNode1].imag(tmpImag3 - tmpImag4);

            data[idxNode2].real(tmpReal1 - tmpReal2);
            data[idxNode2].imag(tmpImag1 - tmpImag2);

            data[idxNode3].real(tmpReal3 - tmpReal4);
            data[idxNode3].imag(tmpImag3 + tmpImag4);
        }
    };

    /** Specialization for case SampleCnt=1.
        \param DirectionFactor ... Specifies the direction of the DFT (forward: 1, backward: -1)
        \param Complex ... The complex type.
    */
    template<typename DirectionFactor,
             typename Complex>
    class Radix4DIT<std::integral_constant<unsigned, 1>, DirectionFactor, Complex>
    {
    public:
		void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex*, unsigned int) const
        {}
    };

    /** Performs a radix 4 decimation in frequency FFT using template metaprogramming.
        \param SampleCnt ... The count of samples to be processed in this recursion level (stage)
        \param DirectionFactor ... Specifies the direction of the DFT (forward: 1, backward: -1)
        \param Complex ... The complex type.
    */
    template<typename SampleCnt,
             typename DirectionFactor,
             typename Complex>
    class Radix4DIF
		: public SubTask<Radix4DIF<SampleCnt, DirectionFactor, Complex>,
                         Complex>
    {
        Radix4DIF<std::integral_constant<unsigned, SampleCnt::value / 4>, DirectionFactor, Complex> recursionLevel_;

    public:
		void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned groupNodeIdx = 0) const
        {
            using ValueType = typename Complex::value_type;

            // quaternaryNodeDistance is the distance between elements (successive nodes) of a
            // quaternary tuple, e.g. ... 16, 4, 1.
            auto quaternaryNodeDistance = SampleCnt::value >> 2;

            // Create twiddle factor multiplier for trigonometric recurrence.
            constexpr Complex twiddleMultiplier(
				static_cast<ValueType>(-2.0 * basic::sine<ValueType>(1.0 / SampleCnt::value * constants::pi<ValueType>()) * basic::sine<ValueType>(1.0 / SampleCnt::value * constants::pi<ValueType>())),
                static_cast<ValueType>(DirectionFactor::value * basic::sine<ValueType>(2.0 / SampleCnt::value * constants::pi<ValueType>())));
            // Create transform factor.
            Complex twiddleFactor(1.0, 0.0);

            // Run through quaternary tuples of the current group. idxNode0 flags the tuple start.
            for (auto idxNode0 = groupNodeIdx, idxEnd = (groupNodeIdx + quaternaryNodeDistance); idxNode0 < idxEnd; ++idxNode0)
            {
                // Create twiddle factors.
                ValueType temp = 1.5 - 0.5 * (twiddleFactor.real() * twiddleFactor.real() + twiddleFactor.imag() * twiddleFactor.imag());
                Complex wn4(twiddleFactor.real() * temp, twiddleFactor.imag() * temp);
                Complex wn2 = wn4 * wn4;
                Complex w3n4 = wn2 * wn4;

                // Create tuple node indexes.
                auto idxNode1 = idxNode0 + quaternaryNodeDistance;
                auto idxNode2 = idxNode1 + quaternaryNodeDistance;
                auto idxNode3 = idxNode2 + quaternaryNodeDistance;

                // DIF radix-4 butterfly:
                // y[l] =  (x[l] + x[l + N/2]) +  (x[l + N/4] + x[l + 3N/4])
                // z[l] = ((x[l] - x[l + N/2]) - j(x[l + N/4] - x[l + 3N/4])) * W^l
                // g[l] = ((x[l] + x[l + N/2]) -  (x[l + N/4] + x[l + 3N/4])) * W^2l
                // h[l] = ((x[l] - x[l + N/2]) + j(x[l + N/4] - x[l + 3N/4])) * W^3l
                //
                // with the following correspondences:
                // data[idxNode0] -> y[l]
                // data[idxNode1] -> z[l]
                // data[idxNode2] -> g[l]
                // data[idxNode3] -> h[l]
                //
                // Attention: In order to use bit reversal instead of digit reversal this original formula has been 
                //            customized slightly (Sidney Burrus).

                // Save temporary results.
                Complex temp1 = data[idxNode0] + data[idxNode2];
                Complex temp2 = data[idxNode0] - data[idxNode2];
                Complex temp3 = data[idxNode1] + data[idxNode3];
                Complex temp4 = data[idxNode1] - data[idxNode3];

                // Conduct radix-4 DIF butterfly.
                data[idxNode0] =         temp1 + temp3;
                data[idxNode1] = wn2  * (temp1 - temp3);
                data[idxNode2] = wn4  * (temp2 + Complex(0, DirectionFactor::value) * temp4);
                data[idxNode3] = w3n4 * (temp2 - Complex(0, DirectionFactor::value) * temp4);

                //data[idxNode0] =         temp1 + temp3;
                //data[idxNode1] = wn4  * (temp2 - Complex(0, DirectionFactor::value) * temp4);
                //data[idxNode2] = wn2  * (temp1 - temp3);
                //data[idxNode3] = w3n4 * (temp2 + Complex(0, DirectionFactor::value) * temp4);

                if ((idxNode0 + 1) < idxEnd)
                    twiddleFactor += twiddleMultiplier * twiddleFactor;
            }

            // Recursion goes down. Calculation starts in the last recursion stage with 4 nodes and goes up: 16, 64, ... .
            recursionLevel_.apply(data, groupNodeIdx);
            recursionLevel_.apply(data, groupNodeIdx + quaternaryNodeDistance);
            recursionLevel_.apply(data, groupNodeIdx + quaternaryNodeDistance * 2);
            recursionLevel_.apply(data, groupNodeIdx + quaternaryNodeDistance * 3);
        }
    };

    /** Specialization for case SampleCnt=4, direction=1 (forward).
        \param Complex ... The complex type.
    */
    template<typename Complex>
    class Radix4DIF<std::integral_constant<unsigned, 4>, std::integral_constant<int, 1>, Complex>
		: public SubTask<Radix4DIF<std::integral_constant<unsigned, 4>, std::integral_constant<int, 1>, Complex>,
                         Complex>
    {
    public:
		void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned int groupNodeIdx = 0) const
        {
            // 1st stage butterfly between sequent nodes - no need for twiddle factor multiplies, since twiddle factors 
            // are 1 when N = 4.
            // Multiplication with j results in a 90� rotation of the complex vector in the complex plane.
            // So we just need need 2 * (N = 4) = 8 complex additions.

            using ValueType = typename Complex::value_type;

            // Create tuple node indexes.
            auto idxNode0 = groupNodeIdx;
            auto idxNode1 = idxNode0 + 1;
            auto idxNode2 = idxNode1 + 1;
            auto idxNode3 = idxNode2 + 1;

            // Temporary results.
            ValueType tmpReal1 = data[idxNode0].real() + data[idxNode2].real();
            ValueType tmpReal2 = data[idxNode1].real() + data[idxNode3].real();
            ValueType tmpImag1 = data[idxNode0].imag() + data[idxNode2].imag();
            ValueType tmpImag2 = data[idxNode1].imag() + data[idxNode3].imag();

            ValueType tmpReal3 = data[idxNode0].real() - data[idxNode2].real();
            ValueType tmpReal4 = data[idxNode1].imag() - data[idxNode3].imag();
            ValueType tmpImag3 = data[idxNode0].imag() - data[idxNode2].imag();
            ValueType tmpImag4 = data[idxNode1].real() - data[idxNode3].real();

            data[idxNode0].real(tmpReal1 + tmpReal2);
            data[idxNode0].imag(tmpImag1 + tmpImag2);

            data[idxNode1].real(tmpReal1 - tmpReal2);
            data[idxNode1].imag(tmpImag1 - tmpImag2);

            data[idxNode2].real(tmpReal3 - tmpReal4);
            data[idxNode2].imag(tmpImag3 + tmpImag4);

            data[idxNode3].real(tmpReal3 + tmpReal4);
            data[idxNode3].imag(tmpImag3 - tmpImag4);
        }
    };

    /** Specialization for case SampleCnt=4, direction=-1 (backward).
        \param Complex ... The complex type.
    */
    template<typename Complex>
    class Radix4DIF<std::integral_constant<unsigned, 4>, std::integral_constant<int, -1>, Complex>
		: public SubTask<Radix4DIF<std::integral_constant<unsigned, 4>, std::integral_constant<int, -1>, Complex>,
                         Complex>
    {
    public:
		void operator()(Complex* data) const
        {
            apply(data);
        }

        void apply(Complex* data, unsigned int groupNodeIdx = 0) const
        {
            // 1st stage butterfly between sequent nodes - no need for twiddle factor multiplies, since twiddle factors 
            // are 1 when N = 4.
            // Multiplication with j results in a 90� rotation of the complex vector in the complex plane.
            // So we just need need 2 * (N = 4) = 8 complex additions.

            using ValueType = typename Complex::value_type;

            // Create tuple node indexes.
            auto idxNode0 = groupNodeIdx;
            auto idxNode1 = idxNode0 + 1;
            auto idxNode2 = idxNode1 + 1;
            auto idxNode3 = idxNode2 + 1;

            // Temporary results.
            ValueType tmpReal1 = data[idxNode0].real() + data[idxNode2].real();
            ValueType tmpReal2 = data[idxNode1].real() + data[idxNode3].real();
            ValueType tmpImag1 = data[idxNode0].imag() + data[idxNode2].imag();
            ValueType tmpImag2 = data[idxNode1].imag() + data[idxNode3].imag();

            ValueType tmpReal3 = data[idxNode0].real() - data[idxNode2].real();
            ValueType tmpReal4 = data[idxNode1].imag() - data[idxNode3].imag();
            ValueType tmpImag3 = data[idxNode0].imag() - data[idxNode2].imag();
            ValueType tmpImag4 = data[idxNode1].real() - data[idxNode3].real();

            data[idxNode0].real(tmpReal1 + tmpReal2);
            data[idxNode0].imag(tmpImag1 + tmpImag2);

            data[idxNode1].real(tmpReal1 - tmpReal2);
            data[idxNode1].imag(tmpImag1 - tmpImag2);

            data[idxNode2].real(tmpReal3 + tmpReal4);
            data[idxNode2].imag(tmpImag3 - tmpImag4);

            data[idxNode3].real(tmpReal3 - tmpReal4);
            data[idxNode3].imag(tmpImag3 + tmpImag4);
        }
    };

    /** Specialization for case SampleCnt=1.
        \param DirectionFactor ... Specifies the direction of the DFT (forward: 1, backward: -1)
        \param Complex ... The complex type.
    */
    template<typename DirectionFactor,
             typename Complex>
    class Radix4DIF<std::integral_constant<unsigned, 1>, DirectionFactor, Complex>
		: public SubTask<Radix4DIF<std::integral_constant<unsigned, 1>, DirectionFactor, Complex>,
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