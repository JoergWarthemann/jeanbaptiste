# JeanBaptiste

JeanBaptiste is an easy to use implementation of FFT (Fast Fourier Transform). It intensively makes use of template metaprogramming and constant expressions to move as much of the necessary computations as possible into compilation time.

## Features

* different cases of the Cooley-Tukey algorithm
    * radix-2
    * radix-4
    * split-radix-2-4
* square root normalization
* different types of window functions (e.g. Bartlett, Hamming, von Hann)
* runtime selection of transform length

## Implementation

All transform algorithms are implemented using their recursive definition based on template metaprogramming. Transform lengths are compile time constants. Although flexibility is given to the user by allowing him to select the right one from a boost::hana powered factory at execution time.

Partially computing twiddle factors at compilation time reduces the calculation cost at runtime.

Bit reversal indexing and window functions are implemented using constant expressions. This again reduces the calculation cost at execution time.

## Compilation

cmake, VSCode

## Usage

## What's next?