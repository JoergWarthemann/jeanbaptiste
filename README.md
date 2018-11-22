# JeanBaptiste

JeanBaptiste is an easy to use implementation of FFT (Fast Fourier Transform). It intensively makes use of template metaprogramming and constant expressions to move as much of the necessary computations as possible into compilation time.

In order to make using compile time generated types easier and to give back some flexibility to the user boost::hana is used.

## Features

* different cases of the Cooley-Tukey algorithm
    * radix-2
    * radix-4
    * split-radix-2-4
* different cases of normalization
    * square root
    * division by sample count
* different types of window functions
    * Bartlett
    * Hamming
    * von Hann
    * ...
* runtime selection of transform length

## Implementation

All transform algorithms are implemented using their recursive definition based on template metaprogramming. Transform lengths are compile time constants. Although flexibility is given to the user by allowing him to select the right one from a boost::hana powered factory at execution time.

Partially computing twiddle factors at compilation time reduces the calculation cost at runtime.

Bit reversal indexing and window functions are implemented using constant expressions. This again reduces the calculation cost at execution time.

## Compilation

cmake, VSCode

## Usage

## What's next?