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

All transform algorithms are implemented using their recursive definition and based on template metaprogramming. Transform lengths are compile time constants. They are used to partially compute twiddle factors which reduces the calculation cost at runtime. Despite the limitation by defining the transform length at compile time, runtime flexibility is given to the user by allowing him to select the right transform length from a [boost::hana](https://www.boost.org/doc/libs/1_68_0/libs/hana/doc/html/index.html) powered factory.

Bit reversal indexing and window functions are implemented using constant expressions. This again reduces the calculation cost at execution time.

## Requirements

JeanBaptiste requires a C++17 aware compiler. It uses [boost::hana](https://www.boost.org/doc/libs/1_68_0/libs/hana/doc/html/index.html).  
The tests are implemented using [boost::test](https://www.boost.org/doc/libs/1_68_0/libs/test/doc/html/index.html).

## Compilation

JeanBaptiste comes with a [CMake](https://cmake.org) build script. CMake works by generating native Makefiles or build projects which can be used in the particular environment.

## Usage

First of all an algorithm factory needs to be created.

```cpp
AlgorithmFactory<  
    Begin,
    End,
    Radix,
    Decimation,
    Direction,
    Window,
    Normalization,
    Complex> algorithmFactory;
```

* `Begin` and `End` define the range of FFT stages for the factory. If runtime transform sample counts of 1024, 2048 and 4096 are expected in a radix-2 use case, `Begin` and `End` should be chosen as 10 and 12. Where 2^stage results into the actual sample count.
* `Radix` defines the radix of the used algorithm. Options: `Radix_2`, `Radix_4`, `Radix_Split_2_4`
* `Decimation` defines the FFT type: `Decimation_In_Time` and `Decimation_In_Frequency`.
* `Direction` defines whether to run a FFT (`Direction_Forward`) or an inverse FFT (`Direction_Backward`).
* `Window` defines whether to use a windowing function before running the actual FFT algorithm. Options: `Window_None`,  `Window_Bartlett`, `Window_BlackmanHarris`, `Window_Blackman`, `Window_Cosine`, `Window_FlatTop`, `Window_Hamming`, `Window_vonHann`, `Window_Welch`
* `Normalization` defines whether a normalization is to be used. Options: `Normalization_No`, `Normalization_Division_By_Length` (result is normalized by a factor of 1/N), `Normalization_Square_Root` (result is normalized by a factor of 1/√N)
* `Complex` defines the type of complex number which is to be used.

Then the factory can be instructed to create the specified algorithm for the desired sample count.

```cpp
auto algorithm = algorithmFactory.getAlgorithm(11);
```

Then sample data is handed over and the results are calculated in place.

```cpp
algorithm->operator()(&sampleData[0]);
```

## Further development

* integrate the real FFT algorithm
* make use of expression templates to remove temporaries and increase efficiency
