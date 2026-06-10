# Jean-baptiste library development

## Next steps

* [x] use native cmake
  `/usr/bin/cmake`
* [x] use native tools (compiler, stdlib)
* [x] use natively installed gtest -> this is not feasable due to problems caused by a proprietary SDK - deactivate
      proprietary SDK in cmake-kits.json ion this case
* [x] unuse proprietary toolchain -> need to add toolchain in order to use the SDKs googletest - deactivate
      proprietary SDK in cmake-kits.json ion this case
* [x] enable use of C++23
* [x] write tests for range reduction
* [x] use modern C++ concept for creation and control of work data in bit reversal test (ranges?)
* [ ] consider working in all algorithms on std::span, i.e. forward std::span only?
* [ ] Replace all #pragma preprocessor commands by #ifdef
* [x] Rewrite SubTask::operator() to use std::span (rewrite operator() as template function)
* [x] modernize Radix2 and Radix2Test
* [x] make Radix2Test go through all relevant test files
* [ ] modernize Radix4 and Radix4Test
* [ ] modernize SplitRadix and SplitRadixTest
* [x] rewrite AlgorithmFixture.hpp using fold expressions and use it in Radix2Test.cpp
* [x] rewrite AlgorithmResultAnalysis.hpp and use it in Radix2Test.cpp to load sample data from files (use ranges)
* [x] rewrite ExecutableAlgorithm.hpp (prefer = default for destructor?)
* [ ] override virtual base class destructors and make them default at least

* [x] use std::span for sample range in SubTask::operator()
* [x] update TestCaseLoader to use std::string_view, std::filesystem
* [x] use std::format instead of boost::format
* [x] use std::span for sample range in Radix2::operator()
* [ ] use concepts
* [x] turn mWorkingSet, mExpectedOutFFT and mExpectedOutIFFT into mInput and mOutput
* [ ] make Algorithm, AlgorithmFactory, SubTask usable to instantiate Radix2 for tests
* [x] replace namespace name "jeanbaptiste" by shorter "jb"
* [ ] unify usage of std::numbers or the more complete boost::math::constants for pi
* [x] use clangd
* [x] update window types

--> rebuild Algorithm.hpp, then AlgorithmFactory - ignore Windowing for now
--> update Radix4 technically like Radix2
--> add radix-4 and split-radix cases to Algorithm

------------------------------------------------------------------------------------------------------------------------

std::vector<float> audio_samples = {1.0f, 2.0f, 3.0f, 4.0f};

auto complex_view = audio_samples
                    | std::views::transform( [&](const float& real){
                        return std::complex<float>(real, 0.0f);
                      });

for (const auto& c : complex_view) {
    std::cout << c << '\n';
}

Lazy Evaluation: The view doesn’t store the transformed elements. Instead, it stores a reference to the original range (audio_samples) and the transformation function (which converts a float to std::complex<float>).
On-Demand Transformation: When you iterate over the view, the transformation function is applied to each element of the original range as you access it. This means that a new std::complex<float> object is created each time you access an element of the view.

So, yes, new std::complex<float> objects are indeed created during iteration. However, these objects are created only when needed and are not stored in a separate container. This approach minimizes memory usage and avoids unnecessary allocations.

------------------------------------------------------------------------------------------------------------------------

#include <vector>
#include <complex>
#include <ranges>
#include <algorithm>

int main() {
    std::vector<float> audio_samples = {1.0f, 2.0f, 3.0f, 4.0f};

    // Create a vector to store the complex numbers
    std::vector<std::complex<float>> complex_samples(audio_samples.size());

    // Transform and store the results in the complex_samples vector
    std::transform(audio_samples.begin(), audio_samples.end(),
                   complex_samples.begin(),
                   [](const float& real) {
                       return std::complex<float>(real, 0.0f);
                   });

    // Create a view on the complex_samples vector
    auto complex_view = complex_samples | std::views::all;

    // Now you can use complex_view as needed
    for (const auto& complex_sample : complex_view) {
        // Use complex_sample here
    }

    return 0;
}
Summary
By transforming the audio samples once and storing the results in a std::vector, you avoid the overhead of creating temporary objects repeatedly.
You can still utilize views for iteration or further processing without incurring additional transformation costs. This approach strikes a balance between using modern C++ features and maintaining performance.