#ifndef JB_OPTIONS_HPP_
#define JB_OPTIONS_HPP_

namespace jb::options {
struct Direction_Forward {}; // Forward transform direction (analysis)
struct Direction_Backward {}; // Backward transform direction (synthesis)
struct Data_Complex {}; // Complex-valued data
struct Data_Real {}; // Real-valued data
struct Decimation_In_Frequency {}; // Decimation in frequency (DIF) FFT
struct Decimation_In_Time {}; // Decimation in time (DIT) FFT
struct Normalization_No {}; // No normalization
struct Normalization_Division_By_Length {}; // Normalize by dividing by the FFT length
struct Normalization_Square_Root {}; // Normalize by the square root of the FFT length
struct Radix_2 {}; // Radix-2 FFT
struct Radix_4 {}; // Radix-4 FFT
struct Radix_Split_2_4 {}; // Split-radix FFT
struct Window_None {}; // No window
struct Window_Bartlett {}; // Bartlett window
struct Window_BlackmanHarris {}; // Blackman-Harris window
struct Window_Blackman {}; // Blackman window
struct Window_Cosine {}; // Cosine window
struct Window_FlatTop {}; // Flat-top window
struct Window_Hamming {}; // Hamming window
struct Window_vonHann {}; // von Hann (Hann) window
struct Window_Welch {}; // Welch window
struct Window_Before {}; // Apply window before the transform (analysis)
struct Window_After {}; // Apply window after the transform (synthesis)
} // namespace jb::options
#endif // JB_OPTIONS_HPP_