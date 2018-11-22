#pragma once

namespace jeanbaptiste::options
{
    struct Direction_Forward {};
    struct Direction_Backward {};
    struct Decimation_In_Frequency {};
    struct Decimation_In_Time {};
    struct Normalization_No {};
    struct Normalization_Division_By_Length {};
    struct Normalization_Square_Root {};
    struct Radix_2 {};
    struct Radix_4 {};
    struct Radix_Split_2_4 {};
    struct Window_None {};
    struct Window_Bartlett {};
    struct Window_BlackmanHarris {};
    struct Window_Blackman {};
    struct Window_Cosine {};
    struct Window_FlatTop {};
    struct Window_Hamming {};
    struct Window_vonHann {};
    struct Window_Welch {};
}