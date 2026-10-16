#include <gmock/gmock.h>
#include <gtest/gtest.h>

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

// #define BOOST_TEST_MODULE jeanbaptiste tests
// #define BOOST_TEST_DYN_LINK

// #include "FixtureRadix2.cpp"
// #include "FixtureRadix4.cpp"
// #include "FixtureRadixSplit24.cpp"
// #include "FixtureBitReversal.cpp"
// #include "FixtureSinCos.cpp"
// #include "FixtureWindowCalculation.cpp"
// #include "FixtureFft.cpp"