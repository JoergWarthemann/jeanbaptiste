#include <array>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "tools/RingBuffer.hpp"

namespace jb::testing {

class RingBufferInsertLValuesTest
    : public ::testing::Test {
protected:
    static constexpr std::size_t sourceSize{20};
    static constexpr std::size_t destinationSize{19};

    using ValueType = std::string;
    using SourceArrayType = std::array<ValueType, sourceSize>;
    using DestinationArrayType = tools::RingBuffer<ValueType>;

    SourceArrayType source{
        "one", "two", "three", "four", "five",
        "six", "seven", "eight", "nine", "ten",
        "eleven", "twelve", "thirteen", "fourteen", "fivteen",
        "sixteen", "seventeen", "eighteen", "nineteen", "twenty"};
};

TEST_F(RingBufferInsertLValuesTest, CopiesSingleElementsWithoutMovingFromSource)
{
    DestinationArrayType destination(destinationSize);

    for (std::size_t i{0}; i < sourceSize; ++i) {
        destination.insert(source[i]);

        EXPECT_FALSE(source[i].empty());
        EXPECT_EQ(std::min(destination.capacity(), i + 1), destination.currentSize());
    }
}

TEST_F(RingBufferInsertLValuesTest, DoesNotExceedCapacityWhenInsertingOversizedBlock)
{
    DestinationArrayType destination(destinationSize);

    destination.insert(source, sourceSize);

    EXPECT_EQ(destinationSize, destination.currentSize());
}

TEST_F(RingBufferInsertLValuesTest, CopiesNewestElementsWhenInsertingOversizedBlock)
{
    DestinationArrayType destination(destinationSize);

    destination.insert(source, sourceSize);

    for (std::size_t i{0}; i < destinationSize; ++i) {
        EXPECT_EQ(source[sourceSize - i - 1], destination.copy(i));
    }
}

TEST_F(RingBufferInsertLValuesTest, CopiesNewestElementsWhenInsertingDifferentSizedBlocks)
{
    DestinationArrayType destination(destinationSize);
    SourceArrayType source1 = source;
    SourceArrayType source2 = source;
    SourceArrayType source3 = source;

    constexpr std::size_t size1{15};
    constexpr std::size_t size2{10};
    constexpr std::size_t size3{5};
    std::vector<ValueType> controlGroup;
    std::copy(source1.begin(), source1.begin() + size1, std::back_inserter(controlGroup));
    std::copy(source2.begin(), source2.begin() + size2, std::back_inserter(controlGroup));
    std::copy(source3.begin(), source3.begin() + size3, std::back_inserter(controlGroup));

    destination.insert(source1, size1);
    destination.insert(source2, size2);
    destination.insert(source3, size3);

    auto currentlyExpected = controlGroup.end();
    for (std::size_t i{0}; i < destinationSize; ++i) {
        --currentlyExpected;
        EXPECT_EQ(*currentlyExpected, destination.copy(i));
    }
}

} // namespace jb::testing