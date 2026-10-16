#include <array>
#include <string>

#include <gtest/gtest.h>

#include "tools/RingBuffer.hpp"

namespace jb::testing {

class RingBufferCapacityTest
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

TEST_F(RingBufferCapacityTest, DoesNotExceedCapacityWhenOversizedBlockIsInserted)
{
    DestinationArrayType destination(destinationSize);

    destination.insert(source, sourceSize);

    EXPECT_EQ(destinationSize, destination.currentSize());
}

TEST_F(RingBufferCapacityTest, ClearsCurrentSizeWhenResetting)
{
    DestinationArrayType destination(destinationSize);
    destination.insert(source, sourceSize);

    destination.reset();

    EXPECT_EQ(0U, destination.currentSize());
}

TEST_F(RingBufferCapacityTest, ClearsPartiallyFilledBufferWhenResetting)
{
    DestinationArrayType destination(destinationSize);
    constexpr std::size_t insertedElements{3};

    destination.insert(source, insertedElements);
    destination.reset();

    EXPECT_EQ(destinationSize, destination.capacity());
    EXPECT_EQ(0U, destination.currentSize());
    EXPECT_EQ(ValueType{}, destination.copy(0));
}

TEST_F(RingBufferCapacityTest, ClearsAndResizesBufferWhenResettingWithNewCapacity)
{
    DestinationArrayType destination(destinationSize);
    destination.insert(source, sourceSize);

    destination.reset(sourceSize * 2);

    EXPECT_EQ(sourceSize * 2, destination.capacity());
    EXPECT_EQ(0U, destination.currentSize());

    destination.insert(source, sourceSize);
    destination.insert(source, sourceSize);

    for (std::size_t i{0}; i < sourceSize; ++i) {
        EXPECT_EQ(source[sourceSize - i % sourceSize - 1], destination.copy(i));
    }
}

TEST_F(RingBufferCapacityTest, KeepsCapacityAndCurrentSizeEqualAfterOversizedBlockIsInserted)
{
    DestinationArrayType destination(destinationSize);

    destination.insert(source, sourceSize);

    EXPECT_NE(sourceSize, destination.currentSize());
    EXPECT_EQ(destination.capacity(), destination.currentSize());
}

} // namespace jb::testing