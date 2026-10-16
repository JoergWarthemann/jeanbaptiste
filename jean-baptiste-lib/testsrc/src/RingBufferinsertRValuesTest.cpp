#include <array>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "tools/RingBuffer.hpp"

namespace jb::testing {

class RingBufferInsertRValuesTest
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

TEST_F(RingBufferInsertRValuesTest, MovesSingleElementsFromSource)
{
    DestinationArrayType destination(destinationSize);

    for (std::size_t i{0}; i < sourceSize; ++i) {
        destination.insert(std::move(source[i]));

        EXPECT_TRUE(source[i].empty());
        EXPECT_EQ(std::min(destination.capacity(), i + 1), destination.currentSize());
    }
}

TEST_F(RingBufferInsertRValuesTest, DoesNotExceedCapacityWhenInsertingOversizedBlock)
{
    DestinationArrayType destination(destinationSize);

    destination.insert(std::move(source), sourceSize);
    EXPECT_EQ(destinationSize, destination.currentSize());
}

TEST_F(RingBufferInsertRValuesTest, MovesNewestElementsOnlyWhenInsertingOversizedBlock)
{
    DestinationArrayType destination(destinationSize);
    SourceArrayType sourceControl = source;

    destination.insert(std::move(source), sourceSize);

    for (std::size_t i{0}; i < destinationSize; ++i) {
        EXPECT_EQ(sourceControl[sourceSize - i - 1], destination.copy(i));
    }

    for (std::size_t i{0}; i < sourceSize - destinationSize; ++i) {
        EXPECT_FALSE(source[i].empty());
    }
}

TEST_F(RingBufferInsertRValuesTest, MovesNewestElementsWhenInsertingDifferentSizedBlocks)
{
    DestinationArrayType destination(destinationSize);
    SourceArrayType sourceControl = source;
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

    destination.insert(std::move(source1), size1);
    destination.insert(std::move(source2), size2);
    destination.insert(std::move(source3), size3);

    auto checkMovedPrefix = [&sourceControl](std::size_t index, const SourceArrayType& sourceArray, std::size_t movedSize) {
        if (index < movedSize)
            EXPECT_TRUE(sourceArray[index].empty());
        else
            EXPECT_EQ(sourceControl[index], sourceArray[index]);
    };

    for (std::size_t i{0}; i < sourceSize; ++i) {
        checkMovedPrefix(i, source1, size1);
        checkMovedPrefix(i, source2, size2);
        checkMovedPrefix(i, source3, size3);
    }

    auto currentlyExpected = controlGroup.end();
    for (std::size_t i{0}; i < destinationSize; ++i) {
        --currentlyExpected;
        EXPECT_EQ(*currentlyExpected, destination.copy(i));
    }
}

} // namespace jb::testing