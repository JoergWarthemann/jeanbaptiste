#include <array>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "tools/RingBuffer.hpp"

namespace jb::testing {

class RingBufferGetValuesTest
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

class ThrowsOnCopy {
public:
    static inline bool throwOnCopy{false};
    static inline bool throwOnCopyConstruction{false};
    static inline bool throwOnCopyAssignment{false};
    static inline int liveCount{0};

    ThrowsOnCopy()
        : mValue(0)
    {
        ++liveCount;
    }

    explicit ThrowsOnCopy(int value)
        : mValue(value)
    {
        ++liveCount;
    }

    ThrowsOnCopy(const ThrowsOnCopy& other)
        : mValue(other.mValue)
    {
        if (throwOnCopy || throwOnCopyConstruction)
            throw std::runtime_error("copy construction failed");
        ++liveCount;
    }

    ThrowsOnCopy(ThrowsOnCopy&& other) noexcept
        : mValue(other.mValue)
    {
        ++liveCount;
    }

    ThrowsOnCopy& operator=(const ThrowsOnCopy& other)
    {
        if (throwOnCopy || throwOnCopyAssignment)
            throw std::runtime_error("copy assignment failed");
        mValue = other.mValue;
        return *this;
    }

    ThrowsOnCopy& operator=(ThrowsOnCopy&& other) noexcept
    {
        mValue = other.mValue;
        return *this;
    }

    ~ThrowsOnCopy()
    {
        --liveCount;
    }

    int value() const
    {
        return mValue;
    }

private:
    int mValue;
};

TEST_F(RingBufferGetValuesTest, CopiesNewestSingleValuesAfterOversizedInsert)
{
    DestinationArrayType destination(destinationSize);
    destination.insert(source, sourceSize);

    for (std::size_t i{0}; i < destinationSize * 2; ++i) {
        EXPECT_EQ(source[sourceSize - (i % destinationSize) - 1], destination.copy(i));
    }
}

TEST_F(RingBufferGetValuesTest, CopiesNewestBlockValuesAfterOversizedInsert)
{
    DestinationArrayType destination(destinationSize);
    destination.insert(source, sourceSize);

    constexpr std::size_t blockCount{4};
    constexpr std::size_t blockSize{sourceSize / blockCount};

    for (std::size_t i{1}; i <= blockCount; ++i) {
        std::array<ValueType, destinationSize> output{};

        const std::size_t currentBlockSize{i * blockSize};
        const std::size_t copiedCount{destination.copy(output.data(), currentBlockSize)};

        for (std::size_t j{0}; j < copiedCount; ++j) {
            EXPECT_EQ(source[sourceSize - j - 1], output[copiedCount - j - 1]);
        }
    }
}

TEST_F(RingBufferGetValuesTest, CopiesNewestBlockValuesAfterVariousInsertSizes)
{
    DestinationArrayType destination(destinationSize);
    constexpr std::size_t size1{5};
    constexpr std::size_t size2{10};
    constexpr std::size_t size3{15};

    std::vector<ValueType> controlGroup;
    std::copy(source.begin(), source.begin() + size1, std::back_inserter(controlGroup));
    std::copy(source.begin(), source.begin() + size2, std::back_inserter(controlGroup));
    std::copy(source.begin(), source.begin() + size3, std::back_inserter(controlGroup));
    std::copy(source.begin(), source.begin() + size1, std::back_inserter(controlGroup));
    std::copy(source.begin(), source.begin() + size1, std::back_inserter(controlGroup));
    std::copy(source.begin(), source.begin() + size2, std::back_inserter(controlGroup));

    destination.insert(source, size1);
    destination.insert(source, size2);
    destination.insert(source, size3);
    destination.insert(source, size1);
    destination.insert(source, size1);
    destination.insert(source, size2);

    std::array<ValueType, destinationSize> output{};
    EXPECT_EQ(destinationSize, destination.copy(output.data(), destinationSize));

    for (std::size_t i{0}; i < destinationSize; ++i) {
        EXPECT_EQ(controlGroup[controlGroup.size() - destinationSize + i], output[i]);
    }
}

TEST_F(RingBufferGetValuesTest, KeepsNewestValuesWhenOverwritingAfterWrap)
{
    tools::RingBuffer<int> destination(3);

    destination.insert(1);
    destination.insert(2);
    destination.insert(3);
    destination.insert(4);
    destination.insert(5);

    EXPECT_EQ(3U, destination.currentSize());
    EXPECT_EQ(5, destination.copy(0));
    EXPECT_EQ(4, destination.copy(1));
    EXPECT_EQ(3, destination.copy(2));
}

TEST_F(RingBufferGetValuesTest, KeepsExistingElementsAliveWhenOverwriteConstructionFails)
{
    ThrowsOnCopy::throwOnCopy = false;
    ThrowsOnCopy::liveCount = 0;

    {
        tools::RingBuffer<ThrowsOnCopy> destination(2);
        const ThrowsOnCopy one(1);
        const ThrowsOnCopy two(2);
        const ThrowsOnCopy three(3);

        destination.insert(one);
        destination.insert(two);

        ThrowsOnCopy::throwOnCopy = true;
        EXPECT_THROW(destination.insert(three), std::runtime_error);
        ThrowsOnCopy::throwOnCopy = false;

        EXPECT_EQ(2U, destination.currentSize());
        EXPECT_EQ(5, ThrowsOnCopy::liveCount);
        EXPECT_EQ(2, destination.copy(0).value());
        EXPECT_EQ(1, destination.copy(1).value());
    }

    EXPECT_EQ(0, ThrowsOnCopy::liveCount);
}

TEST_F(RingBufferGetValuesTest, KeepsExistingElementsAliveWhenBlockOverwriteAssignmentFails)
{
    ThrowsOnCopy::throwOnCopy = false;
    ThrowsOnCopy::throwOnCopyConstruction = false;
    ThrowsOnCopy::throwOnCopyAssignment = false;
    ThrowsOnCopy::liveCount = 0;

    {
        tools::RingBuffer<ThrowsOnCopy> destination(2);
        const ThrowsOnCopy one(1);
        const ThrowsOnCopy two(2);
        std::array<ThrowsOnCopy, 1> replacement{ThrowsOnCopy(3)};

        destination.insert(one);
        destination.insert(two);

        ThrowsOnCopy::throwOnCopy = true;
        EXPECT_THROW(destination.insert(replacement, replacement.size()), std::runtime_error);
        ThrowsOnCopy::throwOnCopy = false;

        EXPECT_EQ(2U, destination.currentSize());
        EXPECT_EQ(5, ThrowsOnCopy::liveCount);
        EXPECT_EQ(2, destination.copy(0).value());
        EXPECT_EQ(1, destination.copy(1).value());
    }

    EXPECT_EQ(0, ThrowsOnCopy::liveCount);
}

TEST_F(RingBufferGetValuesTest, TracksFirstConstructedSegmentWhenWrappedBlockInsertFails)
{
    ThrowsOnCopy::throwOnCopy = false;
    ThrowsOnCopy::throwOnCopyConstruction = false;
    ThrowsOnCopy::throwOnCopyAssignment = false;
    ThrowsOnCopy::liveCount = 0;

    {
        tools::RingBuffer<ThrowsOnCopy> destination(5);
        std::array<ThrowsOnCopy, 4> initial{ThrowsOnCopy(1), ThrowsOnCopy(2), ThrowsOnCopy(3), ThrowsOnCopy(4)};
        std::array<ThrowsOnCopy, 3> replacement{ThrowsOnCopy(5), ThrowsOnCopy(6), ThrowsOnCopy(7)};

        destination.insert(initial, initial.size());

        ThrowsOnCopy::throwOnCopyAssignment = true;
        EXPECT_THROW(destination.insert(replacement, replacement.size()), std::runtime_error);
        ThrowsOnCopy::throwOnCopyAssignment = false;

        EXPECT_EQ(5U, destination.currentSize());
        EXPECT_EQ(12, ThrowsOnCopy::liveCount);
        EXPECT_EQ(5, destination.copy(0).value());
        EXPECT_EQ(4, destination.copy(1).value());
        EXPECT_EQ(1, destination.copy(4).value());
    }

    EXPECT_EQ(0, ThrowsOnCopy::liveCount);
}

TEST_F(RingBufferGetValuesTest, ReturnsZeroWhenCopyingToNullDestination)
{
    tools::RingBuffer<int> destination(2);

    destination.insert(1);

    EXPECT_EQ(0U, destination.copy(nullptr, 1));
    EXPECT_EQ(1U, destination.currentSize());
    EXPECT_EQ(1, destination.copy(0));
}

TEST_F(RingBufferGetValuesTest, WrapsWithinCurrentSizeBeforeBufferIsFullWhenCopyingSingleValue)
{
    tools::RingBuffer<int> destination(5);

    destination.insert(1);
    destination.insert(2);
    destination.insert(3);

    EXPECT_EQ(3, destination.copy(0));
    EXPECT_EQ(2, destination.copy(1));
    EXPECT_EQ(1, destination.copy(2));
    EXPECT_EQ(3, destination.copy(3));
}

TEST_F(RingBufferGetValuesTest, DoesNotChangeBufferWhenInsertingWithEmptyCapacityAndZeroLength)
{
    tools::RingBuffer<int> empty;
    std::array<int, 3> values{1, 2, 3};

    empty.insert(1);
    empty.insert(values, 0);

    EXPECT_EQ(0U, empty.capacity());
    EXPECT_EQ(0U, empty.currentSize());
    EXPECT_EQ(0, empty.copy(0));
}

} // namespace jb::testing