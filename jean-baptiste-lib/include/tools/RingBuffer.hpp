#ifndef JB_TOOLS_RINGBUFFER_HPP_
#define JB_TOOLS_RINGBUFFER_HPP_

#include <algorithm>
#include <assert.h>
#include <cstddef>
#include <memory>
#include <new>
#include <type_traits>
#include <utility>

namespace jb::tools {

/**
 * Does not cast into move iterator. Simply return the passed iterator.
 *
 * \tparam Iterator ... The type of iterator. std::is_lvalue_reference<Type&&>::value returns true.
 * \param iterator ... The iterator.
 * \param std::true_type ... Is always std::true_type.
 * \return auto ... Iterator.
 */
template <typename Iterator>
Iterator make_forward_iterator(Iterator iterator, std::true_type)
{
    return iterator;
}

/**
 * Casts iterator into a move iterator.
 *
 * \tparam Iterator ... The type of iterator. std::is_lvalue_reference<Type&&>::value returns false.
 * \param iterator ... The iterator.
 * \param std::false_type ... Is always std::false_type.
 * \return auto ... The moveable version of iterator.
 */
template <typename Iterator>
auto make_forward_iterator(Iterator iterator, std::false_type)
{
    return std::make_move_iterator(iterator);
}

/**
 * Uses SFINAE to cast an iterator into a move iterator when the forwarded type is moveable.
 *
 * \tparam Type ... The type that is to be iterated.
 * \tparam Iterator ... The type of iterator. std::is_lvalue_reference<Type&&>{} is true or false type.
 * \param iterator ... The iterator.
 * \return auto ... Iterator or iterator casted to be moveable.
 */
template <typename Type, typename Iterator>
auto make_forward_iterator(Iterator iterator)
{
    return make_forward_iterator<Iterator>(iterator, std::is_lvalue_reference<Type&&>{});
}

/**
 * Iterates through an array and calls a function on the element together with an incremental index.
 *
 * \tparam Iterator ... The type of iterator.
 * \tparam Function ... The type of function that is to be called on iterated elements.
 * \param first ... A start iterator.
 * \param last ... An end iterator.
 * \param initial ... An initial value that gets incremented inside.
 * \param func ... The function object that is to be called on the element.
 * \return Function ... The function object.
 */
template <typename Iterator, typename Function>
Function enumerate(Iterator first, Iterator last, typename std::iterator_traits<Iterator>::difference_type initial, Function func)
{
    for (; first != last; ++first, ++initial) {
        func(initial, std::forward<typename std::iterator_traits<Iterator>::value_type>(*first));
    }

    return func;
}

/**
 * This defines a ring buffer. When reaching the buffer end while inserting new elements it will overwrite the oldest elements.
 * Uses aligned storage for performance gain.
 * \tparam Element ... The type of elements that are getting stored.
 */
template <typename Element>
class RingBuffer {
public:
    RingBuffer(const RingBuffer<Element>&) = delete;
    RingBuffer& operator=(const RingBuffer<Element>&) = delete;
    RingBuffer(RingBuffer<Element>&&) = delete;
    RingBuffer& operator=(RingBuffer<Element>&&) = delete;

    RingBuffer()
        : mBuffer(new ElementStorageType[0])
        , mWritePosition(0)
        , mCurrentSize(0)
        , mCapacity(0)
    {}

    explicit RingBuffer(const std::size_t capacity)
        : mBuffer(new ElementStorageType[capacity])
        , mWritePosition(0)
        , mCurrentSize(0)
        , mCapacity(capacity)
    {}

    ~RingBuffer(void)
    {
        reset();
    }

    /**
     * Destructs all buffer elements and resets the buffer size.
     */
    void reset()
    {
        if (mCurrentSize == 0) {
            return;
        }

        auto oldestElement = wrapBackward(mWritePosition, mCurrentSize, mCapacity);
        auto supposedNewestElement = oldestElement + mCurrentSize;
        auto factualNewestElement = std::min(supposedNewestElement, mCapacity);

        destruct(oldestElement, factualNewestElement);

        if (supposedNewestElement > mCapacity) {
            destruct(0, supposedNewestElement - mCapacity);
        }

        mWritePosition = 0;
        mCurrentSize = 0;
    }

    /**
     * Destructs all buffer elements and reinitializes the buffer with a new capacity.
     */
    void reset(std::size_t newCapacity)
    {
        reset();

        mCapacity = newCapacity;
        mBuffer.reset(new ElementStorageType[newCapacity]);
    }

    /**
     * Inserts a sample.
     * Declaring and using T brings us an universal reference and we can profit from type deduction.
     * So lvalue and rvalue references can be passed to this member (pefectly forwarded).
     * \param sample ... Universal reference of an object of type Element.
     */
    template <typename T>
    void insert(T&& sample)
    {
        insertImpl(std::forward<T>(sample));
    }

    /**
     * Inserts a block of length samples.
     * \param block ... Universal reference to an array of Element objects.
     * \param length ... The length of the array of Element objects.
     */
    template <typename T>
    void insert(T&& block, std::size_t blockLength)
    {
        insertBlockImpl(std::forward<T>(block), blockLength);
    }

    /**
     * Copies the sample feeded sampleBackwards samples ago.
     * e.g. get(0) Returns the last feeded sample.
     * get(capacity_ - 1) Returns the oldest sample.
     * get(capacity_ * n) returns the last feeded sample.
     * get(capacity_ * n + 1) returns the sample feeded in prior to the last one.
     * \param samplesBackward ... The reverse sample index.
     * \return The sample feeded sampleBackward samples ago.
     */
    Element copy(const std::size_t samplesBackward)
    {
        return extractImpl(samplesBackward);
    }

    /**
     * Copies the last numberOfElements samples from buffer_ into destination.
     * Attention: destination needs to have enough space for numberOfElements samples.
     * \param destination ... Pointer to an array of Element.
     * \param numberOfElements ... The number of elements that are to be copied into destination.
     * \return The number of samples being copied.
     */
    std::size_t copy(Element* destination, std::size_t numberOfElements)
    {
        return extractBlockElements(mBuffer, destination, numberOfElements);
    }

    /**
     * \return The number of elements that the buffer could contain without overwriting older ones.
     */
    std::size_t capacity(void) const
    {
        return mCapacity;
    }

    /**
     * \return The number of elements in the buffer.
     */
    std::size_t currentSize(void) const
    {
        return mCurrentSize;
    }

private:
    // Create a type which optimally aligns Element for internal use.
    struct ElementStorageType {
        alignas(Element) std::byte storage[sizeof(Element)];
    };

    /**
     * Returns a pointer to the storage of an element.
     * \param storage ... Pointer to the element storage.
     * \return Pointer to the element.
     */
    static Element* storagePointer(ElementStorageType* storage)
    {
        return reinterpret_cast<Element*>(storage->storage);
    }

    /**
     * Returns a const pointer to the storage of an element.
     * \param storage ... Pointer to the element storage.
     * \return Const pointer to the element.
     */
    static const Element* storagePointer(const ElementStorageType* storage)
    {
        return reinterpret_cast<const Element*>(storage->storage);
    }

    /**
     * Returns a pointer to the element.
     * \param storage ... Pointer to the element storage.
     * \return Pointer to the element.
     */
    static Element* elementPointer(ElementStorageType* storage)
    {
        return std::launder(storagePointer(storage));
    }

    /**
     * Returns a const pointer to the element.
     * \param storage ... Pointer to the element storage.
     * \return Const pointer to the element.
     */
    static const Element* elementPointer(const ElementStorageType* storage)
    {
        return std::launder(storagePointer(storage));
    }

    /**
     * Increments an index and wraps it around the capacity if necessary.
     * \param index ... The current index.
     * \param capacity ... The capacity of the buffer.
     * \return The incremented and wrapped index.
     */
    static std::size_t incrementWrapped(std::size_t index, const std::size_t capacity)
    {
        assert(capacity > 0);
        assert(index < capacity);

        ++index;
        return index == capacity ? 0 : index;
    }

    /**
     * Wraps an index forward within the capacity.
     * \param index ... The current index.
     * \param capacity ... The capacity of the buffer.
     * \return The wrapped index.
     */
    static std::size_t wrapForward(std::size_t index, const std::size_t capacity)
    {
        assert(capacity > 0);
        assert(index <= capacity);

        return index >= capacity ? index - capacity : index;
    }
    /**
     * Wraps an index backward within the capacity.
     * \param index ... The current index.
     * \param distance ... The distance to move backward.
     * \param capacity ... The capacity of the buffer.
     * \return The wrapped index.
     */
    static std::size_t wrapBackward(const std::size_t index, const std::size_t distance, const std::size_t capacity)
    {
        assert(capacity > 0);
        assert(index < capacity);
        assert(distance <= capacity);

        return distance <= index ? index - distance : capacity - (distance - index);
    }

    /**
     * Wraps a distance within a given range.
     * \param distance ... The distance to wrap.
     * \param range ... The range to wrap within.
     * \return The wrapped distance.
     */
    static std::size_t wrapDistance(const std::size_t distance, const std::size_t range)
    {
        assert(range > 0);

        return distance - (distance / range) * range;
    }

    /**
     * Destructs a range of elements.
     * \param from ... The index of the first element.
     * \param to ... The indexof the last element.
     */
    void destruct(const std::size_t from, const std::size_t to)
    {
        enumerate(
            make_forward_iterator<ElementStorageType>(&mBuffer.get()[from]),
            make_forward_iterator<ElementStorageType>(&mBuffer.get()[to]),
            0,
            [](std::size_t, ElementStorageType&& element) {
                std::destroy_at(elementPointer(&element));
            });
    }

    /**********************************************************************************************
     * insert
     **********************************************************************************************/

    /**
     * Inserts a sample.
     * \param sample ... Universal reference of an object of type Element.
     */
    template <typename T>
    void insertImpl(T&& sample)
    {
        if (mCapacity == 0) {
            return;
        }

        const auto writePosition = mWritePosition;
        const auto bufferIsFull = mCurrentSize == mCapacity;

        if (bufferIsFull) {
            if constexpr (std::is_assignable_v<Element&, T&&>) {
                *elementPointer(&mBuffer.get()[writePosition]) = std::forward<T>(sample);
            }
            else {
                static_assert(std::is_nothrow_move_constructible_v<Element>, "Overwriting a non-assignable element requires noexcept move construction.");
                Element replacement(std::forward<T>(sample));
                std::destroy_at(elementPointer(&mBuffer.get()[writePosition]));
                std::construct_at(storagePointer(&mBuffer.get()[writePosition]), std::move_if_noexcept(replacement));
            }
        }
        else {
            // Assign memory at buffer_[writePosition_] to a pointer of Element.
            // Since we use placement new no memory needs to be allocated.
            // Sample gets moved (when possible).
            std::construct_at(storagePointer(&mBuffer.get()[writePosition]), std::forward<T>(sample));

            // Update the size of the stored array.
            ++mCurrentSize;
        }
        // Update the position.
        mWritePosition = incrementWrapped(mWritePosition, mCapacity);
    }

    /**
     * Moves or copies number elements from source to buffer_.
     * \param source ... Universal reference to an array of Element objects.
     * \param sourceStart ... Start position of the move/copy operation in source.
     * \param number ... The number of elements that is to be moved/copied from source to buffer_.
     * \param destinationStart ... Start position of the move/copy operation in buffer_.
     */
    template <typename T>
    void insertBlockElements(T&& source, const std::size_t sourceStart, const std::size_t destinationStart, const std::size_t number, const std::size_t capacity, std::size_t& currentSize)
    {
        if (number == 0) {
            return;
        }

        auto sourceIterator = make_forward_iterator<T>(&source[sourceStart]);

        if (currentSize == capacity) {
            if constexpr (std::is_assignable_v<Element&, decltype(*sourceIterator)>) {
                for (std::size_t index = 0; index < number; ++index, ++sourceIterator) {
                    *elementPointer(&mBuffer.get()[destinationStart + index]) = *sourceIterator;
                }
            }
            else {
                static_assert(std::is_nothrow_constructible_v<Element, decltype(*sourceIterator)>, "Overwriting non-assignable block elements requires noexcept construction from the source reference.");
                destruct(destinationStart, destinationStart + number);

                std::uninitialized_copy(
                    sourceIterator,
                    sourceIterator + number,
                    storagePointer(&mBuffer.get()[destinationStart]));
            }

            return;
        }

        std::uninitialized_copy(
            sourceIterator,
            // Create an iterator for the last element. Since uninitialized_copy copies like [first, ..., last)
            // we add 1 to point after the array.
            make_forward_iterator<T>(&source[sourceStart + number - 1]) + 1,
            storagePointer(&mBuffer.get()[destinationStart]));

        // Update the size of the stored array.
        if (currentSize < capacity) {
            currentSize = std::min(capacity, currentSize + number);
        }
    }

    /**
     * Inserts a block of length samples.
     * \param block ... Universal reference to an array of Element objects.
     * \param length ... The length of the array of Element objects.
     */
    template <typename T>
    void insertBlockImpl(T&& block, std::size_t blockLength)
    {
        const auto tmpCapacity = mCapacity;
        if (tmpCapacity == 0 || blockLength == 0) {
            return;
        }

        auto tmpCurrentSize = mCurrentSize;
        auto tmpWritePosition = mWritePosition;

        auto adjustedBlockLength = blockLength;
        // Do we need to crop the block to fit into the buffer?
        if (blockLength > tmpCapacity) {
            adjustedBlockLength = tmpCapacity; // Limit the length to the buffer's length.
        }

        // Do we need to divide the block at the physical end of the buffer?
        if ((tmpWritePosition + adjustedBlockLength) > tmpCapacity) {
            const auto firstSegmentSize = tmpCapacity - tmpWritePosition;
            const auto secondSegmentSize = adjustedBlockLength - firstSegmentSize;

            // Forwards samples to the end of the buffer.
            insertBlockElements(
                std::forward<T>(block),
                blockLength - adjustedBlockLength,
                tmpWritePosition,
                firstSegmentSize,
                tmpCapacity,
                tmpCurrentSize);

            tmpWritePosition = wrapForward(tmpWritePosition + firstSegmentSize, tmpCapacity);
            mCurrentSize = tmpCurrentSize;
            mWritePosition = tmpWritePosition;

            // Forwards samples to the begin of the buffer.
            insertBlockElements(
                std::forward<T>(block),
                blockLength - secondSegmentSize,
                0,
                secondSegmentSize,
                tmpCapacity,
                tmpCurrentSize);

            // Update current position.
            tmpWritePosition = secondSegmentSize;
        }
        else {
            // Forwards the entire block to the first free location.
            insertBlockElements(
                std::forward<T>(block),
                blockLength - adjustedBlockLength,
                tmpWritePosition, adjustedBlockLength,
                tmpCapacity,
                tmpCurrentSize);

            // Update current position.
            tmpWritePosition = tmpWritePosition + adjustedBlockLength;
        }

        // Reset the position of the first free element if buffer is full.
        tmpWritePosition = wrapForward(tmpWritePosition, tmpCapacity);

        mCurrentSize = tmpCurrentSize;
        mWritePosition = tmpWritePosition;
    }

    /**********************************************************************************************
     * extract
     **********************************************************************************************/

    /**
     * Copies the block of the samples first ... last.
     * SourceType ... The type of the source array.
     * DestinationType ... The type of the destination array.
     * \param source ... Universal reference to the buffer array.
     * \param destination ... Pointer to the destination array that is to be filled with source elements.
     * \param first ... Index of the first element that is to be copied.
     * \param last ... Index of the last element that is to be copied.
     * \param std::true_type ... Is always std::true_type.
     */
    template <typename SourceType, typename DestinationType>
    void extractBlockElementsImpl(SourceType&& source, DestinationType* destination, std::size_t first, std::size_t last, std::true_type)
    {
        enumerate(
            make_forward_iterator<SourceType>(&source.get()[first]),
            make_forward_iterator<SourceType>(&source.get()[last]),
            0,
            [&destination](std::size_t index, ElementStorageType&& element) {
                destination[index] = *elementPointer(&element);
            });
    }

    /**
     * Moves the block of the samples first ... last.
     * SourceType ... The type of the source array.
     * DestinationType ... The type of the destination array.
     * \param source ... Universal reference to the buffer array.
     * \param destination ... Pointer to the destination array that is to be filled with source elements.
     * \param first ... Index of the first element that is to be moved.
     * \param last ... Index of the last element that is to be moved.
     * \param std::false_type ... Is always std::false_type.
     */
    template <typename SourceType, typename DestinationType>
    void extractBlockElementsImpl(SourceType&& source, DestinationType* destination, std::size_t first, std::size_t last, std::false_type)
    {
        enumerate(
            make_forward_iterator<SourceType>(&source.get()[first]),
            make_forward_iterator<SourceType>(&source.get()[last]),
            0,
            [&destination](std::size_t index, ElementStorageType&& element) {
                destination[index] = std::move(*elementPointer(&element));
            });
    }

    /**
     * Extracts the block of numberOfElements samples.
     * Uses SFINAE (based on the forwarded type of the source array) to get the right method for moving/copying.
     * SourceType ... The type of the source array.
     * DestinationType ... The type of the destination array.
     * \param source ... Universal reference to the buffer array.
     * \param destination ... Pointer to the destination array that is to be filled with source elements.
     * \param numberOfElements ... The number of elements that are to be moved/copied from source to destinatin.
     */
    template <typename SourceType, typename DestinationType>
    std::size_t extractBlockElements(SourceType&& source, DestinationType* destination, std::size_t numberOfElements)
    {
        if (destination == nullptr) {
            return 0;
        }

        const auto tmpCapacity = mCapacity;
        const auto tmpCurrentSize = mCurrentSize;
        const auto tmpWritePosition = mWritePosition;

        // We cannot give back more than we have.
        if (numberOfElements > tmpCurrentSize) {
            numberOfElements = tmpCurrentSize;
        }

        if (numberOfElements > 0) {
            auto start = wrapBackward(tmpWritePosition, numberOfElements, tmpCapacity);
            auto size = ((tmpCapacity - start) < numberOfElements) ? (tmpCapacity - start) : numberOfElements;

            extractBlockElementsImpl(
                std::forward<SourceType>(source),
                destination,
                start,
                start + size,
                std::is_lvalue_reference<SourceType&&>{});

            auto numberOfRemaining = static_cast<int>(numberOfElements - size);
            auto startOfRemaining = wrapForward(start + size, tmpCapacity);

            if (numberOfRemaining > 0) {
                extractBlockElementsImpl(
                    std::forward<SourceType>(source),
                    destination + size,
                    startOfRemaining,
                    startOfRemaining + numberOfRemaining,
                    std::is_lvalue_reference<SourceType&&>{});
            }
        }

        return numberOfElements;
    }

    /**
     * Copies the sample feeded sampleBackwards samples ago.
     * \param samplesBackward ... The reverse sample index.
     * \return The sample feeded sampleBackward samples ago.
     */
    Element extractImpl(const std::size_t samplesBackward)
    {
        if (mCurrentSize == 0) {
            return Element();
        }

        const auto wrappedSamplesBackward = wrapDistance(samplesBackward, mCurrentSize);
        const auto newestIndex = mWritePosition == 0 ? mCapacity - 1 : mWritePosition - 1;
        const auto index = wrapBackward(newestIndex, wrappedSamplesBackward, mCapacity);

        return *elementPointer(&mBuffer.get()[index]);
    }

    std::unique_ptr<ElementStorageType[]> mBuffer;
    std::size_t mWritePosition;
    std::size_t mCurrentSize;
    std::size_t mCapacity;
};

} // namespace jb::tools

#endif // JB_TOOLS_RINGBUFFER_HPP_
