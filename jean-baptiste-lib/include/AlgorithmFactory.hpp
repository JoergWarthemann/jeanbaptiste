#ifndef JB_ALGORITHM_FACTORY_HPP_
#define JB_ALGORITHM_FACTORY_HPP_

#include <cassert>
#include <functional>
#include <memory>
#include <unordered_map>

#include <boost/hana.hpp>

#include "Algorithm.hpp"
#include "IExecutableAlgorithm.hpp"

namespace hana = boost::hana;

namespace jb {

/**
 * A factory for FFT algorithms of different stage. Each stage is used for a certain count of data samples.
 * \tparam Begin ... The starting index of supported FFT algorithm stages.
 * \tparam End ... The end index of supported FFT algorithm stages.
 * \tparam Radix ... The radix used in the FFT algorithm.
 * \tparam Decimation ... The decimation type used in the FFT algorithm.
 * \tparam Direction ... The direction of the FFT (forward or backward).
 * \tparam Window ... The windowing function applied to the input data.
 * \tparam Normalization ... The normalization method applied to the FFT results.
 * \tparam Complex ... The complex data type.
 * \tparam Data ... The data type used in the FFT algorithm.
 */
template <std::size_t Begin,
    std::size_t End,
    typename Radix,
    typename Decimation,
    typename Direction,
    typename Window,
    typename Normalization,
    typename Complex,
    typename Data = options::Data_Complex>
class AlgorithmFactory {
    /**
     * Create a map of FFT algorithm stages at compile time.
     * The key of a map element is the stage of an FFT algorithm - used as its ID.
     * The value of a map element is the FFT algorithm which contains a tuple of sub tasks.
     * \return hana::map ... A map of stages and FFT algorithms.
     */
    static constexpr auto createMapOfAlgorithms(void)
    {
        auto stages = hana::make_range(hana::int_c<Begin>, hana::int_c<End>);

        return hana::unpack(stages, [](auto... stage) {
            return hana::make_map(hana::make_pair(
                stage,
                hana::type_c<Algorithm<decltype(stage)::value,
                    Radix,
                    Decimation,
                    Direction,
                    Window,
                    Normalization,
                    Complex,
                    Data>>)...);
        });
    }

    static constexpr auto algorithmMap_ = createMapOfAlgorithms();

    using Callback = std::function<std::unique_ptr<IExecutableAlgorithm<Complex>>()>;
    std::unordered_map<std::size_t, Callback> dynamicAlgorithmMap_;

public:
    /**
     * Creates a runtime map containing concrete FFT algorithm initializations.
     */
    AlgorithmFactory(void)
    {
        // Create a runtime map of executable algorithms.
        hana::for_each(algorithmMap_, [&](const auto constantPair) {
            dynamicAlgorithmMap_.insert({decltype(+hana::first(constantPair))::value,
                [&]() -> std::unique_ptr<IExecutableAlgorithm<Complex>> {
                    using AlgorithmType = typename decltype(+hana::second(constantPair))::type;
                    return std::make_unique<AlgorithmType>();
                }});
        });
    }

    /**
     * Creates a pointer to an FFT algorithm instantiation.
     * \param[in] stage ... The stage of the FFT algorithm which is to be returned.
     * \return std::unique_ptr ... Pointer to the FFT algorithm instantiation.
     */
    std::unique_ptr<IExecutableAlgorithm<Complex>> getAlgorithm(const std::size_t stage)
    {
        auto callback = dynamicAlgorithmMap_.find(stage);
        assert(callback != dynamicAlgorithmMap_.end() && "Trying to find algorithm of unknown stage.");

        return callback->second();
    }
};

} // namespace jb

#endif // JB_ALGORITHM_FACTORY_HPP_