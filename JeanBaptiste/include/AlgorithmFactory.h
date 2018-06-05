#pragma once

#include "Algorithm.h"
#include <boost/hana.hpp>
#include <cassert>
#include <functional>
#include <iostream>
#include <memory>
#include "SubTask.h"
#include <unordered_map>

namespace hana = boost::hana;

namespace jeanbaptiste
{
    /** A factory for FFT algorithms of different stage. Each stage is used for a certain count of data samples.
        \param Begin ... The starting index of supported FFT algorithm stages.
        \param End ... The end index of supported FFT algorithm stages.
        \param Complex ... The complex data type.
    */
    template <std::size_t Begin,
              std::size_t End,
              typename Radix,
              typename Decimation,
              typename Direction,
              typename Window,
              typename Complex>
    class AlgorithmFactory
    {
        /** Create a map of FFT algorithm stages at compile time.
            The key of a map element is the stage of an FFT algorithm - used as its ID.
            The value of a map element is the FFT algorithm which contains a tuple of sub tasks.
            \return hana::map ... A map of stages and FFT algorithms.
        */
        static constexpr auto createMapOfAlgorithms(void)
        {
            auto stages = hana::make_range(hana::int_c<Begin>, hana::int_c<End>);

            return hana::unpack(stages, [](auto... stage)
            {
                return hana::make_map(hana::make_pair(stage, hana::template_<Algorithm>(stage, hana::type_c<Radix>,
                    hana::type_c<Decimation>, hana::type_c<Direction>, hana::type_c<Window>, hana::type_c<Complex>))...);
            });
        }

        static constexpr auto algorithmMap_ = createMapOfAlgorithms();

        using Callback = std::function<std::unique_ptr<ExecutableAlgorithm<Complex>>()>;
        std::unordered_map<std::size_t, Callback> dynamicAlgorithmMap_;

    public:
        /** Creates a runtime map containing concrete FFT algorithm initializations.
        */
        AlgorithmFactory(void)
        {
            // Create a runtime map of executable algorithms.
            hana::for_each(algorithmMap_, [&](const auto constantPair)
            {
                dynamicAlgorithmMap_.insert(
                {
                    decltype(+hana::first(constantPair))::value,
                    [&]() -> std::unique_ptr<ExecutableAlgorithm<Complex>>
                    {
                        using AlgorithmType = typename decltype(+hana::second(constantPair))::type;
                        return std::make_unique<AlgorithmType>();
                    }
                });
            });
        }

        /** Creates a pointer to an FFT algorithm instantiation.
            \param[in] stage ... The stage of the FFT algorithm which is to be returned.
            \return std::unique_ptr ... Pointer to the FFT algorithm instantiation.
        */
        std::unique_ptr<ExecutableAlgorithm<Complex>> getAlgorithm(const std::size_t stage)
        {
            auto callback = dynamicAlgorithmMap_.find(stage);
            assert(callback != dynamicAlgorithmMap_.end() && "Trying to find algorithm of unknown stage.");

            return callback->second();
        }
    };
}