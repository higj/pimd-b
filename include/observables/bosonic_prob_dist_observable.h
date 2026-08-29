#pragma once

#include "observables/observable.h"

#include <memory>

class BosonicProbabilityStrategy;

/**
 * @brief Observable for the bosonic probability of the distinct (non-exchanged) permutation sector.
 *
 * Reads a pre-computed value from BosonicProbabilityStrategy; the underlying
 * exchange amplitude is evaluated during propagation, so this observable is
 * a lightweight reader with no expensive computation of its own.
 *
 * The strategy object is a thin wrapper around the shared BosonicExchangeBase
 * instance held by StatisticsManager.
 */
class BosonicProbDistObservable final : public Observable {
public:
    BosonicProbDistObservable(int out_freq, const std::string& out_unit);
    ~BosonicProbDistObservable() override;

    void calculate() override;

private:
    std::unique_ptr<BosonicProbabilityStrategy> m_strategy;
};