#pragma once

#include "observables/observable.h"

#include <memory>

class BosonicProbabilityStrategy;

/**
 * @brief Observable for the bosonic probability of the longest-cycle (all-exchanged) permutation sector.
 *
 * Reads a pre-computed value from BosonicProbabilityStrategy; the underlying
 * exchange amplitude is evaluated during propagation, so this observable is
 * a lightweight reader with no expensive computation of its own.
 *
 * The strategy object is a thin wrapper around the shared BosonicExchangeBase
 * instance held by StatisticsManager.
 */
class BosonicProbAllObservable final : public Observable {
public:
    BosonicProbAllObservable(int out_freq, const std::string& out_unit);
    ~BosonicProbAllObservable() override;

    void calculate() override;

private:
    std::unique_ptr<BosonicProbabilityStrategy> m_strategy;
};