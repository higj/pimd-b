#pragma once

#include "observables/observable.h"

#include <memory>

class BosonicProbabilityStrategy;

/**
 * @brief Observable for the sign - the ratio of the fermionic to bosonic statistical weights.
 *
 * The strategy object is a thin wrapper around the shared BosonicExchangeBase
 * instance held by StatisticsManager.
 */
class BosonicSignObservable final : public Observable {
public:
    BosonicSignObservable(int out_freq, const std::string& out_unit);
    ~BosonicSignObservable() override;

    void calculate() override;

private:
    std::unique_ptr<BosonicProbabilityStrategy> m_strategy;
};