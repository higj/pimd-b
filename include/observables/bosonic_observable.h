#pragma once

#include "observables/observable.h"

#include <memory>

class BosonicProbabilityStrategy;

class BosonicObservable final : public Observable
{
public:
    /**
     * @brief Bosonic observable class constructor.
     */
    BosonicObservable(int out_freq, const std::string& out_unit);

    ~BosonicObservable() override;

    /**
     * @brief Calculates quantities pertaining to bosonic exchange.
     */
    void calculate() override;

private:
    std::unique_ptr<BosonicProbabilityStrategy> m_bosonic_prob_strategy;
};
