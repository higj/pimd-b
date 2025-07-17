#pragma once

#include "observables/observable.h"
#include "contexts/bead_context.h"
#include "strategies/observables/bosonic_probability_strategy.h"

#include <memory>

class BosonicExchangeBase;

class BosonicObservable final : public Observable
{
public:
    /**
     * @brief Bosonic observable class constructor.
     */
    BosonicObservable(
        const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange,
        const BeadContext& bead_ctx,
        bool bosonic,
        int out_freq,
        const std::string& out_unit
    );

    /**
     * @brief Calculates quantities pertaining to bosonic exchange.
     */
    void calculate() override;

private:
    std::unique_ptr<BosonicProbabilityStrategy> m_bosonic_prob_strategy;
};
