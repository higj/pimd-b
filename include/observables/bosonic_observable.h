#pragma once

#include "observables/observable.h"

#include <memory>

class ExchangeState;

class BosonicObservable final : public Observable
{
public:
    /**
     * @brief Bosonic observable class constructor.
     */
    BosonicObservable(
        const std::shared_ptr<ExchangeState>& exchange_state,
        int out_freq,
        const std::string& out_unit
    );

    /**
     * @brief Calculates quantities pertaining to bosonic exchange.
     */
    void calculate() override;

private:
    std::shared_ptr<ExchangeState> m_exchange_state;
};
