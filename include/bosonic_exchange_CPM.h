#pragma once

#include "common.h"
#include "bosonic_exchange.h"

class BosonicExchangeCPM final : public BosonicExchange {
public:
    BosonicExchangeCPM(const Simulation& _sim, const int physical_exchange_time);
    ~BosonicExchangeCPM() override = default;

protected:
    void evaluateConnectionProbabilities() override;
    const int physical_exchange_time;
};
