#pragma once

#include "common.h"
#include "bosonic_exchange.h"

class BosonicExchangeShuffle final : public BosonicExchange {
public:
    BosonicExchangeShuffle(const Simulation& _sim);
    ~BosonicExchangeShuffle() override = default;

protected:
    void assignIndirectionCoords() override;
};
