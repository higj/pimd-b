#pragma once

#include "common.h"
#include "bosonic_exchange.h"

class BosonicExchangeCarrousel final : public BosonicExchange {
public:
    BosonicExchangeCarrousel(const Simulation& _sim);
    ~BosonicExchangeCarrousel() override = default;

protected:
    void assignIndirectionCoords() override;
};
