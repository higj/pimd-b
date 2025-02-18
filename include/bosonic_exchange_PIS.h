#pragma once

#include "common.h"
#include "bosonic_exchange.h"

class BosonicExchangePIS final : public BosonicExchange {
public:
    BosonicExchangePIS(const Simulation& _sim);
    ~BosonicExchangePIS() override = default;

protected:
    void assignIndirectionCoords() override;
};
