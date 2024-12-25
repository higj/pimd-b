#pragma once

#include "common.h"
#include "bosonic_exchange.h"

class BosonicExchangeCarrousel final : public BosonicExchange {
public:
    BosonicExchangeCarrousel(const Simulation& _sim);
    ~BosonicExchangeCarrousel() override = default;
    void exteriorSpringForce(dVec& f) override;
    void prepare() override;

protected:
    void assignFirstLast(dVec& x_first_bead, dVec& x_last_bead) const override;

private:
    dVec temp_x;
    dVec temp_x_prev;
    dVec temp_x_next;
    void assignTempCoords();
};
