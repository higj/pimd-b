#pragma once

#include "common.h"
#include "bosonic_exchange.h"
#include "random_mars.h"
#include <memory>

class BosonicExchangeShuffle : public BosonicExchange {
public:
    BosonicExchangeShuffle(Params& param_obj, const dVec& coord, const dVec& prev_coord, const dVec& next_coord, const int this_bead);
    ~BosonicExchangeShuffle() override = default;

protected:
    void assignIndirectionCoords() override;
private:
    int timer;
    int step;
    std::vector<int> indexes;
    std::unique_ptr<RanMars> mars_gen;
};
