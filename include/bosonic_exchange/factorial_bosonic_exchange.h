#pragma once

#include "common.h"
#include "bosonic_exchange/bosonic_exchange_base.h"

class FactorialBosonicExchange final : public BosonicExchangeBase {
public:
    FactorialBosonicExchange(Params& param_obj, const dVec& coord, const dVec& prev_coord, const dVec& next_coord, const int this_bead);
    ~FactorialBosonicExchange() override = default;

    double effectivePotential() override;
    void prepare() override;
    double primEstimator() override;

    double getDistinctProbability() override;
    double getLongestProbability() override;

    void printBosonicDebug() override;
protected:
    void springForceFirstBead(dVec& f) override;
    void springForceLastBead(dVec& f) override;

private:
    int firstBeadNeighbor(int ptcl_idx) const;
    int lastBeadNeighbor(int ptcl_idx) const;

    double getMinExteriorSpringEnergy();

    std::vector<int> labels;  // Particle labels

    double e_shift;
};
