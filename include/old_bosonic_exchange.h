#pragma once

#include "common.h"
#include "bosonic_exchange_base.h"

class OldBosonicExchange final : public BosonicExchangeBase {
public:
    OldBosonicExchange(const Simulation& _sim);
    ~OldBosonicExchange() override = default;

    double effectivePotential() override;
    void prepare() override;
    double primEstimator() override;

    double getDistinctProbability() override;
    double getLongestProbability() override;
    void exteriorSpringForce(dVec& f) override;
    void printBosonicDebug() override;

private:
    int firstBeadNeighbor(int ptcl_idx) const;
    int lastBeadNeighbor(int ptcl_idx) const;
    void springForceFirstBead(dVec& f, const dVec& x, const dVec& x_prev, const dVec& x_next);
    void springForceLastBead(dVec& f, const dVec& x, const dVec& x_prev, const dVec& x_next);
    double getMinExteriorSpringEnergy();

    std::vector<int> labels;  // Particle labels

    double e_shift;
};
