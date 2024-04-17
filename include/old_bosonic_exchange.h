#pragma once

#include "common.h"
#include "bosonic_exchange_base.h"

class OldBosonicExchange final : public BosonicExchangeBase {
public:
    OldBosonicExchange(int nbosons_, int np_, int bead_num_, double beta_, double spring_constant_,
        const dVec& x_, const dVec& x_prev_, const dVec& x_next_, bool pbc_, double size_);
    ~OldBosonicExchange() override = default;

    double effectivePotential() override;
    void prepare() override;
    double primEstimator() override;

protected:
    void springForceFirstBead(dVec& f) override;
    void springForceLastBead(dVec& f) override;

private:
    int firstBeadNeighbor(int ptcl_idx) const;
    int lastBeadNeighbor(int ptcl_idx) const;

    double getMaxExteriorSpringEnergy();

    std::vector<int> labels;  // Particle labels

    double e_shift;
};