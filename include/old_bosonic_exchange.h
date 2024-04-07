#pragma once

#include "common.h"
#include "bosonic_exchange_base.h"

class OldBosonicExchange : public BosonicExchangeBase {
public:
    OldBosonicExchange(int nbosons, int np, int bead_num, double beta, double spring_constant,
        const dVec x, const dVec x_prev, const dVec x_next, bool pbc, double L);
    ~OldBosonicExchange();

    double effectivePotential() override;
    void updateCoordinates(const dVec new_x, const dVec new_x_prev, const dVec new_x_next) override;

    double primEstimator();

protected:
    void springForceFirstBead(dVec& f) override;
    void springForceLastBead(dVec& f) override;

private:
    int neighbor_of_first(int ptcl_idx);
    int neighbor_of_last(int ptcl_idx);

    double get_elongest();

    std::vector<int> labels;  // Particle labels

    double e_longest;
};