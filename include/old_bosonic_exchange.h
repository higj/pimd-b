#pragma once

#include "common.h"
#include "bosonic_exchange_base.h"

class OldBosonicExchange : public BosonicExchangeBase {
public:
    OldBosonicExchange(int nbosons, int np, int bead_num, double beta, double spring_constant,
        const dVec x, const dVec x_prev, const dVec x_next, bool pbc, double L);
    ~OldBosonicExchange();

    double prim_estimator();

protected:
    void spring_force_first_bead(dVec& f) override;
    void spring_force_last_bead(dVec& f) override;

private:
    int neighbor_of_first(int ptcl_idx);
    int neighbor_of_last(int ptcl_idx);

    double get_elongest();

    std::vector<int> labels;  // Particle labels
    double e_longest;
};