#pragma once

#include "common.h"

class OldBosonicExchange {
public:
    OldBosonicExchange(int nbosons, int np, int bead_num, double beta, double spring_constant,
        const dVec x, const dVec x_prev, const dVec x_next, bool pbc, double L);
    ~OldBosonicExchange();

    double get_potential() const;

    void spring_force(dVec &f);

    double prim_estimator();

private:
    void diff_two_beads(const dVec x1, int l1, const dVec x2, int l2, double diff[NDIM]);
    void spring_force_last_bead(dVec& f);
    void spring_force_first_bead(dVec& f);
    void spring_force_interior_bead(dVec& f);

    int neighbor_of_first(int ptcl_idx);
    int neighbor_of_last(int ptcl_idx);

    double get_elongest();

    const int nbosons;
    const int np;
    const int bead_num;

    double spring_constant;
    double beta;
    const dVec x;
    const dVec x_prev;
    const dVec x_next;

    bool pbc;
    double L; // Linear size of the system

    std::vector<int> labels;  // Particle labels
    double e_longest;
};