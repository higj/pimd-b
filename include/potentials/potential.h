#pragma once

#include "common.h"

/* -------------- Basic potential class -------------- */
class Potential {
public:
    Potential(int start_potential_activation, int finish_potential_activation);
    virtual ~Potential() = default;

    // Potential
    virtual double V(const dVec& x) {
        return 0.0;
    }
    double getV(const dVec& x, int md_step);

    // Potential gradient
    virtual dVec gradV(const dVec& x) {
        return dVec(x.len()); // Zero vector of the same length as x
    }
    dVec getGradV(const dVec& x, int md_step);

    // Potential laplacian
    virtual double laplacianV(const dVec& x) {
        return 0.0;
    }

    // Tail correction
    /// @todo Implement the tail correction
    double tailV;
private:
    int start_potential_activation, finish_potential_activation;
};