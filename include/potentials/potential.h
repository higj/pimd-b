#pragma once

#include "common.h"

/* -------------- Basic potential class -------------- */
class Potential {
public:
    Potential();
    virtual ~Potential() = default;

    // Potential
    virtual double V(const dVec& x) {
        return 0.0;
    }

    // Potential gradient
    virtual dVec gradV(const dVec& x) {
        return dVec(x.len()); // Zero vector of the same length as x
    }

    // Potential laplacian
    virtual double laplacianV(const dVec& x) {
        return 0.0;
    }

    // Tail correction
    /// @todo Implement the tail correction
    double tailV;
};