#pragma once

#include "common.h"
#include "potentials/potential.h"

/* -------------- Isotropic double well potential -------------- */
class DoubleWellPotential : public Potential {
public:
    DoubleWellPotential(double mass, double strength, double loc);
    ~DoubleWellPotential() override = default;

    // Potential
    double V(const SingleVec& x) override;

    // Potential gradient
    SingleVec gradV(const SingleVec& x) override;

    // Potential laplacian
    double laplacianV(const SingleVec& x) override;

    [[nodiscard]] bool isFree() const override { return false; }

private:
    double mass;      // Mass of the particle experiencing the double-well potential
    double strength;  // Height of the well
    double loc;       // Distance of the wells from the origin
};