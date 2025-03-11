#pragma once

#include "common.h"
#include "potentials/potential.h"

/* -------------- Dipole potential -------------- */
class DipolePotential : public Potential {
public:
    DipolePotential(double strength);
    ~DipolePotential() override = default;

    // Potential
    double V(const dVec& x) override;

    // Potential gradient
    dVec gradV(const dVec& x) override;

    // Potential laplacian
    double laplacianV(const dVec& x) override;

private:
    double strength;
};