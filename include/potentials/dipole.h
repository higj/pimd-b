#pragma once

#include "common.h"
#include "potentials/potential.h"

/* -------------- Dipole potential -------------- */
class DipolePotential : public Potential {
public:
    DipolePotential(int start_potential_activation, int finish_potential_activation, double strength);
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