#pragma once

#include "common.h"
#include "potentials/potential.h"

/* -------------- Dipole potential -------------- */
class DipolePotential : public Potential {
public:
    DipolePotential(double strength);
    ~DipolePotential() override = default;

    // Potential
    double V(const Vec& x) override;

    // Potential gradient
    Vec gradV(const Vec& x) override;

    // Potential laplacian
    double laplacianV(const Vec& x) override;



private:
    double strength;
};