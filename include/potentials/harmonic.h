#pragma once

#include "common.h"
#include "potentials/potential.h"

/* -------------- Isotropic harmonic potential -------------- */
class HarmonicPotential : public Potential {
public:
    HarmonicPotential(int start_potential_activation, int finish_potential_activation, double mass, double omega);
    ~HarmonicPotential() override = default;

    // Potential
    double V(const dVec& x) override;

    // Potential gradient
    dVec gradV(const dVec& x) override;

    // Potential laplacian
    double laplacianV(const dVec& x) override;

private:
    double mass;  // Mass of the particle experiencing the harmonic potential
    double omega; // Angular frequency of the oscillator
    double k;     // Force constant of the oscillator (k=mw^2)
};
