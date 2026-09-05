#pragma once

#include "common.h"
#include "potentials/potential.h"

/* -------------- Isotropic harmonic potential -------------- */
class HarmonicPotential : public Potential {
public:
    HarmonicPotential(double mass, double omega);
    ~HarmonicPotential() override = default;

    // Potential
    double V(const Vec& x) override;

    // Potential gradient
    Vec gradV(const Vec& x) override;

    // Potential laplacian
    double laplacianV(const Vec& x) override;



private:
    double mass;  // Mass of the particle experiencing the harmonic potential
    double omega; // Angular frequency of the oscillator
    double k;     // Force constant of the oscillator (k=mw^2)
};
