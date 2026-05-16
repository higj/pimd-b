#pragma once

#include "common.h"
#include "potentials/potential.h"

/* -------------- Isotropic harmonic potential -------------- */
class HarmonicPotential : public Potential {
public:
    HarmonicPotential(double mass, double omega);
    ~HarmonicPotential() override = default;

    // Potential
    double V(const SingleVec& x) override;

    // Potential gradient
    SingleVec gradV(const SingleVec& x) override;

    // Potential laplacian
    double laplacianV(const SingleVec& x) override;

    [[nodiscard]] bool isFree() const override { return false; }

private:
    double mass;  // Mass of the particle experiencing the harmonic potential
    double omega; // Angular frequency of the oscillator
    double k;     // Force constant of the oscillator (k=mw^2)
};
