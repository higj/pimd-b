#pragma once

#include "common.h"
#include "potentials/potential.h"

/* -------------- Isotropic anharmonic potential -------------- */
class AnharmonicPotential : public Potential {
public:
    AnharmonicPotential(double mass, double omega, double cubic_const, double quart_const);
    ~AnharmonicPotential() override = default;

    // Potential
    double V(const Vec& x) override;

    // Potential gradient
    Vec gradV(const Vec& x) override;

    // Potential laplacian
    double laplacianV(const Vec& x) override;

private:
    double mass;         // Mass of the particle
    double omega;        // Angular frequency of the harmonic term
    double k;            // Force constant of the harmonic term (k=mw^2)
    double cubic_const;  // Constant for the cubic term
    double quart_const;  // Constant for the quartic term
};
