#pragma once

#include "common.h"
#include "potentials/potential.h"

/* -------------- Isotropic cosine potential -------------- */
class CosinePotential : public Potential {
public:
    CosinePotential(double amplitude, double wavelength, double phase);
    ~CosinePotential() override = default;

    // Potential
    double V(const Vec& x) override;

    // Potential gradient
    Vec gradV(const Vec& x) override;

    // Potential laplacian
    double laplacianV(const Vec& x) override;

    [[nodiscard]] bool isFree() const override { return false; }

private:
    double amplitude;   // Amplitude of the potential (V_0)
    double wavelength;  // Wavelength of the potential (L)
    double phase;       // Phase shift of the potential (phi)
    double k;           // Wavenumber of the potential (k = 2*pi/L)
};