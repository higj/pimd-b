#pragma once

#include "common.h"
#include "potentials/potential.h"

/* -------------- Isotropic cosine potential -------------- */
class CosinePotential : public Potential {
public:
    CosinePotential(double amplitude, double wavelength, double phase);
    ~CosinePotential() override = default;

    // Potential
    double V(const SingleVec& x) override;

    // Potential gradient
    SingleVec gradV(const SingleVec& x) override;

    // Potential laplacian
    double laplacianV(const SingleVec& x) override;

    [[nodiscard]] bool isFree() const override { return false; }

private:
    double amplitude;   // Amplitude of the potential (V_0)
    double wavelength;  // Wavelength of the potential (L)
    double phase;       // Phase shift of the potential (phi)
    double k;           // Wavenumber of the potential (k = 2*pi/L)
};