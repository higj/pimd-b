#pragma once

#include "common.h"
#include "potentials/potential.h"

/* -------------- Isotropic cosine potential -------------- */
class CosinePotential : public Potential {
public:
    CosinePotential(int start_potential_activation, int finish_potential_activation, double amplitude, double wavelength, double phase);
    ~CosinePotential() override = default;

    // Potential
    double V(const dVec& x) override;

    // Potential gradient
    dVec gradV(const dVec& x) override;

    // Potential laplacian
    double laplacianV(const dVec& x) override;

private:
    double amplitude;   // Amplitude of the potential (V_0)
    double wavelength;  // Wavelength of the potential (L)
    double phase;       // Phase shift of the potential (phi)
    double k;           // Wavenumber of the potential (k = 2*pi/L)
};