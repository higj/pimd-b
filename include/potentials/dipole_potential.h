#pragma once

#include "common.h"
#include "potentials/potential.h"

/* -------------- Dipole potential -------------- */
class DipolePotential : public Potential {
public:
    DipolePotential(double strength);
    ~DipolePotential() override = default;

    // Potential
    double V(const SingleVec& x) override;

    // Potential gradient
    SingleVec gradV(const SingleVec& x) override;

    // Potential laplacian
    double laplacianV(const SingleVec& x) override;

    [[nodiscard]] bool isFree() const override { return false; }

private:
    double strength;
};