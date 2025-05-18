#pragma once

#include "common.h"
#include "potentials/potential.h"

/* -------------- Aziz potential -------------- */
class AzizPotential : public Potential {
public:
    AzizPotential();
    ~AzizPotential() override = default;

    // Potential
    double V(const dVec& x) override;

    // Potential gradient
    dVec gradV(const dVec& x) override;

    // Potential laplacian
    double laplacianV(const dVec& x) override;

private:
    double rm, A, epsilon, alpha, D, C6, C8, C10;

    // The auxiliary F-function for the Aziz potential
    double F(const double x) const {
        return (x < D ? exp(-(D / x - 1.0) * (D / x - 1.0)) : 1.0);
    }

    // The derivative of the F-function
    double dF(const double x) const {
        double ix = 1.0 / x;
        double r = 2.0 * D * ix * ix * (D * ix - 1.0) * exp(-(D * ix - 1.0) * (D * ix - 1.0));
        return (x < D ? r : 0.0);
    }

    // The 2nd derivative of the F-function
    double d2F(const double x) const {
        double ix = 1.0 / x;
        double r = 2.0 * D * ix * ix * ix * (2.0 * D * D * D * ix * ix * ix - 4.0 * D * D * ix * ix
            - D * ix + 2.0) * exp(-(D * ix - 1.0) * (D * ix - 1.0));
        return (x < D ? r : 0.0);
    }
};