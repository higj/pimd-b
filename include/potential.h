#pragma once

#include "common.h"
#include "units.h"

/* -------------- Basic potential class -------------- */
class Potential {
public:
    Potential();
    virtual ~Potential();

    // Potential
    virtual double V(const dVec& x) { return 0.0; }

    // Potential gradient
    virtual dVec gradV(const dVec& x) {
        dVec temp(x.len());
        return temp;
    }

    // Potential laplacian
    virtual double laplacianV(const dVec& x) { return 0.0; }

    // Tail correction
    double tailV;
};

/* -------------- Isotropic harmonic potential -------------- */
class HarmonicPotential : public Potential {
public:
    HarmonicPotential(double mass, double omega);
    ~HarmonicPotential();

    // Potential
    double V(const dVec& x);

    // Potential gradient
    dVec gradV(const dVec& x);

    // Potential laplacian
    double laplacianV(const dVec& x);

private:
    double mass;  // Mass of the particle experiencing the harmonic potential
    double omega; // Angular frequency of the oscillator
    double k;     // Force constant of the oscillator (k=mw^2)
};

/* -------------- Isotropic double well potential -------------- */
class DoubleWellPotential : public Potential {
public:
    DoubleWellPotential(double mass, double strength, double loc);
    ~DoubleWellPotential();

    // Potential
    double V(const dVec& x);

    // Potential gradient
    dVec gradV(const dVec& x);

    // Potential laplacian
    double laplacianV(const dVec& x);

private:
    double mass;      // Mass of the particle experiencing the double-well potential
    double strength;  // Height of the well
    double loc;       // Distance of the wells from the origin
};

/* -------------- Dipole potential -------------- */
class DipolePotential : public Potential {
public:
    DipolePotential(double strength);
    ~DipolePotential();

    // Potential
    double V(const dVec& x);

    // Potential gradient
    dVec gradV(const dVec& x);

    // Potential laplacian
    double laplacianV(const dVec& x);

private:
    double strength;
};