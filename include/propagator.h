#pragma once

#include <vector>
#include "mpi.h"
#include "common.h"

class Simulation;
class NormalModes;

class Propagator {
public:
    explicit Propagator(Simulation& _sim);
    virtual ~Propagator() = default;
    
    virtual void step() = 0;
    void momentStep();
    void coordsStep();

protected:
    Simulation& sim; // Reference to the simulation object
};

/* -------------------------------- */

class VelocityVerletPropagator : public Propagator {
public:
    VelocityVerletPropagator(Simulation& _sim);
    ~VelocityVerletPropagator() override = default;

    void step() override;
};

/* -------------------------------- */

class NormalModesPropagator : public Propagator {
public:
    NormalModesPropagator(Simulation& _sim);
    ~NormalModesPropagator() override = default;

    void step() override;

private:
    double freq, c, s, m_omega;
    dVec ext_forces, spring_forces;
    void momentaExternalForces();
};
