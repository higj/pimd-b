#pragma once

#include "propagators/propagator.h"

class Simulation;
class NormalModes;

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
