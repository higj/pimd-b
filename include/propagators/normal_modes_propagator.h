#pragma once

#include "propagators/propagator.h"

class Simulation;
class NormalModes;

class NormalModesPropagator : public Propagator {
public:
    NormalModesPropagator(Simulation& _sim, Params& param_obj, dVec& coord, dVec& momenta, dVec& forces);
    ~NormalModesPropagator() override = default;

    void preForceStep() override;
    void postForceStep() override;

private:
    double freq, c, s, m_omega;
    dVec ext_forces, spring_forces;
    void momentaExternalForces();
};
