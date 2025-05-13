#pragma once

#include "propagators/propagator.h"

class Simulation;
class NormalModes;

class NormalModesPropagator : public Propagator {
public:
    NormalModesPropagator(Params& param_obj, dVec& coord, dVec& momenta, dVec& forces,
                          dVec& spring_forces, dVec& physical_forces, dVec& prev_coord, dVec& next_coord,
                          int this_bead, NormalModes& normal_modes, bool bosonic);
    ~NormalModesPropagator() override = default;

    virtual void preForceStep();
    virtual void postForceStep();
private:
    double freq, c, s, m_omega;
    dVec& ext_forces, spring_forces, prev_coord, next_coord;
    void momentaExternalForces();
    int this_bead, nbeads;
    double spring_constant;
    NormalModes& normal_modes;
    bool bosonic;
};