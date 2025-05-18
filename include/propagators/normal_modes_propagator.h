#pragma once

#include "propagators/propagator.h"

class NormalModes;

class NormalModesPropagator : public Propagator {
public:
    NormalModesPropagator(Params& param_obj, dVec& coord, dVec& momenta, dVec& forces, 
                          dVec& physical_forces, dVec& spring_forces, dVec& prev_coord, dVec& next_coord,
                          int this_bead, NormalModes& normal_modes);
    ~NormalModesPropagator() override = default;

    void preForceStep() override;
    void postForceStep() override;

private:
    double freq, c, s, m_omega;
    dVec& physical_forces;
    dVec& spring_forces;
    dVec& prev_coord;
    dVec& next_coord;

    void momentaExternalForces();
    int this_bead, nbeads;
    double spring_constant;
    NormalModes& normal_modes;
    bool bosonic;
};