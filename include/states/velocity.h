#pragma once

#include "states/state.h"

class VelocityState : public State {
public:
    VelocityState(Params& param_obj, int _freq, const std::string& _out_unit, dVec& momenta);

    void initialize(int this_bead) override;
    void output(int step) override;
private:
    double mass;
    dVec& momenta;
};