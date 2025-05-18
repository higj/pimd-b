#pragma once

#include "states/state.h"

class Params; // Forward declaration

class ForceState : public State {
public:
    ForceState(Params& param_obj, int _freq, const std::string& _out_unit, dVec& _forces);

    void initialize(int this_bead) override;
    void output(int step) override;
private:
    dVec& forces;
};