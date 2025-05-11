#pragma once

#include "states/state.h"

class PositionState : public State {
public:
    PositionState(Params& param_obj, int _freq, const std::string& _out_unit, dVec& coords);

    void initialize(int this_bead) override;
    void output(int step) override;
private:
    dVec& coord;
};