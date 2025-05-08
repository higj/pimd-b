#pragma once

#include "states/state.h"

class Simulation; // Forward declaration

class PositionState : public State {
public:
    PositionState(const Simulation& _sim, int _freq, const std::string& _out_unit);

    void initialize() override;
    void output(int step) override;
};