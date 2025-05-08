#pragma once

#include "states/state.h"

class Simulation; // Forward declaration

class VelocityState : public State {
public:
    VelocityState(const Simulation& _sim, int _freq, const std::string& _out_unit);

    void initialize() override;
    void output(int step) override;
};