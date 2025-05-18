#pragma once

#include <memory>

class SystemState;
class NormalModes;

struct ThermostatContext {
    std::shared_ptr<SystemState> state;
    std::shared_ptr<NormalModes> normal_modes;
    bool couple_to_nm;
    double beta;
    double thermo_beta;
    int natoms;
    int nbeads;
    double dt;
    double mass;
};