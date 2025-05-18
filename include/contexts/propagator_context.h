#pragma once

#include <memory>

class SystemState;
class ForceFieldManager;
struct ExchangeState;

struct PropagatorContext {
    // Propagators need access to the system state (including updateNeighboringCoordinates)
    std::shared_ptr<SystemState> state;
    // Propagators need access to the force field manager (including updateSpringForces)
    std::shared_ptr<ForceFieldManager> force_mgr;
    // Propagators need access to bosonic exchange methods
    std::shared_ptr<ExchangeState> exchange_state;
    double dt;
    int natoms;
    int nbeads;
    double mass;
    double omega_p;
    double spring_constant;
    int this_bead;
    bool bosonic;
};