#include "thermostats/thermostat.h"
#include "simulation.h"
#include "thermostats/thermostat_coupling.h"

Thermostat::Thermostat(Simulation& _sim, bool normal_modes) : sim(_sim) {
    // Choose coupling (Cartesian coords or normal modes of distinguishable ring polymers)
    if (normal_modes) {
        coupling = std::make_unique<NMCoupling>(_sim);
    } else {
        coupling = std::make_unique<CartesianCoupling>(_sim);
    }
}

// This is the step function of a general thermostat, called in the simulation's run loop
void Thermostat::step() {
    coupling->mpiCommunication();
    momentaUpdate();
    coupling->updateCoupledMomenta();
}

// This is an update of the momenta within the thermostat step, unique for each thermostat
void Thermostat::momentaUpdate() {}

double Thermostat::getAdditionToH() {
    return 0;
}