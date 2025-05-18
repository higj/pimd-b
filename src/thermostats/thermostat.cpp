#include "thermostats/thermostat.h"

#include "core/system_state.h"
#include "thermostats/thermostat_coupling.h"

Thermostat::Thermostat(const ThermostatContext& context) : m_context(context) {
    const auto& momenta_ptr = std::shared_ptr<dVec>(context.state, &context.state->momenta);
    // Choose coupling (Cartesian coords or normal modes of distinguishable ring polymers)
    if (m_context.couple_to_nm) {
        coupling = std::make_unique<NormalModesCoupling>(
            momenta_ptr,
            context.normal_modes,
            context.state->currentBead()
        );
    } else {
        coupling = std::make_unique<CartesianCoupling>(
            momenta_ptr
        );
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