#include "thermostats/thermostat.h"
#include "simulation.h"
#include "thermostats/thermostat_coupling.h"

Thermostat::Thermostat(Coupling& coupling, Params& param_obj) : coupling(coupling) {
    getVariant(param_obj.sys["natoms"], natoms);
    getVariant(param_obj.sim["nbeads"], nbeads);
    getVariant(param_obj.sys["mass"], mass);
    getVariant(param_obj.sim["dt"], dt);
    
    double temperature;
    getVariant(param_obj.sys["temperature"], temperature);
    beta = 1.0 / (Constants::kB * temperature);

#if IPI_CONVENTION
    // i-Pi convention [J. Chem. Phys. 133, 124104 (2010)]
    thermo_beta = beta / nbeads;
#else
    // Tuckerman convention
    thermo_beta = beta;
#endif
}

// This is the step function of a general thermostat, called in the simulation's run loop
void Thermostat::step() {
    coupling.mpiCommunication();
    momentaUpdate();
    coupling.updateCoupledMomenta();
}

// This is an update of the momenta within the thermostat step, unique for each thermostat
void Thermostat::momentaUpdate() {}

double Thermostat::getAdditionToH() {
    return 0;
}