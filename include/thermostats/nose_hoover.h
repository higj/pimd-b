#pragma once

#include <vector>
#include <memory>

#include "thermostats/thermostat.h"

class Simulation;
class Coupling;

class NoseHooverThermostat : public Thermostat {
public:
    NoseHooverThermostat(Coupling& coupling, Params& param_obj);
    ~NoseHooverThermostat() override = default;

    void momentaUpdate() override;
    double getAdditionToH() override; // The equations of motion conserve H + additionToH, where H is the Hamiltonian of the physical system
protected:
    double singleChainStep(const double& current_energy, const int& index);
    int nchains; // Number of components in each Nose-Hoover chain
    double Q1, Qi; // The masses of the eta
    double dt2, dt4, dt8;
    double required_energy; // A measure of the required preserved energy (kT times the number of degrees of freedom associated with the chain)
    std::vector<double> eta_dot, eta_dot_dot; // The first and second derivatives of eta
    std::vector<double> eta; 
    double singleChainGetAdditionToH(const int& ndof, const int& index);
};

/* -------------------------------- */

class NoseHooverNpThermostat : public NoseHooverThermostat {
public:
    NoseHooverNpThermostat(Coupling& coupling, Params& param_obj);
    ~NoseHooverNpThermostat() override = default;

    void momentaUpdate() override;
    double getAdditionToH() override; // The equations of motion conserve H + additionToH, where H is the Hamiltonian of the physical system
};

/* -------------------------------- */

class NoseHooverNpDimThermostat : public NoseHooverThermostat {
public:
    NoseHooverNpDimThermostat(Coupling& coupling, Params& param_obj);
    ~NoseHooverNpDimThermostat() override = default;

    void momentaUpdate() override;
    double getAdditionToH() override; // The equations of motion conserve H + additionToH, where H is the Hamiltonian of the physical system
};
