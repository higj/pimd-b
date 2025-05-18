#pragma once

#include <vector>

#include "thermostats/thermostat.h"
#include "contexts/thermostats/nose_hoover_thermostat_context.h"

class Simulation;
class Coupling;

class NoseHooverThermostat : public Thermostat {
public:
    NoseHooverThermostat(const ThermostatContext& context, const NoseHooverThermostatContext& nh_context);
    ~NoseHooverThermostat() override = default;

    void momentaUpdate() override;
    double getAdditionToH() override; // The equations of motion conserve H + additionToH, where H is the Hamiltonian of the physical system
protected:
    double singleChainStep(const double& current_energy, const int& index);
    double Q1, Qi; // The masses of the eta
    double dt2, dt4, dt8;
    double required_energy; // A measure of the required preserved energy (kT times the number of degrees of freedom associated with the chain)
    std::vector<double> eta_dot, eta_dot_dot; // The first and second derivatives of eta
    std::vector<double> eta; 
    double singleChainGetAdditionToH(const int& ndof, const int& index);
    NoseHooverThermostatContext m_nh_context;
};

/* -------------------------------- */

class NoseHooverNpThermostat final : public NoseHooverThermostat {
public:
    NoseHooverNpThermostat(const ThermostatContext& context, const NoseHooverThermostatContext& nh_context);
    ~NoseHooverNpThermostat() override = default;

    void momentaUpdate() override;
    double getAdditionToH() override; // The equations of motion conserve H + additionToH, where H is the Hamiltonian of the physical system
};

/* -------------------------------- */

class NoseHooverNpDimThermostat final : public NoseHooverThermostat {
public:
    NoseHooverNpDimThermostat(const ThermostatContext& context, const NoseHooverThermostatContext& nh_context);
    ~NoseHooverNpDimThermostat() override = default;

    void momentaUpdate() override;
    double getAdditionToH() override; // The equations of motion conserve H + additionToH, where H is the Hamiltonian of the physical system
};
