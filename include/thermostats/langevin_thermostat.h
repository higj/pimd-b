#pragma once

#include "thermostats/thermostat.h"

class RandomGenerators;

class LangevinThermostat final : public Thermostat
{
public:
    LangevinThermostat(
        const ThermalContext& thermal_ctx,
        const NormalModesContext& nm_ctx,
        const std::shared_ptr<SystemState>& state,
        const std::shared_ptr<RandomGenerators>& rng,
        double gamma,
        double dt,
        double mass
    );
    ~LangevinThermostat() override = default;

    void momentaUpdate() override;

private:
    std::shared_ptr<RandomGenerators> m_rng;  // Generator for the random noise
    double m_gamma;                           // Damping rate (TODO: Used only in the constructor)
    double m_dt;                              // Time step (TODO: dt needed only in the constructor)
    double m_mass;                            // Mass of the particles (TODO: Used only in the constructor)
    int m_natoms;                             // Number of atoms in the system

    double m_friction_coefficient;
    double m_noise_coefficient;
};
