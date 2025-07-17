#pragma once

#include "observables/observable.h"
#include "contexts/velocity_context.h"
#include "contexts/thermostat_context.h"
#include "contexts/bead_context.h"
#include "contexts/spring_context.h"
#include "contexts/box_context.h"
#include "strategies/observables/classical_spring_energy_strategy.h"

#include <memory>

class BosonicExchangeBase;
class Thermostat;

class ClassicalObservable final : public Observable {
public:
    /**
     * @brief Classical observable class constructor.
     */
    ClassicalObservable(
        const std::shared_ptr<const dVec>& coord,
        const std::shared_ptr<const dVec>& prev_coord,
        const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange,
        bool bosonic,
        const VelocityContext& vel_ctx,
        const ThermostatContext& thermostat_ctx,
        const BeadContext& bead_ctx,
        const SpringContext& spring_ctx,
        const BoxContext& box_ctx,
        int out_freq, 
        const std::string& out_unit
    );

    void calculate() override;

private:
    std::shared_ptr<const dVec> m_coord;
    std::shared_ptr<const dVec> m_prev_coord;
    std::unique_ptr<ClassicalSpringEnergyStrategy> m_spring_energy_strategy;

    VelocityContext m_vel_ctx;
    ThermostatContext m_thermostat_ctx;
    BeadContext m_bead_ctx;
    SpringContext m_spring_ctx;
    BoxContext m_box_ctx;

    bool m_is_nose_hoover = false;  // Is the thermostat a Nose-Hoover?

    /**
     * @brief Calculates the classical kinetic energy of the ring polymers, as well as the temperature of the system.
     */
    void calculateKineticEnergy();

    /**
     * @brief Calculates the spring energy of the classical ring-polymer system.
     * In the bosonic case, one must use the effective bosonic spring potential
     * for the exterior connection.
     */
    void calculateSpringEnergy();

    /**
     * @brief Calculates the energy associated with the thermal fluctuations.
     */
    void calculateThermostatEnergy();
};