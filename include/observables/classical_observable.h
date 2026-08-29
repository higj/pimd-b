#pragma once

#include "observables/observable.h"
#include "contexts/velocity_context.h"
#include "contexts/thermostat_context.h"
#include "contexts/bead_context.h"
#include "contexts/spring_context.h"
#include "contexts/box_context.h"

#include <memory>

class BosonicExchangeBase;
class Thermostat;
class ClassicalSpringEnergyStrategy;

class ClassicalObservable final : public Observable {
public:
    /**
     * @brief Classical observable class constructor.
     */
    ClassicalObservable(
        const std::shared_ptr<const VecArray>& coord,
        const std::shared_ptr<const VecArray>& prev_coord,
        const ThermostatContext& thermostat_ctx,
        const SpringContext& spring_ctx,
        const BoxContext& box_ctx,
        int out_freq, 
        const std::string& out_unit
    );

    ~ClassicalObservable() override;

    void calculate() override;

private:
    std::shared_ptr<const VecArray> m_coord;
    std::shared_ptr<const VecArray> m_prev_coord;
    std::unique_ptr<ClassicalSpringEnergyStrategy> m_spring_energy_strategy;

    ThermostatContext m_thermostat_ctx;
    SpringContext m_spring_ctx;
    BoxContext m_box_ctx;

    bool m_is_nose_hoover = false;  // Is the thermostat a Nose-Hoover?

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