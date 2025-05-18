#pragma once

#include "observables/observable.h"
#include "contexts/observables/classical_observable_context.h"

class ClassicalObservable final : public Observable {
public:
    /**
     * @brief Classical observable class constructor.
     */
    ClassicalObservable(const ClassicalObservableContext& obs_context, int out_freq, const std::string& out_unit);

    void calculate() override;

private:
    ClassicalObservableContext m_context;
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