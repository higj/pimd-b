#pragma once

#include "observables/observable.h"
#include "contexts/thermostat_context.h"

/**
 * @brief Observable for the Nose-Hoover thermostat contribution to the conserved energy.
 *
 * The Nose-Hoover chain contributes an additional term to the classical Hamiltonian
 * that is needed to verify energy conservation.
 *
 * This observable must only be instantiated when a Nose-Hoover thermostat is
 * active. The constructor validates this and throws std::invalid_argument if
 * the thermostat type does not contain "nose_hoover".
 */
class NoseHooverEnergyObservable final : public Observable {
public:
    /**
     * @throws std::invalid_argument if thermostat_ctx does not describe a
     *         Nose-Hoover thermostat.
     */
    NoseHooverEnergyObservable(
        const ThermostatContext& thermostat_ctx,
        int                      out_freq,
        const std::string&       out_unit
    );

    void calculate() override;

private:
    ThermostatContext m_thermostat_ctx;
};