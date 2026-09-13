#include "observables/nose_hoover_energy_observable.h"
#include "thermostats/thermostat.h"
#include "units.h"

#include <stdexcept>

NoseHooverEnergyObservable::NoseHooverEnergyObservable(
    const ThermostatContext& thermostat_ctx,
    int                      out_freq,
    const std::string& out_unit
) : Observable("nh_energy", out_freq, out_unit),
m_thermostat_ctx(thermostat_ctx) {
    if (thermostat_ctx.thermostat_type.find("nose_hoover") == std::string::npos) {
        throw std::invalid_argument(
            "NoseHooverEnergyObservable requires a Nose-Hoover thermostat, got: '"
            + thermostat_ctx.thermostat_type + "'"
        );
    }

    initializeLabel("nh_energy");
}

void NoseHooverEnergyObservable::calculate() {
    quantities["nh_energy"] = Units::convertToUser(
        "energy", m_out_unit, m_thermostat_ctx.thermostat->getAdditionToH()
    );
}