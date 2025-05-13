#include "observables/classical.h"
#include "thermostats.h"
#include "bosonic_exchange.h"
#include "units.h"
#include <ranges>
#include <iostream>

/**
 * @brief Classical observable class constructor.
 */
ClassicalObservable::ClassicalObservable(Params& param_obj, int _freq, const std::string& _out_unit, int this_bead, 
                                         dVec& prev_coord, dVec& coord, BosonicExchangeBase& bosonic_exchange, Thermostat& thermostat, dVec& momenta) :
    EnergyObservable(param_obj, _freq, _out_unit, this_bead, prev_coord, coord, bosonic_exchange), 
    thermostat(thermostat), 
    momenta(momenta) {
    thermostat_type = std::get<std::string>(param_obj.sim["thermostat_type"]);
    if (thermostat_type == "nose_hoover" || thermostat_type == "nose_hoover_np" || thermostat_type == "nose_hoover_np_dim") {
        initialize({ "temperature", "cl_kinetic", "cl_spring", "nh_energy" });
    }
    else {
        initialize({ "temperature", "cl_kinetic", "cl_spring" });
    }
}

void ClassicalObservable::calculate() {
    calculateKineticEnergy();
    calculateSpringEnergy();
    if (thermostat_type == "nose_hoover" || thermostat_type == "nose_hoover_np" || thermostat_type == "nose_hoover_np_dim") {
        quantities["nh_energy"] = Units::convertToUser("energy", out_unit, thermostat.getAdditionToH());
    }
}

/**
 * @brief Calculates the classical kinetic energy of the ring polymers, as well as the temperature of the system.
 */
void ClassicalObservable::calculateKineticEnergy() {
    double kinetic_energy = 0.0;

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            kinetic_energy += momenta(ptcl_idx, axis) * momenta(ptcl_idx, axis);
        }
    }

    kinetic_energy *= 0.5 / mass;
    quantities["cl_kinetic"] = Units::convertToUser("energy", out_unit, kinetic_energy);

    // Temperature is calculated according to Tolman's equipartition theorem as the average kinetic 
    // energy per degree of freedom. This might not be accurate for systems with constraints.
    // See [J. Chem. Theory Comput. 2019, 15, 1, 84-94.] for a discussion on the topic.

    /// @todo When zeroing the center of mass motion, the number of degrees of freedom must be reduced by NDIM.
    double dof = NDIM * natoms * nbeads;
    double temperature = 2.0 * kinetic_energy / (dof * Constants::kB);

    // In the i-Pi convention, the ring-polymer simulation is performed at a temperature that is P times higher
    // than the actual (quantum) temperature. Therefore, to ensure the quantum temperature is calculated correctly,
    // one must divide the classical temperature by the number of beads.
#if IPI_CONVENTION
    temperature /= nbeads;
#endif

    /// @todo Allow conversion to different temperature units
    quantities["temperature"] = Units::convertToUser("temperature", "kelvin", temperature);
}

/**
 * @brief Calculates the spring energy of the classical ring-polymer system.
 * In the bosonic case, one must use the effective bosonic spring potential
 * for the exterior connection.
 */
void ClassicalObservable::calculateSpringEnergy() {
    double spring_energy;

    if (this_bead == 0 && bosonic) {
        spring_energy = bosonic_exchange.effectivePotential();
    } else {
        spring_energy = classicalSpringEnergy();
    }
    quantities["cl_spring"] = Units::convertToUser("energy", out_unit, spring_energy);
}