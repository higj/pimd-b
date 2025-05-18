#include "observables/classical_observable.h"
#include "thermostats/thermostat.h"
#include "core/exchange_state.h"
#include "bosonic_exchange/bosonic_exchange_base.h"
#include "ring_polymer_utils.h"

ClassicalObservable::ClassicalObservable(const ClassicalObservableContext& obs_context, int out_freq, const std::string& out_unit) :
    Observable(out_freq, out_unit), m_context(obs_context)
{
    m_is_nose_hoover = obs_context.thermostat_type.find("nose_hoover") != std::string::npos;
    //const bool is_nose_hoover = obs_context.thermostat_type == "nose_hoover" || obs_context.thermostat_type == "nose_hoover_np" || obs_context.thermostat_type == "nose_hoover_np_dim";

    initialize({ "temperature", "cl_kinetic", "cl_spring" });

    if (m_is_nose_hoover)
        initializeLabel("nh_energy");
}

void ClassicalObservable::calculate() {
    calculateKineticEnergy();
    calculateSpringEnergy();
    calculateThermostatEnergy();
}

void ClassicalObservable::calculateKineticEnergy() {
    double kinetic_energy = 0.0;

    const auto& momenta = *m_context.momenta;

    for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            kinetic_energy += momenta(ptcl_idx, axis) * momenta(ptcl_idx, axis);
        }
    }

    kinetic_energy *= 0.5 / m_context.mass;
    quantities["cl_kinetic"] = Units::convertToUser("energy", m_out_unit, kinetic_energy);

    // Temperature is calculated according to Tolman's equipartition theorem as the average kinetic 
    // energy per degree of freedom. This might not be accurate for systems with constraints.
    // See [J. Chem. Theory Comput. 2019, 15, 1, 84-94.] for a discussion on the topic.

    /// @todo When zeroing the center of mass motion, the number of degrees of freedom must be reduced by NDIM.
    double dof = NDIM * m_context.natoms * m_context.nbeads;
    double temperature = 2.0 * kinetic_energy / (dof * Constants::kB);

    // In the i-Pi convention, the ring-polymer simulation is performed at a temperature that is P times higher
    // than the actual (quantum) temperature. Therefore, to ensure the quantum temperature is calculated correctly,
    // one must divide the classical temperature by the number of beads.
#if IPI_CONVENTION
    temperature /= m_context.nbeads;
#endif

    /// @todo Allow conversion to different temperature units
    quantities["temperature"] = Units::convertToUser("temperature", "kelvin", temperature);
}

void ClassicalObservable::calculateSpringEnergy() {
    double spring_energy;

    if (m_context.this_bead == 0 && m_context.bosonic) {
        spring_energy = m_context.exchange_state->bosonic_exchange->effectivePotential();
    } else {
        spring_energy = RingPolymerUtils::classicalSpringEnergy(
            *m_context.coord,
            *m_context.prev_coord,
            m_context.spring_constant,
            MINIM, /// TODO: MINIM should become a parameter (mic_spring and mic_potential)
            m_context.box_size
        );
    }

    quantities["cl_spring"] = Units::convertToUser("energy", m_out_unit, spring_energy);
}

void ClassicalObservable::calculateThermostatEnergy() {
    if (m_is_nose_hoover) {
        quantities["nh_energy"] = Units::convertToUser("energy", m_out_unit, m_context.thermostat->getAdditionToH());
    }
}