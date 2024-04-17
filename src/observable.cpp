#include "observable.h"
#include "simulation.h"
#include "units.h"

/**
 * @brief Generic observable class constructor
 */
Observable::Observable(const Simulation& _sim, int _freq, const std::string& _out_unit) : 
    sim(_sim), freq(_freq), out_unit(_out_unit) {
}

/**
 * Initializes observables with the given labels.
 * 
 * @param _labels Labels of the quantities to be calculated
 */
void Observable::initialize(const std::vector<std::string>& _labels) {
    for (const std::string& _label : _labels) {
        quantities.insert({ _label, 0.0 });
    }
}

/**
 * @brief Resets the values of the observable to zero.
 * Useful for clearing results from the previous moleclar dynamics step.
 */
void Observable::resetValues() {
    for (auto it = quantities.begin(); it != quantities.end(); ++it) {
        it.value() = 0.0;
    }
}

/**
 * @brief Creates an observable object based on the observable type.
 *
 * @param observable_type Type of observable to create
 * @param _sim Simulation object
 * @param _freq Frequency at which the observable is calculated
 * @param _out_unit Output unit of the observable
 * @return std::unique_ptr<Observable> Pointer to the created observable object
 * @throw std::invalid_argument If the observable type is unknown
 */
std::unique_ptr<Observable> ObservableFactory::createQuantity(const std::string& observable_type, const Simulation& _sim, int _freq, const std::string& _out_unit) {
    if (observable_type == "energy") {
        return std::make_unique<EnergyObservable>(_sim, _freq, _out_unit);
    } else if (observable_type == "classical") {
        return std::make_unique<ClassicalObservable>(_sim, _freq, _out_unit);
    } else {
        throw std::invalid_argument("Unknown observable type.");
    }
}

/**
 * @brief Energy observable class constructor.
 */
EnergyObservable::EnergyObservable(const Simulation& _sim, int _freq, const std::string& _out_unit) : Observable(_sim, _freq, _out_unit) {
    initialize({ "kinetic", "potential", "ext_pot", "int_pot", "virial" });
}

void EnergyObservable::calculate() {
    calculateKinetic();
    calculatePotential();
}

/** 
 * @brief Calculates the contribution of the current imaginary time-slice to
 * the primitive kinetic energy estimator of distinguishable particles.
 */
double EnergyObservable::primitiveKineticDistinguishable() const
{
    double spring_energy = 0.0;

    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double diff = sim.coord(ptcl_idx, axis) - sim.prev_coord(ptcl_idx, axis);
#if MINIM
            if (pbc)
                applyMinimumImage(diff, sim.size);
#endif
            spring_energy += diff * diff;
        }
    }

    spring_energy = 0.5 * sim.mass * sim.omega_p * sim.omega_p * spring_energy;

#if IPI_CONVENTION
    spring_energy /= sim.nbeads;
#endif

    return spring_energy;
}

/**
 * @brief Calculates the quantum kinetic energy of the system using the primitive kinetic energy estimator.
 * Works both for distinguishable particles and bosons.
 */
void EnergyObservable::calculateKinetic() {
    double prefactor = 0.5 * NDIM * sim.natoms / sim.beta;

    if (sim.bosonic) {
#if OLD_BOSONIC_ALGORITHM
        if (sim.this_bead == 0) {
            quantities["kinetic"] = prefactor - sim.bosonic_exchange->primEstimator();
        } else {
            quantities["kinetic"] = prefactor - primitiveKineticDistinguishable();
        }
#else
        if (sim.this_bead == 0) {
            quantities["kinetic"] = prefactor * sim.nbeads + sim.bosonic_exchange->primEstimator();
        }
#endif
    } else {
        quantities["kinetic"] = prefactor - primitiveKineticDistinguishable();
    }

    quantities["kinetic"] = Units::convertToUser("energy", out_unit, quantities["kinetic"]);
}

/**
 * @brief Calculates the quantum potential energy of the system, based on the potential energy estimator.
 * The potential energy is the sum of the external potential energy and the interaction potential energy 
 * across all time-slices, divided by the number of beads. In addition, the method calculates the virial
 * kinetic energy of the system.
 */
void EnergyObservable::calculatePotential() {
    double potential = 0.0; // Total potential energy
    double ext_pot = 0.0;   // Potential energy due to external field
    double int_pot = 0.0;   // Potential energy due to interactions
    double virial = 0.0;    // Virial kinetic energy
    
    dVec physical_forces(sim.natoms);
    physical_forces = (-1.0) * sim.ext_potential->gradV(sim.coord);

    double ext_pot_val = sim.ext_potential->V(sim.coord);
    potential += ext_pot_val;
    ext_pot = ext_pot_val;

    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            virial -= sim.coord(ptcl_idx, axis) * physical_forces(ptcl_idx, axis);
        }
    }

    if (sim.int_pot_cutoff != 0.0) {
        for (int ptcl_one = 0; ptcl_one < sim.natoms; ++ptcl_one) {
            for (int ptcl_two = ptcl_one + 1; ptcl_two < sim.natoms; ++ptcl_two) {
                dVec diff = sim.getSeparation(ptcl_one, ptcl_two);  // Vectorial distance

                if (const double distance = diff.norm(); distance < sim.int_pot_cutoff || sim.int_pot_cutoff < 0.0) {
                    dVec force_on_one;
                    force_on_one = (-1.0) * sim.int_potential->gradV(diff);

                    double int_pot_val = sim.int_potential->V(diff);
                    potential += int_pot_val;
                    int_pot += int_pot_val;

                    for (int axis = 0; axis < NDIM; ++axis) {
                        virial -= sim.coord(ptcl_one, axis) * force_on_one(0, axis);
                    }
                }
            }
        }
    }

    potential /= sim.nbeads;
    int_pot /= sim.nbeads;
    ext_pot /= sim.nbeads;
    virial *= 0.5 / sim.nbeads;

    quantities["potential"] = Units::convertToUser("energy", out_unit, potential);
    quantities["ext_pot"] = Units::convertToUser("energy", out_unit, ext_pot);
    quantities["int_pot"] = Units::convertToUser("energy", out_unit, int_pot);
    quantities["virial"] = Units::convertToUser("energy", out_unit, virial);
}

/**
 * @brief Classical observable class constructor.
 */
ClassicalObservable::ClassicalObservable(const Simulation& _sim, int _freq, const std::string& _out_unit) : Observable(_sim, _freq, _out_unit) {
    initialize({ "temperature", "cl_kinetic", "cl_spring" });
}

void ClassicalObservable::calculate() {
    calculateKineticEnergy();
    calculateSpringEnergy();
}

/**
 * @brief Calculates the classical kinetic energy of the ring polymers, as well as the temperature of the system.
 */
void ClassicalObservable::calculateKineticEnergy() {
    double kinetic_energy = 0.0;

    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double mom = sim.momenta(ptcl_idx, axis);
            kinetic_energy += mom * mom;
        }
    }

    kinetic_energy *= 0.5 / sim.mass;
    quantities["cl_kinetic"] = Units::convertToUser("energy", out_unit, kinetic_energy);

    // Temperature is calculated according to Tolman's equipartition theorem as the average kinetic 
    // energy per degree of freedom. This might not be accurate for systems with constraints.
    // See [J. Chem. Theory Comput. 2019, 15, 1, 84-94.] for a discussion on the topic.

    /// @todo When zeroing the center of mass motion, the number of degrees of freedom must be reduced by NDIM.
    double dof = NDIM * sim.natoms * sim.nbeads;
    double temperature = 2.0 * kinetic_energy / (dof * Constants::kB);

#if IPI_CONVENTION
    temperature /= sim.nbeads;
#endif

    /// @todo Allow conversion to different temperature units
    quantities["temperature"] = Units::convertToUser("temperature", "kelvin", temperature);
}

/**
 * @brief Calculates the spring energy of the system. In the bosonic case, the spring energy is given by the
 * effective bosonic spring potential.
 */
void ClassicalObservable::calculateSpringEnergy() {
    double spring_energy = 0.0;

    if (sim.bosonic) {
#if OLD_BOSONIC_ALGORITHM
        spring_energy = sim.bosonic_exchange->effectivePotential();
#else
        if (sim.this_bead == 0) {
            spring_energy = sim.bosonic_exchange->effectivePotential();
        }
#endif
    } else {
        for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                double diff = sim.coord(ptcl_idx, axis) - sim.prev_coord(ptcl_idx, axis);
                spring_energy += diff * diff;
            }
        }

        spring_energy = 0.5 * sim.mass * sim.omega_p * sim.omega_p * spring_energy;
    }

    quantities["cl_spring"] = Units::convertToUser("energy", out_unit, spring_energy);
}