#include "observable.h"
#include "simulation.h"
#include "units.h"

Observable::Observable(const Simulation& _sim, int _freq, const std::string& _out_unit) : 
    sim(_sim), freq(_freq), out_unit(_out_unit) {
}

void Observable::initialize(std::vector<std::string> _labels) {
    for (const std::string& _label : _labels) {
        quantities.insert({ _label, 0.0 });
    }
}

void Observable::resetValues() {
    for (auto it = quantities.begin(); it != quantities.end(); ++it) {
        it.value() = 0.0;
    }
}

Observable::~Observable() = default;

EnergyObservable::EnergyObservable(const Simulation& _sim, int _freq, const std::string& _out_unit) : Observable(_sim, _freq, _out_unit) {
    initialize({ "classical_kinetic", "classical_potential", "kinetic", "potential", "ext_pot", "int_pot", "virial" });
}

void EnergyObservable::calculate() {
    calculateClassicalKinetic();
    calculateKinetic();
    calculatePotential();
    calculateClassicalPotential();
}

// Calculates the contribution of the current imaginary time-slice to
// the primitive kinetic energy estimator of distinguishable particles.
double EnergyObservable::primitiveKineticDistinguishable() {
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

void EnergyObservable::calculateClassicalKinetic() {
    double prefactor =  1 / (2 * sim.mass), kinetic = 0.0, mom;
    for (int axis = 0; axis < NDIM; ++axis) {
        for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
            mom = sim.momenta(ptcl_idx, axis);
            kinetic += prefactor * mom * mom;
        }
    }
    quantities["classical_kinetic"] = Units::unit_to_user("energy", out_unit, kinetic);
}

void EnergyObservable::calculateClassicalPotential() {
    double prefactor = 0.5 * sim.mass * sim.omega_p *sim.omega_p, potential = 0.0, dx;
    if (sim.bosonic) {
        potential = sim.bosonic_exchange->classical_potential();  // not inside if statement because some calculations might require multithreading
        if (sim.this_bead != 0)  // Trust only the root process calculation
            potential = 0.0;
    } else {
        for (int axis = 0; axis < NDIM; ++axis) {
            for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
                dx = sim.coord(ptcl_idx, axis) - sim.prev_coord(ptcl_idx, axis);
                potential += prefactor * (dx * dx);
            }
        }
    }
    quantities["classical_potential"] = Units::unit_to_user("energy", out_unit, potential) + sim.nbeads * quantities["potential"];  // requires 'quantities["potential"]' to be calculated first
}

void EnergyObservable::calculateKinetic() {
    double prefactor = 0.5 * NDIM * sim.natoms / sim.beta;

    if (sim.bosonic) {
#if OLD_BOSONIC_ALGORITHM
        if (sim.this_bead == 0) {
            quantities["kinetic"] = prefactor - sim.bosonic_exchange->prim_estimator();
        } else {
            quantities["kinetic"] = prefactor - primitiveKineticDistinguishable();
        }
#else
        if (sim.this_bead == 0) {
            quantities["kinetic"] = prefactor * sim.nbeads + sim.bosonic_exchange->prim_estimator();
        }
#endif
    } else {
        quantities["kinetic"] = prefactor - primitiveKineticDistinguishable();
    }

    quantities["kinetic"] = Units::unit_to_user("energy", out_unit, quantities["kinetic"]);
}

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
                double distance = diff.norm(0);                     // Scalar distance

                if (distance < sim.int_pot_cutoff || sim.int_pot_cutoff < 0.0) {
                    dVec force_on_one(1);
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

    quantities["potential"] = Units::unit_to_user("energy", out_unit, potential);
    quantities["ext_pot"] = Units::unit_to_user("energy", out_unit, ext_pot);
    quantities["int_pot"] = Units::unit_to_user("energy", out_unit, int_pot);
    quantities["virial"] = Units::unit_to_user("energy", out_unit, virial);
}

std::unique_ptr<Observable> ObservableFactory::createQuantity(const std::string& observable_type, const Simulation& _sim, int _freq, const std::string& _out_unit) {
    if (observable_type == "energy") {
        return std::make_unique<EnergyObservable>(_sim, _freq, _out_unit);
    }
    else {
        throw std::invalid_argument("Unknown observable type.");
    }
}