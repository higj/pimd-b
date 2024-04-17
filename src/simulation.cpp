#include <ranges>
#include <fstream>
#include <filesystem>

#include "observable.h"
#include "simulation.h"

#if OLD_BOSONIC_ALGORITHM
#include "old_bosonic_exchange.h"
#else
#include "bosonic_exchange.h"
#endif

Simulation::Simulation(int& rank, int& nproc, Params& param_obj, unsigned int seed) :
    bosonic_exchange(nullptr),
    rand_gen(seed + rank),
    this_bead(rank),
    nproc(nproc)
{
    getVariant(param_obj.sim["nbeads"], nbeads);
    getVariant(param_obj.sim["dt"], dt);
    getVariant(param_obj.sim["sfreq"], sfreq);
    getVariant(param_obj.sim["steps"], steps);
    getVariant(param_obj.sim["enable_t"], enable_t);

    getVariant(param_obj.sim["threshold"], threshold);
    threshold *= steps;  // Threshold in terms of steps

    getVariant(param_obj.sim["gamma"], gamma);
    getVariant(param_obj.sim["bosonic"], bosonic);
    getVariant(param_obj.sim["fixcom"], fixcom);
    getVariant(param_obj.sim["pbc"], pbc);

    getVariant(param_obj.out["positions"], out_pos);
    getVariant(param_obj.out["velocities"], out_vel);
    getVariant(param_obj.out["forces"], out_force);

    getVariant(param_obj.sys["temperature"], temperature);
    getVariant(param_obj.sys["natoms"], natoms);
    getVariant(param_obj.sys["size"], size);

    beta = 1.0 / (Constants::kB * temperature);

    getVariant(param_obj.sys["mass"], mass);

#if IPI_CONVENTION
    // i-Pi convention [J. Chem. Phys. 133, 124104 (2010)]
    omega_p = nbeads / (beta * Constants::hbar);
#else
    // Tuckerman convention
    omega_p = sqrt(nbeads) / (beta * Constants::hbar);
#endif

    spring_constant = mass * omega_p * omega_p;

    // Get the seed from the config file
    getVariant(param_obj.sim["seed"], params_seed);

    // Based on the provided seed, make sure each process has a unique seed
    rand_gen.seed(params_seed + rank);
    mars_gen = std::make_unique<RanMars>(params_seed + rank);

    init_pos_type = std::get<std::string>(param_obj.sim["init_pos_type"]);
    init_vel_type = std::get<std::string>(param_obj.sim["init_vel_type"]);

    // Initialize the coordinate, momenta, and force arrays
    coord = dVec(natoms);
    prev_coord = dVec(natoms);
    next_coord = dVec(natoms);
    momenta = dVec(natoms);
    forces = dVec(natoms);

    initializePositions(coord, param_obj.sim);
    initializeMomenta(momenta);

    // Initialize the potential based on the input
    external_potential_name = std::get<std::string>(param_obj.external_pot["name"]);
    interaction_potential_name = std::get<std::string>(param_obj.interaction_pot["name"]);
    int_pot_cutoff = std::get<double>(param_obj.interaction_pot["cutoff"]);

    // For cubic cells with PBC, the cutoff distance must be no greater than L/2 for consistency with
    // the minimum image convention (see 1.6.3 in Allen & Tildesley).
    if (pbc)
        int_pot_cutoff = std::min(int_pot_cutoff, 0.5 * size);

    ext_potential = initializePotential(external_potential_name, param_obj.external_pot);
    int_potential = initializePotential(interaction_potential_name, param_obj.interaction_pot);

    // Observables
    /// @todo Add the ability to specify the observables in the input file
    ObservableFactory obs_factory;
    observables.push_back(obs_factory.createQuantity("energy", *this, sfreq, "kelvin"));
    observables.push_back(obs_factory.createQuantity("classical", *this, sfreq, "kelvin"));

    // Update the coordinate arrays of neighboring particles
    updateNeighboringCoordinates();

    if (bosonic) {
#if OLD_BOSONIC_ALGORITHM
        bosonic_exchange = std::make_unique<OldBosonicExchange>(natoms, nbeads, this_bead, beta, spring_constant, coord, prev_coord, next_coord, pbc, size);
#else
        bosonic_exchange = std::make_unique<BosonicExchange>(natoms, nbeads, this_bead, beta, spring_constant, coord, prev_coord, next_coord, pbc, size);
#endif
    }
}

Simulation::~Simulation() = default;

/**
 * Generates random positions in the interval [-L/2, L/2] for each particle along each axis.
 * 
 * @param pos_arr Array to store the generated positions.
 */
void Simulation::genRandomPositions(dVec &pos_arr) {
    /// @todo Add ability to generate non-random positions (e.g., lattice)
    std::uniform_real_distribution<double> u_dist(-0.5 * size, 0.5 * size);

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            pos_arr(ptcl_idx, axis) = u_dist(rand_gen);
        }
    }
}

/**
 * Generates momenta according to the Maxwell-Boltzmann distribution.
 * 
 * @param momenta_arr Array to store the generated momenta.
 */
void Simulation::genMomentum(dVec &momenta_arr) {
    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            momenta_arr(ptcl_idx, axis) = mass * sampleMaxwellBoltzmann();
        }
    }
}

/**
 * Samples velocities from the Maxwell-Boltzmann distribution.
 * 
 * @return Sampled velocity from the Maxwell-Boltzmann distribution.
 */
double Simulation::sampleMaxwellBoltzmann() {
#if IPI_CONVENTION
    std::normal_distribution<double> normal(0.0, 1/sqrt(beta * mass / nbeads));
#else
    std::normal_distribution<double> normal(0.0, 1 / sqrt(beta * mass));
#endif

    return normal(rand_gen);    
}

/**
 * @brief Perform a single Langevin step for the molecular dynamics simulation.
 */
void Simulation::langevinStep() {
    //std::normal_distribution<double> normal; // At the moment we use the Marsaglia generator

    double a = exp(-0.5 * gamma * dt);

#if IPI_CONVENTION
    double b = sqrt((1 - a * a) * mass * nbeads / beta);
#else
    double b = sqrt((1 - a * a) * mass / beta);
#endif

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double noise = mars_gen->gaussian(); // LAMMPS pimd/langevin random numbers

            // Perturb the momenta with a Langevin thermostat
            momenta(ptcl_idx, axis) = a * momenta(ptcl_idx, axis) + b * noise;
        }
    }
}

/**
 * @brief Velocity-Verlet step for the molecular dynamics simulation.
 */
void Simulation::velocityVerletStep() {
    // First step: momenta are propagated by half a step ("B" step)
    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            momenta(ptcl_idx, axis) += 0.5 * dt * forces(ptcl_idx, axis);
        }
    }

    // Second step: positions are propagated using the new momenta ("A" step)
    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            coord(ptcl_idx, axis) += dt * momenta(ptcl_idx, axis) / mass;
        }
    }

#if WRAP
    if (pbc) {
        for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                periodicWrap(coord(ptcl_idx, axis), size);
            }
        }
    }
#endif

#if RECENTER
    // Recentering should be attempted only in the case of periodic boundary conditions.
    // Also, with the current implementation, it can only work for distinguishable particles.
    if (pbc && !bosonic) {
        // If the initial bead has moved outside the fundamental cell,
        // then rigidly translate the whole polymer such that
        // the polymer will start in the fundamental cell.

        // Vector that stores by how much the polymers should be translated (same as Delta(w)*L)
        dVec shift(natoms);

        if (this_bead == 0) {
            // In the distinguishable case, the outer loop is the loop over all the ring polymers.
            for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
                for (int axis = 0; axis < NDIM; ++axis) {
                    shift(ptcl_idx, axis) = std::nearbyint(coord(ptcl_idx, axis) / size);
                    // Shift the initial bead back to the fundamental cell
                    coord(ptcl_idx, axis) -= size * shift(ptcl_idx, axis);
                }
            }
        }

        // Broadcast the shift vector from process 0 (first bead) to all other processes
        const int shift_vec_size = shift.size();
        MPI_Bcast(shift.data(), shift_vec_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // For other time-slices, we receive info from the first bead and make a decision whether to move the polymer
        if (this_bead != 0) {
            for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
                for (int axis = 0; axis < NDIM; ++axis) {
                    // Shift the current bead according to the shift of the initial bead
                    coord(ptcl_idx, axis) -= size * shift(ptcl_idx, axis);
                }
            }
        }
    }
#endif

    // Remember to update the neighboring coordinates after every coordinate propagation
    updateNeighboringCoordinates();

    // Third step: forces are updated using the new positions
    updateForces();

    // Fourth step: momenta are propagated once more ("B" step)
    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            momenta(ptcl_idx, axis) += 0.5 * dt * forces(ptcl_idx, axis);
        }
    }
}

/**
 * @brief Perform a molecular dynamics run using the OBABO scheme.
 */
void Simulation::run() {
    std::filesystem::create_directory(Output::FOLDER_NAME);
    std::ofstream out_file;

    if (this_bead == 0) {
        out_file.open(std::format("{}/{}", Output::FOLDER_NAME, Output::MAIN_FILENAME), std::ios::out | std::ios::app);

        // Header of the output file
        out_file << std::format("{:^16s}", "step");

        for (const auto& observable : observables) {
            for (const auto& key : observable->quantities | std::views::keys) {
                out_file << std::vformat(" {:^16s}", std::make_format_args(key));
            }
        }

        out_file << "\n";
    }
    
    // Main loop performing molecular dynamics steps
    for (int step = 0; step <= steps; ++step) {
        for (const auto& observable : observables) {
            observable->resetValues();
        }

        if (step % sfreq == 0) {
            if (out_pos)
                outputTrajectories(step);
            if (out_vel)
                outputVelocities(step);
            if (out_force)
                outputForces(step);
        }

        // "O" step
        if (enable_t) 
            langevinStep();

        // If fixcom=true, the center of mass of the ring polymers is fixed during the simulation
        if (fixcom)
            zeroMomentum();

        velocityVerletStep();

        // "O" step
        if (enable_t)
            langevinStep();

        // Zero momentum after every Langevin step (if needed)
        if (fixcom)
            zeroMomentum();

#if PROGRESS
            printProgress(step, steps, this_bead);
#endif

        // If we have not reached the thermalization threshold, skip to the next step (thermalization stage)
        if (step < threshold)
            continue;

        // Calculate and print the observables (production stage)
        for (const auto& observable : observables) {
            observable->calculate();
        }

        if (step % sfreq == 0) {
            if (this_bead == 0)
                out_file << std::format("{:^16.8e}", static_cast<double>(step));

            for (const auto& observable : observables) {
                // The inner loop is necessary because some observable classes can calculate
                // more than one observable (e.g., "energy" calculates both the kinetic and potential energies).
                for (const double& val : observable->quantities | std::views::values) {
                    double quantity_value = 0.0;
                    double local_quantity_value = val;

                    // Sum the results from all processes (beads)
                    MPI_Allreduce(&local_quantity_value, &quantity_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                    if (this_bead == 0)
                        out_file << std::format(" {:^16.8e}", quantity_value);
                }
            }
            
            if (this_bead == 0)
                out_file << "\n";
        }
    }

    if (this_bead == 0) {
        out_file.close();

        std::ofstream report_file;
        report_file.open(std::format("{}/report.txt", Output::FOLDER_NAME), std::ios::out | std::ios::app);
        printReport(report_file);
        report_file.close();
    }
}

/**
 * Obtains the coordinates from the previous timeslice.
 * 
 * @param prev Vector to store the previous coordinates.
 */
void Simulation::getPrevCoords(dVec &prev) {
    const int coord_size = coord.size();

    MPI_Sendrecv(
        coord.data(),
        coord_size,
        MPI_DOUBLE,
        (this_bead + 1) % nproc,
        0,
        prev.data(),
        coord_size,
        MPI_DOUBLE,
        (this_bead - 1 + nproc) % nproc,
        0,
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE
    );

    // Ensure all processes have completed the neighbor communication
    MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * Obtains the coordinates from the next timeslice.
 * 
 * @param next Vector to store the next coordinates.
 */
void Simulation::getNextCoords(dVec& next) {
    const int coord_size = coord.size();

    MPI_Sendrecv(
        coord.data(),
        coord_size,
        MPI_DOUBLE,
        (this_bead - 1 + nproc) % nproc,
        0,
        next.data(),
        coord_size,
        MPI_DOUBLE,
        (this_bead + 1) % nproc,
        0,
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE
    );

    // Ensure all processes have completed the neighbor communication
    MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * @brief Updates the forces acting on the particles. Evaluates both the spring forces
 * and the external forces.
 */
void Simulation::updateForces() {
    // We distinguish between two types of forces: spring forces between the beads (due to the 
    // classical isomorphism) and the non-spring forces (due to either an external potential 
    // or interactions).
    dVec spring_forces(natoms);
    dVec physical_forces(natoms);

    updateSpringForces(spring_forces);
    updatePhysicalForces(physical_forces);

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
#if IPI_CONVENTION
            // i-Pi convention; exp[-(beta/P)*H_cl]
            forces(ptcl_idx, axis) = spring_forces(ptcl_idx, axis) + physical_forces(ptcl_idx, axis);
#else
            // Corresponds to Eqn. (12.6.4) from Tuckerman; exp[-beta*H_cl]
            forces(ptcl_idx, axis) = spring_forces(ptcl_idx, axis) + physical_forces(ptcl_idx, axis) / nbeads;
#endif
        }
    }
}

/**
 * @brief Updates the neighboring coordinates arrays.
 */
void Simulation::updateNeighboringCoordinates() {
    getPrevCoords(prev_coord);
    getNextCoords(next_coord);
}

/**
 * Updates the spring force vector exerted on a specific bead by the two neighboring beads.
 * In the distinguishable case, the force is given by Eqn. (12.6.4) in Tuckerman (1st ed).
 * In the bosonic case, by default, the forces are evaluated using the algorithm 
 * described in https://doi.org/10.1063/5.0173749. It is also possible to perform the
 * bosonic simulation using the original (ineffcient) algorithm, that takes into
 * account all the N! permutations, by setting OLD_BOSONIC_ALGORITHM to true.
 * 
 * @param spring_force_arr Vector to store the spring forces.
 */
void Simulation::updateSpringForces(dVec& spring_force_arr) const {
    if (bosonic) {
        bosonic_exchange->prepare();
        bosonic_exchange->springForce(spring_force_arr);
    } else {
        for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                double d_prev = prev_coord(ptcl_idx, axis) - coord(ptcl_idx, axis);
                double d_next = next_coord(ptcl_idx, axis) - coord(ptcl_idx, axis);

#if MINIM
                if (pbc) {
                    applyMinimumImage(d_prev, size);
                    applyMinimumImage(d_next, size);
                }
#endif

                spring_force_arr(ptcl_idx, axis) = spring_constant * (d_prev + d_next);
            }
        }
    }
}

/**
 * Updates the physical forces acting on the particles. This includes both the forces
 * due external potentials and the interaction forces between the particles.
 * 
 * @param physical_force_arr Vector to store the physical forces.
 */
void Simulation::updatePhysicalForces(dVec& physical_force_arr) const {
    // Calculate the external forces acting on the particles
    physical_force_arr = (-1.0) * ext_potential->gradV(coord);

    if (int_pot_cutoff != 0.0) {
        for (int ptcl_one = 0; ptcl_one < natoms; ++ptcl_one) {
            for (int ptcl_two = ptcl_one + 1; ptcl_two < natoms; ++ptcl_two) {
                dVec diff = getSeparation(ptcl_one, ptcl_two);  // Vectorial distance

                // If the distance between the particles exceeds the cutoff length
                // then we assume the interaction is negligible and do not bother
                // calculating the force.
                // We use the convention that when cutoff < 0 then the interaction is
                // calculated for all distances.
                if (const double distance = diff.norm(); distance < int_pot_cutoff || int_pot_cutoff < 0.0) {
                    dVec force_on_one;
                    force_on_one = (-1.0) * int_potential->gradV(diff);

                    for (int axis = 0; axis < NDIM; ++axis) {
                        physical_force_arr(ptcl_one, axis) += force_on_one(0, axis);
                        physical_force_arr(ptcl_two, axis) -= force_on_one(0, axis);
                    }
                }
            }
        }
    }
}

/**
 * Returns the vectorial distance between two particles at the same imaginary timeslice.
 * Mathematically equivalent to r1-r2, where r1 and r2 are the position-vectors
 * of the first and second particles, respectively.
 * 
 * @param first_ptcl Index of the first particle.
 * @param second_ptcl Index of the second particle.
 * @return Vectorial distance between the two particles.
 */
dVec Simulation::getSeparation(int first_ptcl, int second_ptcl) const {
    dVec sep;

    for (int axis = 0; axis < NDIM; ++axis) {
        double diff = coord(first_ptcl, axis) - coord(second_ptcl, axis);
#if MINIM
        if (pbc)
            applyMinimumImage(diff, size);
#endif
        sep(0, axis) = diff;
    }
    
    return sep;
}

/**
 * @brief Prints a summary of the simulation parameters at the end of the simulation.
 */
void Simulation::printReport(std::ofstream& out_file) const {
    out_file << "---------\nParameters\n---------\n";

    if (bosonic) {
        out_file << formattedReportLine("Statistics", "Bosonic");
        std::string bosonic_alg_name = "Feldman-Hirshberg [O(N^2) scaling]";

#if OLD_BOSONIC_ALGORITHM
        bosonic_alg_name = "Primitive [O(N!) scaling]";
#endif

        out_file << formattedReportLine("Bosonic algorithm", bosonic_alg_name);
    }
    else {
        out_file << formattedReportLine("Statistics", "Boltzmannonic");
    }

    out_file << formattedReportLine("PBC", pbc);
    out_file << formattedReportLine("Dimension", NDIM);
    out_file << formattedReportLine("Seed", params_seed);
    out_file << formattedReportLine("Coordinate initialization method", init_pos_type);
    out_file << formattedReportLine("Number of atoms", natoms);
    out_file << formattedReportLine("Number of beads", nbeads);

    double out_temperature = Units::convertToUser("temperature", "kelvin", temperature);
    out_file << formattedReportLine("Temperature", std::format("{} kelvin", out_temperature));

    double out_sys_size = Units::convertToUser("length", "angstrom", size);
    out_file << formattedReportLine("Linear size of the system", std::format("{} angstroms", out_sys_size));

    double out_mass = Units::convertToUser("mass", "dalton", mass);
    out_file << formattedReportLine("Mass", std::format("{} amu", out_mass));

    out_file << formattedReportLine("Total number of MD steps", steps);
    out_file << formattedReportLine("Interaction potential name", interaction_potential_name);
    out_file << formattedReportLine("External potential name", external_potential_name);

    out_file << "---------\nFeatures\n---------\n";
    out_file << formattedReportLine("Minimum image convention", MINIM);
    out_file << formattedReportLine("Wrapping of coordinates", WRAP);
    out_file << formattedReportLine("Polymer recentering", RECENTER);
    out_file << formattedReportLine("Using i-Pi convention", IPI_CONVENTION);
}

/**
 * Outputs the trajectories of the particles to a xyz file.
 * 
 * @param step Current step of the simulation.
 */
void Simulation::outputTrajectories(int step) {
    std::ofstream xyz_file;
    xyz_file.open(std::format("{}/position_{}.xyz", Output::FOLDER_NAME, this_bead), std::ios::out | std::ios::app);

    xyz_file << std::format("{}\n", natoms);
    xyz_file << std::format(" Atoms. MD step: {}\n", step);

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        xyz_file << "1";

        for (int axis = 0; axis < NDIM; ++axis)
            /// @todo Make the units configurable
            xyz_file << std::format(" {:^20.12e}", Units::convertToUser("length", "angstrom", coord(ptcl_idx, axis)));
#if NDIM == 1
        xyz_file << " 0.0 0.0";
#elif NDIM == 2
        xyz_file << " 0.0";
#endif
        xyz_file << "\n";
    }

    xyz_file.close();
}

/**
 * Outputs the velocities of the particles to a dat file.
 * 
 * @param step Current step of the simulation.
 */
void Simulation::outputVelocities(int step) {
    std::ofstream vel_file;
    vel_file.open(std::format("{}/velocity_{}.dat", Output::FOLDER_NAME, this_bead), std::ios::out | std::ios::app);

    vel_file << std::format("{}\n", natoms);
    vel_file << std::format(" Atoms. Timestep: {}\n", step);

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        vel_file << (ptcl_idx + 1) << " 1";

        for (int axis = 0; axis < NDIM; ++axis)
            /// @todo Make the units configurable
            vel_file << std::format(" {:^20.12e}", Units::convertToUser("velocity", "angstrom/ps", momenta(ptcl_idx, axis) / mass));
#if NDIM == 1
        vel_file << " 0.0 0.0";
#elif NDIM == 2
        vel_file << " 0.0";
#endif
        vel_file << "\n";
    }

    vel_file.close();
}

/**
 * Outputs the forces acting on the particles to a dat file.
 *
 * @param step Current step of the simulation.
 */
void Simulation::outputForces(int step) {
    std::ofstream force_file;
    force_file.open(std::format("{}/force_{}.dat", Output::FOLDER_NAME, this_bead), std::ios::out | std::ios::app);

    force_file << std::format("{}\n", natoms);
    force_file << std::format(" Atoms. Timestep: {}\n", step);

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        force_file << (ptcl_idx + 1) << " 1";

        for (int axis = 0; axis < NDIM; ++axis)
            force_file << std::format(" {:^20.12e}", Units::convertToUser("force", "ev/ang", forces(ptcl_idx, axis)));
#if NDIM == 1
        force_file << " 0.0 0.0";
#elif NDIM == 2
        force_file << " 0.0";
#endif
        force_file << "\n";
    }

    force_file.close();
}

/**
 * @brief Zero the linear momentum of a group of atoms by subtracting the velocity
 * of the center of mass from the velocity of each atom.
 * The calculation assumes that all atoms have the same mass, in which case the
 * center of mass momentum is given by p_c=m*v_c=(p_1+..+p_n)/n, where n=N*P
 * is the total number of beads in the system.
 */
void Simulation::zeroMomentum()
{
    dVec momentum_cm;           // Resulting center of mass momentum vector
    dVec momentum_cm_per_bead;  // Contribution of the current time-slice to the center of mass momentum vector

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            momentum_cm_per_bead(0, axis) += momenta(ptcl_idx, axis);
        }
    }

    for (int axis = 0; axis < NDIM; ++axis) {
        momentum_cm_per_bead(0, axis) /= (natoms * nbeads);
    }

    MPI_Allreduce(momentum_cm_per_bead.data(), momentum_cm.data(), momentum_cm.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            momenta(ptcl_idx, axis) -= momentum_cm(0, axis);
        }
    }
}

/**
 * Initializes the potential based on the input parameters.
 *
 * @param potential_name Name of the potential.
 * @param potential_options Physical parameters of the potential. 
 * @return Pointer to the initialized potential.
 */
std::unique_ptr<Potential> Simulation::initializePotential(const std::string& potential_name, const VariantMap& potential_options)
{
    if (potential_name == "free") {
        return std::make_unique<Potential>();
    }

    if (potential_name == "harmonic") {
        double omega = std::get<double>(potential_options.at("omega"));
        return std::make_unique<HarmonicPotential>(mass, omega);
    }
    
    if (potential_name == "double_well") {
        double strength = std::get<double>(potential_options.at("strength"));
        double loc = std::get<double>(potential_options.at("location"));
        return std::make_unique<DoubleWellPotential>(mass, strength, loc);
    }

    if (potential_name == "dipole") {
        double strength = std::get<double>(potential_options.at("strength"));
        return std::make_unique<DipolePotential>(strength);
    }

    if (potential_name == "aziz") {
        return std::make_unique<AzizPotential>();
    }
    
    return std::make_unique<Potential>();
}

/**
 * Initializes the positions of the particles based on the input parameters.
 *
 * @param[out] coord_arr Array to store the generated positions.
 * @param sim_params Simulation parameters object containing information about the init trajectories file. 
 */
void Simulation::initializePositions(dVec& coord_arr, const VariantMap& sim_params) {
    if (init_pos_type == "xyz") {
        // Load initial Cartesian coordinates from the provided .xyz file
        const std::string xyz_filename = std::get<std::string>(sim_params.at("init_pos_xyz_filename"));

        /// @todo Remove the bead_num==0 condition later.
        /// This is placed for testing purposes, to compare the
        /// simulation results against i-Pi.
        if (this_bead == 0)
            loadTrajectories(xyz_filename, coord_arr);
    } else {
        // Sample positions from a uniform distribution
        genRandomPositions(coord_arr);
    }
}

/**
 * Initializes the momenta of the particles based on the input parameters.
 *
 * @param[out] momentum_arr Array to store the generated momenta.
 */
void Simulation::initializeMomenta(dVec& momentum_arr) {
    if (init_vel_type == "manual") {
        // If loading momenta from elsewhere, do not zero momentum after that
        loadMomenta(std::format("init/vel_{:02}.dat", this_bead + 1), mass, momentum_arr); // LAMMPS convention
        //loadMomenta(std::format("init/vel_{:02}.dat", this_bead), mass, momenta); // i-Pi convention

        /// @todo Handle the cases when the file does not exist or has a wrong format
    } else {
        // If generating momenta from the Maxwell-Boltzmann distribution, zero the total momentum
        genMomentum(momentum_arr);
        zeroMomentum();
    }
}

/**
 * Prints information for debugging purposes.
 * 
 * @param text Text to print.
 */
void Simulation::printDebug(const std::string& text) const
{
    if (this_bead == 0) {
        std::ofstream debug;
        debug.open(std::format("{}/debug.log", Output::FOLDER_NAME), std::ios::out | std::ios::app);
        debug << text;
        debug.close();
    }
}
