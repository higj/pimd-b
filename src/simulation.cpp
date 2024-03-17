#include "observable.h"
#include "simulation.h"

#if OLD_BOSONIC_ALGORITHM
#include "old_bosonic_exchange.h"
#else
#include "bosonic_exchange.h"
#endif

Simulation::Simulation(int& rank, int& nproc, Params& paramObj, unsigned int seed) :
    this_bead(rank),
    nproc(nproc),
    rand_gen(seed + rank),
    bosonic_exchange(nullptr)
{
    getVariant(paramObj.sim["nbeads"], nbeads);
    getVariant(paramObj.sim["dt"], dt);
    getVariant(paramObj.sim["sfreq"], sfreq);
    getVariant(paramObj.sim["steps"], steps);
    getVariant(paramObj.sim["enable_t"], enable_t);

    getVariant(paramObj.sim["threshold"], threshold);
    threshold *= steps;  // Threshold in terms of steps

    getVariant(paramObj.sim["gamma"], gamma);
    getVariant(paramObj.sim["bosonic"], bosonic);
    getVariant(paramObj.sim["fixcom"], fixcom);
    getVariant(paramObj.sim["pbc"], pbc);

    getVariant(paramObj.out["positions"], out_pos);
    getVariant(paramObj.out["velocities"], out_vel);
    getVariant(paramObj.out["forces"], out_force);

    getVariant(paramObj.sys["temperature"], temperature);
    getVariant(paramObj.sys["natoms"], natoms);
    getVariant(paramObj.sys["size"], size);

    beta = 1.0 / (Constants::kB * temperature);

    getVariant(paramObj.sys["mass"], mass);

#if IPI_CONVENTION
    // i-Pi convention [J. Chem. Phys. 133, 124104 (2010)]
    omega_p = nbeads / (beta * Constants::hbar);
#else
    // Tuckerman convention
    omega_p = sqrt(nbeads) / (beta * Constants::hbar);
#endif

    spring_constant = mass * omega_p * omega_p;

    // Get the seed from the config file
    getVariant(paramObj.sim["seed"], params_seed);

    // Based on the provided seed, make sure each process has a unique seed
    rand_gen.seed(params_seed + rank);
    mars_gen = std::make_unique<RanMars>(params_seed + rank);

    init_pos_type = std::get<std::string>(paramObj.sim["init_pos_type"]);
    init_vel_type = std::get<std::string>(paramObj.sim["init_vel_type"]);
    propagator_type = std::get<std::string>(paramObj.sim["propagator_type"]);
    
    // Choose time propagation scheme
    if (propagator_type == "cartesian") {
        propagator = std::make_unique<VelocityVerletPropagator>(*this);
    } else if (propagator_type == "normal_modes") {
        propagator = std::make_unique<NormalModesPropagator>(*this);
    }
    
    // Generate initial conditions
    coord = dVec(natoms);
    prev_coord = dVec(natoms);
    next_coord = dVec(natoms);
    momenta = dVec(natoms);

    if (init_pos_type == "xyz") {
        // Load initial Cartesian coordinates from the provided .xyz file
        std::string xyz_filename = std::get<std::string>(paramObj.sim["init_pos_xyz_filename"]);
        loadXYZ(xyz_filename, coord);
    } else {
        // Generate random positions (sample from a uniform distribution)
        genRandomPositions(coord);
    }

    if (init_vel_type == "manual") {
        // If loading momenta from elsewhere, do not zero momentum after that
        loadMomenta(std::format("init/vel_{:02}.dat", this_bead + 1), mass, momenta); // LAMMPS convention
        //loadMomenta(std::format("init/vel_{:02}.dat", this_bead), mass, momenta); // i-Pi convention

        // TODO: Throw error if the file doesn't exist or has a wrong format
    } else {
        // If generating momenta from the Maxwell-Boltzmann distribution, zero the total momentum
        genMomentum(momenta);
        zeroMomentum();
    }

    // Initialize forces to zeros
    forces = dVec(natoms);

    // Initialize the potential based on the input
    // TODO: Make it customizable from config
    external_potential_name = std::get<std::string>(paramObj.external_pot["name"]);
    interaction_potential_name = std::get<std::string>(paramObj.interaction_pot["name"]);

    if (external_potential_name == "free") {
        ext_potential = std::make_unique<Potential>();
    } else if (external_potential_name == "harmonic") {
        double omega = std::get<double>(paramObj.external_pot["omega"]);
        ext_potential = std::make_unique<HarmonicPotential>(mass, omega);
    } else if (external_potential_name == "double_well") {
        double strength = std::get<double>(paramObj.external_pot["strength"]);
        double loc = std::get<double>(paramObj.external_pot["location"]);
        ext_potential = std::make_unique<DoubleWellPotential>(mass, strength, loc);
    } else {
        ext_potential = std::make_unique<Potential>();
    }

    int_pot_cutoff = std::get<double>(paramObj.interaction_pot["cutoff"]);

    // For cubic cells with PBC, the cutoff distance must be no greater than L/2 for consistency with
    // the minimum image convention (see 1.6.3 in Allen & Tildesley).
    if (pbc)
        int_pot_cutoff = std::min(int_pot_cutoff, 0.5 * size);

    if (interaction_potential_name == "free") {
        int_potential = std::make_unique<Potential>();
        int_pot_cutoff = 0.0;
    }
    else if (interaction_potential_name == "harmonic") {
        double omega = std::get<double>(paramObj.interaction_pot["omega"]);
        int_potential = std::make_unique<HarmonicPotential>(mass, omega);
    }
    else if (interaction_potential_name == "double_well") {
        double strength = std::get<double>(paramObj.interaction_pot["strength"]); // TODO: Check the units of "strength"
        double loc = std::get<double>(paramObj.interaction_pot["location"]);
        int_potential = std::make_unique<DoubleWellPotential>(mass, strength, loc);
    }
    else if (interaction_potential_name == "dipole") {
        double strength = std::get<double>(paramObj.interaction_pot["strength"]); // TODO: Check the units of "strength"
        int_potential = std::make_unique<DipolePotential>(strength);
    }
    else {
        int_potential = std::make_unique<Potential>();
    }

    /***** Observables *****/
    ObservableFactory obs_factory;
    observables.push_back(obs_factory.createQuantity("energy", *this, sfreq, "kelvin"));

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

// Generate random positions in the interval [-L/2, L/2], along each axis.
// TODO: Add ability to generate non-random positions (e.g., lattice)
void Simulation::genRandomPositions(dVec &pos_arr) {
    std::uniform_real_distribution<double> u_dist(-0.5 * size, 0.5 * size);

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            pos_arr(ptcl_idx, axis) = u_dist(rand_gen);
        }
    }
}

// Generate momenta from the Maxwell-Boltzmann distribution
void Simulation::genMomentum(dVec &momenta_arr) {
    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            momenta_arr(ptcl_idx, axis) = mass * sampleMB();
        }
    }
}

// Sample the Maxwell-Boltzmann distribution
double Simulation::sampleMB() {
#if IPI_CONVENTION
    std::normal_distribution<double> normal(0.0, 1/sqrt(beta * mass / nbeads));
#else
    std::normal_distribution<double> normal(0.0, 1 / sqrt(beta * mass));
#endif

    return normal(rand_gen);    
}

// Perform a single Langevin step
void Simulation::langevinStep() {
    std::normal_distribution<double> normal;

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

void Simulation::vvStep() {
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

#if RECENTER
    // Recentering should be attempted only in the case of periodic boundary conditions.
    // Also, with the current implementation, it can only work for distinguishable particles.
    if (pbc && !bosonic) {
        // If the initial bead has moved outside of the fundamental cell,
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

// Perform an MD run using the OBABO scheme
void Simulation::run() {
    std::filesystem::create_directory(OUTPUT_FOLDER_NAME);
    std::ofstream out_file;

    if (this_bead == 0) {
        out_file.open(std::format("{}/simulation.out", OUTPUT_FOLDER_NAME, this_bead), std::ios::out | std::ios::app);

        /**** Header ****/
        out_file << std::format("{:^16s}", "step");

        for (const auto& observable : observables) {
            for (const auto& elem : observable->quantities) {
                out_file << std::vformat(" {:^16s}", std::make_format_args(elem.first));
            }
        }

        out_file << "\n";
    }
    
    /**** Main loop performing molecular dynamics steps ****/
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

        //vvStep();
        propagator->step();

        // "O" step
        if (enable_t)
            langevinStep();

        // If fixcom=true, the center of mass of the ring polymers is fixed during the simulation
        if (fixcom)
            zeroMomentum();

#if PROGRESS
            printProgress(step, steps, this_bead);
#endif

        // If we have not reached the thermalization threshold, skip to the next step
        if (step < threshold)
            continue;

        // Calculate and print the observables
        for (const auto& observable : observables) {
            observable->calculate();
        }

        if (step % sfreq == 0) {
            if (this_bead == 0)
                out_file << std::format("{:^16.8e}", static_cast<double>(step));

            for (const auto& observable : observables) {
                // The inner loop is necessary because some observable classes can calculate
                // more than one observable (e.g., "energy" calculates both the kinetic and potential energies).
                for (const auto& elem : observable->quantities) {
                    double quantity_value = 0.0;
                    double local_quantity_value = elem.second;

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
        report_file.open(std::format("{}/report.txt", OUTPUT_FOLDER_NAME), std::ios::out | std::ios::app);
        printReport(report_file);
        report_file.close();
    }
}

// Force the provided initial conditions
void Simulation::forceIniCond(dVec pos_arr, dVec momentum_arr) {
    coord = pos_arr;
    momenta = momentum_arr;
}

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

void Simulation::updateForces() {
    // We distinguish between two types of forces: spring forces between the beads (due to the 
    // classical isomorphism) and the non-spring forces (due to either an external potential 
    // or interactions).
    dVec spring_forces(natoms), physical_forces(natoms);

    updateSpringForces(spring_forces);

    // Calculate the external forces acting on the particles
    physical_forces = (-1.0) * ext_potential->gradV(coord);

    if (int_pot_cutoff != 0.0) {
        for (int ptcl_one = 0; ptcl_one < natoms; ++ptcl_one) {
            for (int ptcl_two = ptcl_one + 1; ptcl_two < natoms; ++ptcl_two) {
                dVec diff = getSeparation(ptcl_one, ptcl_two);  // Vectorial distance
                double distance = diff.norm(0);                 // Scalar distance

                // If the distance between the particles exceeds the cutoff length
                // then we assume the interaction is negligible and do not bother
                // calculating the force.
                // We use the convention that when cutoff < 0 then the interaction is
                // calculated for all distances.
                if (distance < int_pot_cutoff || int_pot_cutoff < 0.0) {
                    dVec force_on_one(1);
                    force_on_one = (-1.0) * int_potential->gradV(diff);

                    for (int axis = 0; axis < NDIM; ++axis) {
                        physical_forces(ptcl_one, axis) += force_on_one(0, axis);
                        physical_forces(ptcl_two, axis) -= force_on_one(0, axis);
                    }
                }
            }
        }
    }

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

void Simulation::updateNeighboringCoordinates() {
    getPrevCoords(prev_coord);
    getNextCoords(next_coord);
}

// Updates the force vector exerted on a specific bead by the two neighboring beads.
// In the distinguishable case, the force is given by Eqn. (12.6.4) in Tuckerman (1st ed).
// In the bosonic case, by default, the forces are evaluated using the algorithm 
// described in https://doi.org/10.1063/5.0173749. It is also possible to perform the
// bosonic simulation using the original (ineffcient) algorithm, that takes into
// account all the N! permutations, by setting OLD_BOSONIC_ALGORITHM to true.
void Simulation::updateSpringForces(dVec& spring_force_arr) {
    if (bosonic) {
        bosonic_exchange->updateCoordinates(coord, prev_coord, next_coord);
        bosonic_exchange->spring_force(spring_force_arr);
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

// Return the vectorial distance between two particles at the same imaginary timeslice.
// Mathematically equivalent to r1-r2, where r1 and r2 are the position-vectors
// of the first and second particles, respectively.
dVec Simulation::getSeparation(int first_ptcl, int second_ptcl) const {
    dVec sep(1);
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

// Print a summary of simulation parameters at the end of the simulation.
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
    
    out_file << formattedReportLine("Time propagation algorithm", propagator_type);
    out_file << formattedReportLine("PBC", pbc);
    out_file << formattedReportLine("Dimension", NDIM);
    out_file << formattedReportLine("Seed", params_seed);
    out_file << formattedReportLine("Coordinate initialization method", init_pos_type);
    out_file << formattedReportLine("Number of atoms", natoms);
    out_file << formattedReportLine("Number of beads", nbeads);

    double out_temperature = Units::unit_to_user("temperature", "kelvin", temperature);
    out_file << formattedReportLine("Temperature", std::format("{} kelvin", out_temperature));

    double out_sys_size = Units::unit_to_user("length", "angstrom", size);
    out_file << formattedReportLine("Linear size of the system", std::format("{} angstroms", out_sys_size));

    double out_mass = Units::unit_to_user("mass", "dalton", mass);
    out_file << formattedReportLine("Mass", std::format("{} amu", out_mass));

    out_file << formattedReportLine("Total number of MD steps", steps);
    out_file << formattedReportLine("Interaction potential name", interaction_potential_name);
    out_file << formattedReportLine("External potential name", external_potential_name);
    //out_file << "---------\n";
    //out_file << "Note: All dimensional quantites are printed in atomic units.";

    out_file << "---------\nFeatures\n---------\n";
    out_file << formattedReportLine("Minimum image convention", MINIM);
    out_file << formattedReportLine("Wrapping of coordinates", WRAP);
    out_file << formattedReportLine("Polymer recentering", RECENTER);
    out_file << formattedReportLine("Using i-Pi convention", IPI_CONVENTION);
}

void Simulation::outputTrajectories(int step) {
    std::ofstream xyz_file;
    xyz_file.open(std::format("{}/position_{}.xyz", OUTPUT_FOLDER_NAME, this_bead), std::ios::out | std::ios::app);

    xyz_file << std::format("{}\n", natoms);
    xyz_file << std::format(" Atoms. MD step: {}\n", step);

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        xyz_file << "1";

        for (int axis = 0; axis < NDIM; ++axis)
            // TODO: Make the units configurable
            xyz_file << " " << Units::unit_to_user("length", "angstrom", coord(ptcl_idx, axis));
#if NDIM == 1
        xyz_file << " 0.0 0.0";
#elif NDIM == 2
        xyz_file << " 0.0";
#endif
        xyz_file << "\n";
    }

    xyz_file.close();
}

void Simulation::outputVelocities(int step) {
    std::ofstream vel_file;
    vel_file.open(std::format("{}/velocity_{}.dat", OUTPUT_FOLDER_NAME, this_bead), std::ios::out | std::ios::app);

    vel_file << std::format("{}\n", natoms);
    vel_file << std::format(" Atoms. Timestep: {}\n", step);

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        vel_file << (ptcl_idx + 1) << " 1";

        for (int axis = 0; axis < NDIM; ++axis)
            // TODO: Make the units configurable
            vel_file << " " << Units::unit_to_user("velocity", "angstrom/ps", momenta(ptcl_idx, axis) / mass);
#if NDIM == 1
        vel_file << " 0.0 0.0";
#elif NDIM == 2
        vel_file << " 0.0";
#endif
        vel_file << "\n";
    }

    vel_file.close();
}

void Simulation::outputForces(int step) {
    std::ofstream force_file;
    force_file.open(std::format("{}/force_{}.dat", OUTPUT_FOLDER_NAME, this_bead), std::ios::out | std::ios::app);

    force_file << std::format("{}\n", natoms);
    force_file << std::format(" Atoms. Timestep: {}\n", step);

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        force_file << (ptcl_idx + 1) << " 1";

        for (int axis = 0; axis < NDIM; ++axis)
            force_file << " " << Units::unit_to_user("force", "ev/ang", forces(ptcl_idx, axis));
#if NDIM == 1
        force_file << " 0.0 0.0";
#elif NDIM == 2
        force_file << " 0.0";
#endif
        force_file << "\n";
    }

    force_file.close();
}

// Zero the linear momentum of a group of atoms by subtracting the velocity
// of the center of mass from the velocity of each atom.
void Simulation::zeroMomentum() {
    dVec momentum_cm(1);
    const int cm_size = momentum_cm.size();

    if (this_bead == 0) {
        double denom = natoms;
        dVec momentum_cm_per_bead(1);

        for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                momentum_cm(0, axis) += momenta(ptcl_idx, axis);
            }
        }

        for (int axis = 0; axis < NDIM; ++axis) {
            momentum_cm(0, axis) /= denom;
        }
    }

    MPI_Bcast(momentum_cm.data(), cm_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            momenta(ptcl_idx, axis) -= momentum_cm(0, axis);
        }
    }
}

void Simulation::printDebug(const std::string& text) {
    if (this_bead == 0) {
        std::ofstream debug;
        debug.open(std::format("{}/debug.log", OUTPUT_FOLDER_NAME), std::ios::out | std::ios::app);
        debug << text;
        debug.close();
    }
}