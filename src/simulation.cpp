#include <ranges>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <array>
#include <cassert>
#include "mpi.h"

#include "state.h"
#include "observable.h"
#include "simulation.h"

#if OLD_BOSONIC_ALGORITHM
#include "old_bosonic_exchange.h"
#else
#include "bosonic_exchange.h"
#endif

Simulation::Simulation(const int& rank, const int& nproc, Params& param_obj, unsigned int seed) :
    bosonic_exchange(nullptr),
    rand_gen(seed + rank),
    this_bead(rank),
    nproc(nproc) {
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
    getVariant(param_obj.sim["max_wind"], max_wind);

    getVariant(param_obj.sim["apply_mic_spring"], apply_mic_spring);
    getVariant(param_obj.sim["apply_mic_potential"], apply_mic_potential);
    getVariant(param_obj.sim["apply_wrap"], apply_wrap);
    getVariant(param_obj.sim["apply_wrap_first"], apply_wrap_first);
    getVariant(param_obj.sim["apply_wind"], apply_wind);

    /*getVariant(param_obj.out["positions"], out_pos);
    getVariant(param_obj.out["velocities"], out_vel);
    getVariant(param_obj.out["forces"], out_force);
    getVariant(param_obj.out["wind_prob"], out_wind_prob);*/

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

    beta_half_k = beta * 0.5 * spring_constant;

#if IPI_CONVENTION
    beta_half_k /= nbeads;
#endif

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
    initializeMomenta(momenta, param_obj.sim);

    // Initialize the potential based on the input
    external_potential_name = std::get<std::string>(param_obj.external_pot["name"]);
    interaction_potential_name = std::get<std::string>(param_obj.interaction_pot["name"]);
    int_pot_cutoff = std::get<double>(param_obj.interaction_pot["cutoff"]);

    // For cubic cells with PBC, the cutoff distance must be no greater than L/2 for consistency with
    // the minimum image convention (see 1.6.3 in Allen & Tildesley).
    if (pbc) {
        int_pot_cutoff = std::min(int_pot_cutoff, 0.5 * size);
        initializeWindingVectors(wind, max_wind);
    }

    include_wind_corr = pbc && apply_wind && (max_wind > 0);

    ext_potential = initializePotential(external_potential_name, param_obj.external_pot);
    int_potential = initializePotential(interaction_potential_name, param_obj.interaction_pot);

    // Update the coordinate arrays of neighboring particles
    updateNeighboringCoordinates();

    bosonic = bosonic && (nbeads > 1); // Bosonic exchange is only possible for P>1
    is_bosonic_bead = bosonic && (this_bead == 0 || this_bead == nbeads - 1);

    if (is_bosonic_bead) {
#if OLD_BOSONIC_ALGORITHM
        bosonic_exchange = std::make_unique<OldBosonicExchange>(*this);
#else
        bosonic_exchange = std::make_unique<BosonicExchange>(*this);
#endif
    }

    initializeStates(param_obj.states);
    initializeObservables(param_obj.observables);
}

Simulation::~Simulation() = default;

int Simulation::getStep() const {
    return md_step;
}

void Simulation::setStep(int step) {
    md_step = step;
}

/**
 * Generates random positions in the interval [-L/2, L/2] for each particle along each axis.
 * 
 * @param[out] pos_arr Array to store the generated positions.
 */
void Simulation::genRandomPositions(dVec& pos_arr) {
    /// @todo Add ability to generate non-random positions (e.g., lattice)
    std::uniform_real_distribution<double> u_dist(-0.5 * size, 0.5 * size);

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            pos_arr(ptcl_idx, axis) = u_dist(rand_gen);
        }
    }
}

/**
 * Places particles in a uniform grid according to the specified box size.
 *
 * @param[out] pos_arr Array to store the generated positions.
 */
void Simulation::uniformParticleGrid(dVec& pos_arr) const {
    const double volume = std::pow(size, NDIM);

    // Get the linear size per particle, and the number of particles
    const double init_side = std::pow((1.0 * natoms / volume), -1.0 / (1.0 * NDIM));

    // Determine the number and the size of initial grid boxes in each dimension
    int tot_num_grid_boxes = 1;
    iVec num_nn_grid;
    dVec size_nn_grid;

    for (int i = 0; i < NDIM; i++) {
        num_nn_grid(0, i) = static_cast<int>(std::ceil((size / init_side) - EPS));

        // Make sure we have at least one grid box
        if (num_nn_grid(0, i) < 1)
            num_nn_grid(0, i) = 1;

        // Compute the actual size of the grid
        size_nn_grid(0, i) = size / (1.0 * num_nn_grid(0, i));

        // Determine the total number of grid boxes
        tot_num_grid_boxes *= num_nn_grid(0, i);
    }

    // Place the particles at the middle of each box
    if (tot_num_grid_boxes < natoms) {
        throw std::runtime_error("Number of grid boxes is less than the number of particles");
    }

    dVec pos;
    for (int n = 0; n < tot_num_grid_boxes; n++) {
        iVec grid_index;

        for (int i = 0; i < NDIM; i++) {
            int scale = 1;
            for (int j = i + 1; j < NDIM; j++)
                scale *= num_nn_grid(0, j);

            grid_index(0, i) = (n / scale) % num_nn_grid(0, i);
        }

        for (int axis = 0; axis < NDIM; ++axis) {
            pos(0, axis) = (grid_index(0, axis) + 0.5) * size_nn_grid(0, axis) - 0.5 * size;
        }

        /// @todo If wrapping should be applied regardless of PBC, then it can be moved to the previous loop
        if (pbc) {
            for (int axis = 0; axis < NDIM; ++axis) {
                applyMinimumImage(pos(0, axis), size);
            }
        }

        if (n >= natoms) {
            break;
        }

        for (int axis = 0; axis < NDIM; axis++) {
            pos_arr(n, axis) = pos(0, axis);
        }
    }
}

/**
 * Generates momenta according to the Maxwell-Boltzmann distribution.
 * 
 * @param momenta_arr Array to store the generated momenta.
 */
void Simulation::genMomentum(dVec& momenta_arr) {
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
    std::normal_distribution<double> normal(0.0, 1 / sqrt(beta * mass / nbeads));
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

    if (pbc && apply_wrap) {
        for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                applyMinimumImage(coord(ptcl_idx, axis), size);
            }
        }
    }

    if (pbc && apply_wrap_first && this_bead == 0) {
        for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                applyMinimumImage(coord(ptcl_idx, axis), size);
            }
        }
    }

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
    printStatus("Running the simulation", this_bead);

    MPI_Barrier(MPI_COMM_WORLD);
    const double sim_exec_time_start = MPI_Wtime();

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

    for (const auto& state : states) {
        state->initialize();
    }

    // Main loop performing molecular dynamics steps
    for (int step = 0; step <= steps; ++step) {
        setStep(step);

        for (const auto& observable : observables) {
            observable->resetValues();
        }

        for (const auto& state : states) {
            state->output(step);
        }

        // "O" step
        if (enable_t) {
            langevinStep();
        }

        // If fixcom=true, the center of mass of the ring polymers is fixed during the simulation
        if (fixcom) {
            zeroMomentum();
        }

        velocityVerletStep();

        // "O" step
        if (enable_t) {
            langevinStep();
        }

        // Zero momentum after every Langevin step (if needed)
        if (fixcom) {
            zeroMomentum();
        }


#if PROGRESS
            printProgress(step, steps, this_bead);
#endif

        // If we have not reached the thermalization threshold, skip to the next step (thermalization stage)
        if (step < threshold) {
            continue;
        }

        // Calculate and print the observables (production stage)
        for (const auto& observable : observables) {
            observable->calculate();
        }

        /*
        if (is_bosonic_bead) {
            bosonic_exchange->printBosonicDebug();
        }
        */

        if (step % sfreq == 0) {
            if (this_bead == 0) {
                out_file << std::format("{:^16.8e}", static_cast<double>(step));
            }

            for (const auto& observable : observables) {
                // The inner loop is necessary because some observable classes can calculate
                // more than one observable (e.g., "energy" calculates both the kinetic and potential energies).
                for (const double& val : observable->quantities | std::views::values) {
                    double quantity_value = 0.0;
                    double local_quantity_value = val;

                    // Sum the results from all processes (beads)
                    MPI_Allreduce(&local_quantity_value, &quantity_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                    if (this_bead == 0) {
                        out_file << std::format(" {:^16.8e}", quantity_value);
                    }
                }
            }

            if (this_bead == 0) {
                out_file << '\n';
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    const double sim_exec_time_end = MPI_Wtime();

    const double wall_time = sim_exec_time_end - sim_exec_time_start;

    printStatus(std::format("Simulation finished running successfully (Runtime = {:.3} sec)", wall_time), this_bead);

    for (const auto& state : states) {
        state->finalize();
    }

    if (this_bead == 0) {
        out_file.close();

        std::ofstream report_file;
        report_file.open(std::format("{}/report.txt", Output::FOLDER_NAME), std::ios::out | std::ios::app);
        printReport(report_file, wall_time);
        report_file.close();
    }
}

/**
 * Obtains the coordinates from the previous timeslice.
 * 
 * @param prev Vector to store the previous coordinates.
 */
void Simulation::getPrevCoords(dVec& prev) {
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
 * bosonic simulation using the original (inefficient) algorithm, that takes into
 * account all the N! permutations, by setting OLD_BOSONIC_ALGORITHM to true.
 * 
 * @param spring_force_arr Vector to store the spring forces.
 */
void Simulation::updateSpringForces(dVec& spring_force_arr) const {
    if (is_bosonic_bead) {
        // If the simulation is bosonic and the current bead is either 1 or P, we calculate
        // the exterior spring forces in the appropriate bosonic class.
        // Any winding effects that pertain to the exterior springs are taken into account
        // inside the bosonic exchange class.
        bosonic_exchange->prepare();
        bosonic_exchange->exteriorSpringForce(spring_force_arr);
    } else {
        // If particles are distinguishable, or if the current bead is an interior bead,
        // the force is calculated based on the standard expression for distinguishable particles.
        for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                double diff_prev = prev_coord(ptcl_idx, axis) - coord(ptcl_idx, axis);
                double diff_next = next_coord(ptcl_idx, axis) - coord(ptcl_idx, axis);

                if (pbc && apply_mic_spring) {
                    applyMinimumImage(diff_prev, size);
                    applyMinimumImage(diff_next, size);
                }

                spring_force_arr(ptcl_idx, axis) = spring_constant * (diff_prev + diff_next);

                // Winding effects due to interior neighbors are taken into account for all beads
                if (include_wind_corr) {
                    // Only the nonzero winding numbers contribute to the force
                    for (int wind_idx = 1; wind_idx <= max_wind; ++wind_idx) {
                        //spring_force_arr(ptcl_idx, axis) -= spring_constant * size * 
                        //    wind_idx * (getWindingProbability(-diff_next, wind_idx) - getWindingProbability(-diff_next, -wind_idx));

                        // This is equivalent to the above because the winding probability satisfies p(-diff, w) = p(diff, -w)
                        spring_force_arr(ptcl_idx, axis) += spring_constant * size * 
                                wind_idx * (getWindingProbability(diff_next, wind_idx) - getWindingProbability(diff_next, -wind_idx));

                        spring_force_arr(ptcl_idx, axis) += spring_constant * size * 
                                wind_idx * (getWindingProbability(diff_prev, wind_idx) - getWindingProbability(diff_prev, -wind_idx));
                    }
                }
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
    /// @todo Think if MIC is appropriate when the external potential is periodic
    physical_force_arr = (-1.0) * ext_potential->gradV(coord);

    if (int_pot_cutoff != 0.0) {
        for (int ptcl_one = 0; ptcl_one < natoms; ++ptcl_one) {
            for (int ptcl_two = ptcl_one + 1; ptcl_two < natoms; ++ptcl_two) {
                // Get the vector distance between the two particles.
                // Here "diff" contains just one vector of dimension NDIM.
                dVec diff = getSeparation(ptcl_one, ptcl_two, apply_mic_potential);

                // If the distance between the particles exceeds the cutoff length
                // then we assume the interaction is negligible and do not bother
                // calculating the force.
                // We use the convention that when cutoff < 0 then the interaction is
                // calculated for all distances.
                if (const double distance = diff.norm(); distance < int_pot_cutoff || int_pot_cutoff < 0.0) {
                    dVec force_on_one = (-1.0) * int_potential->gradV(diff);

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
 * Calculates the spring energy contribution of the current and the previous time-slice,
 * provided they are classical, i.e., are not affected by bosonic exchange.
 * If the simulation is bosonic, the function is callable only for the interior connections.
 *
 * @return Classical spring energy contribution of the current and the previous time-slice.
 */
double Simulation::classicalSpringEnergy() const {
    assert(!bosonic || (bosonic && this_bead != 0));

    double interior_spring_energy = 0.0;

    if (include_wind_corr) {
        for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
            interior_spring_energy += getWindingEnergyExpectation(prev_coord, ptcl_idx, coord, ptcl_idx);
        }
    } else {
        for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                double diff = prev_coord(ptcl_idx, axis) - coord(ptcl_idx, axis);

                if (pbc && apply_mic_spring) {
                    applyMinimumImage(diff, size);
                }

                interior_spring_energy += diff * diff;
            }
        }

        interior_spring_energy *= 0.5 * spring_constant;
    }

    return interior_spring_energy;
}

/**
 * Returns the vectorial distance between two particles at the same imaginary timeslice.
 * Mathematically equivalent to r1-r2, where r1 and r2 are the position-vectors
 * of the first and second particles, respectively. In the case of PBC, there is an
 * option to return the minimal distance between the two particles.
 * 
 * @param first_ptcl Index of the first particle.
 * @param second_ptcl Index of the second particle.
 * @param minimum_image Flag determining whether the minimum image convention should be applied.
 * @return Vectorial distance between the two particles.
 */
dVec Simulation::getSeparation(int first_ptcl, int second_ptcl, bool minimum_image) const {
    dVec sep;

    for (int axis = 0; axis < NDIM; ++axis) {
        double diff = coord(first_ptcl, axis) - coord(second_ptcl, axis);

        if (pbc && minimum_image)
            applyMinimumImage(diff, size);
        sep(0, axis) = diff;
    }

    return sep;
}

/**
 * @brief Prints a summary of the simulation parameters at the end of the simulation.
 */
void Simulation::printReport(std::ofstream& out_file, double wall_time) const {
    out_file << "---------\nParameters\n---------\n";

    if (bosonic) {
        out_file << formattedReportLine("Statistics", "Bosonic");
        std::string bosonic_alg_name = "Feldman-Hirshberg [O(N^2) scaling]";

#if OLD_BOSONIC_ALGORITHM
        bosonic_alg_name = "Primitive [O(N!) scaling]";
#endif

        out_file << formattedReportLine("Bosonic algorithm", bosonic_alg_name);
    } else {
        out_file << formattedReportLine("Statistics", "Boltzmannonic");
    }

    out_file << formattedReportLine("Periodic boundary conditions", pbc);
    if (pbc) {
        out_file << formattedReportLine("Winding cutoff", max_wind);
    }
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
    //out_file << formattedReportLine("Minimum image convention", MINIM);
    //out_file << formattedReportLine("Wrapping of coordinates", WRAP);
    //out_file << formattedReportLine("Wrapping of the first time-slice", WRAP_FIRST);
    out_file << formattedReportLine("Minimum image convention (springs)", apply_mic_spring);
    out_file << formattedReportLine("Minimum image convention (potential)", apply_mic_potential);
    out_file << formattedReportLine("Coordinate wrap (all time-slices)", apply_wrap);
    out_file << formattedReportLine("Coordinate wrap (1st time-slice only)", apply_wrap_first);
    out_file << formattedReportLine("Winding effects", apply_wind);

    out_file << formattedReportLine("Polymer recentering", RECENTER);
    out_file << formattedReportLine("Using i-Pi convention", IPI_CONVENTION);

    out_file << "---------\n";
    out_file << formattedReportLine("Wall time", std::format("{:%T}", 
        std::chrono::duration<double>(wall_time)
    ));
}

/**
 * @brief Zero the linear momentum of a group of atoms by subtracting the velocity
 * of the center of mass from the velocity of each atom.
 * The calculation assumes that all atoms have the same mass, in which case the
 * center of mass momentum is given by p_c=m*v_c=(p_1+...+p_n)/n, where n=N*P
 * is the total number of beads in the system.
 */
void Simulation::zeroMomentum() {
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

    MPI_Allreduce(momentum_cm_per_bead.data(), momentum_cm.data(), momentum_cm.size(), MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

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
std::unique_ptr<Potential> Simulation::initializePotential(const std::string& potential_name,
                                                           const VariantMap& potential_options) {
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

    if (potential_name == "cosine") {
        double amplitude = std::get<double>(potential_options.at("amplitude"));
        double phase = std::get<double>(potential_options.at("phase"));
        return std::make_unique<CosinePotential>(amplitude, size, phase);
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

        // Uncomment the next line to restrict the coordinate initialization to the first time-slice
        //if (this_bead == 0)
        loadTrajectories(xyz_filename, coord_arr);
    } else if (init_pos_type == "xyz_formatted") {
        const std::string xyz_filename_format = std::get<std::string>(sim_params.at("init_pos_xyz_filename_format"));
        const int arg = this_bead + std::get<int>(sim_params.at("init_pos_first_index"));
        
        loadTrajectories(std::vformat(xyz_filename_format, std::make_format_args(arg)), coord_arr);
    } else if (init_pos_type == "grid") {
        // Generate a grid of particles
        uniformParticleGrid(coord_arr);
    } else {
        // Sample positions from a uniform distribution
        genRandomPositions(coord_arr);
    }
}

/**
 * Initializes the momenta of the particles based on the input parameters.
 *
 * @param[out] momentum_arr Array to store the generated momenta.
 * @param sim_params Simulation parameters object containing information about the init velocities file. 
 */
void Simulation::initializeMomenta(dVec& momentum_arr, const VariantMap& sim_params) {
    if (init_vel_type == "manual") {
        // If loading momenta from elsewhere, do not zero momentum after that
        //loadMomenta(std::format("init/vel_{:02}.dat", this_bead + 1), mass, momentum_arr); // LAMMPS convention
        loadMomenta(std::format("init/vel_{:02}.dat", this_bead), mass, momenta); // i-Pi convention

        /// @todo Handle the cases when the file does not exist or has a wrong format
    } else if (init_vel_type == "manual_formatted") {
        const std::string vel_filename_format = std::get<std::string>(sim_params.at("init_vel_manual_filename_format"));
        const int arg = this_bead + std::get<int>(sim_params.at("init_vel_first_index"));
        
        loadMomenta(std::vformat(vel_filename_format, std::make_format_args(arg)), mass, momentum_arr);
    } else {
        // If generating momenta from the Maxwell-Boltzmann distribution, zero the total momentum
        genMomentum(momentum_arr);
        zeroMomentum();
    }
}

/**
 * Initializes a state based on the input parameters. Initialization occurs only if the
 * state is enabled, i.e., if the units are not set to "off" or "false".
 *
 * @param sim_params Simulation parameters object containing information about the states.
 * @param param_key Key of the parameter in the simulation parameters object.
 * @param state_name Name of the state.
 */
void Simulation::addStateIfEnabled(const StringMap& sim_params, const std::string& param_key, const std::string& state_name) {
    if (const std::string& units = sim_params.at(param_key); units != "off" && units != "false") {
        if (units == "none") {
            states.push_back(StateFactory::createQuantity(state_name, *this, sfreq, ""));
        } else if (units == "on" || units == "true" || units == "auto") {
            /// @todo What if the state corresponds to dimensionless quantity? Shouldn't auto be the same as none?
            states.push_back(StateFactory::createQuantity(state_name, *this, sfreq, "atomic_unit"));
        } else {
            states.push_back(StateFactory::createQuantity(state_name, *this, sfreq, units));
        }
    }
}

/**
 * Method for initializing all the requested states.
 *
 * @param sim_params Simulation parameters object containing information about the states.
 */
void Simulation::initializeStates(const StringMap& sim_params) {
    addStateIfEnabled(sim_params, "positions", "position");
    addStateIfEnabled(sim_params, "velocities", "velocity");
    addStateIfEnabled(sim_params, "forces", "force");
    addStateIfEnabled(sim_params, "wind_prob", "wind_prob");
}

/**
 * Initializes an observable based on the input parameters. Initialization occurs only if the
 * observable is enabled (the units are not set to "off").
 *
 * @param sim_params Simulation parameters object containing information about the observables.
 * @param param_key Key of the parameter in the simulation parameters object.
 * @param observable_name Name of the observable.
 */
void Simulation::addObservableIfEnabled(const StringMap& sim_params, const std::string& param_key, const std::string& observable_name) {
    if (const std::string& units = sim_params.at(param_key); units != "off") {
        if (units == "none") {
            observables.push_back(ObservableFactory::createQuantity(observable_name, *this, sfreq, ""));
        } else {
            observables.push_back(ObservableFactory::createQuantity(observable_name, *this, sfreq, units));
        }
    }
}

/**
 * Method for initializing all the requested observables.
 *
 * @param sim_params Simulation parameters object containing information about the observables.
 */
void Simulation::initializeObservables(const StringMap& sim_params) {
    addObservableIfEnabled(sim_params, "energy", "energy");
    addObservableIfEnabled(sim_params, "classical", "classical");

    if (bosonic) {
        addObservableIfEnabled(sim_params, "bosonic", "bosonic");
    }
}

/**
 * Returns the natural logarithm of the sum of Boltzmann exponents over all the winding vectors for a pair of beads.
 *
 * @param left_x Coordinate vectors of left time-slice (e.g., P).
 * @param left_idx Particle index at the left time-slice.
 * @param right_x Coordinate vectors of right time-slice (e.g., 1).
 * @param right_idx Particle index at the right time-slice.
 * @return Log of sum of Boltzmann exponents.
 */
double Simulation::getLogWindingWeight(const dVec& left_x, int left_idx, const dVec& right_x, int right_idx) const {
    double result = 0.0;

    for (int axis = 0; axis < NDIM; ++axis) {
        const double diff = left_x(left_idx, axis) - right_x(right_idx, axis);
        const double shift = getWindingShift(diff);
        double weight = exp(-beta_half_k * (diff * diff - shift)); // Zero winding contribution

        for (int wind_num = 1; wind_num <= max_wind; ++wind_num) {
            const double diff_plus = diff + wind_num * size;
            const double diff_minus = diff - wind_num * size;
            weight += exp(-beta_half_k * (diff_plus * diff_plus - shift)) + exp(-beta_half_k * (diff_minus * diff_minus - shift));
        }

        result += log(weight) - beta_half_k * shift;
    }

    return result;
}

/**
 * Returns the expectation of the spring energy with respect to the winding probability.
 *
 * @param left_x Coordinate vectors of left time-slice (e.g., P).
 * @param left_idx Particle index at the left time-slice.
 * @param right_x Coordinate vectors of right time-slice (e.g., 1).
 * @param right_idx Particle index at the right time-slice.
 * @return Expectation of the spring energy.
 */
double Simulation::getWindingEnergyExpectation(const dVec& left_x, int left_idx, const dVec& right_x, int right_idx) const {
    double total = 0.0;

    for (int axis = 0; axis < NDIM; ++axis) {
        const double diff = left_x(left_idx, axis) - right_x(right_idx, axis);
        total += diff * diff * getWindingProbability(diff, 0);

        for (int wind_num = 1; wind_num <= max_wind; ++wind_num) {
            const double diff_plus = diff + wind_num * size;
            const double diff_minus = diff - wind_num * size;

            total += diff_plus * diff_plus * getWindingProbability(diff, wind_num);
            total += diff_minus * diff_minus * getWindingProbability(diff, -wind_num);
        }
    }

    return 0.5 * spring_constant * total;
}

/**
 * Calculates the position squared shift for numerical stability of the winding probabilities
 *
 * @param diff Difference between the positions of the two particles.
 * @return Position squared shift.
 */
double Simulation::getWindingShift(const double diff) const {
    // Start from the zero winding
    double shift = std::min(std::numeric_limits<double>::max(), diff * diff);

    for (int wind_idx = 1; wind_idx <= max_wind; ++wind_idx) {
        const double diff_plus_wind = diff + wind_idx * size;
        const double diff_minus_wind = diff - wind_idx * size;

        shift = std::min(shift, diff_plus_wind * diff_plus_wind);
        shift = std::min(shift, diff_minus_wind * diff_minus_wind);
    }

    return shift;
}

/**
 * Returns the probability of a particular winding number along some axis for a pair of beads.
 *
 * @param diff Distance between the beads along some axis.
 * @param winding_number Value of the winding number.
 * @return Probability of attaing a configuration with the specified winding number.
 */
double Simulation::getWindingProbability(const double diff, const int winding_number) const {
    // Shift for numerical stability
    const double shift = getWindingShift(diff);

    // Start with the zero winding contribution
    double denominator = exp(-beta_half_k * (diff * diff - shift));

    // Add the contribution of the nonzero winding numbers
    for (int wind_idx = 1; wind_idx <= max_wind; ++wind_idx) {
        const double diff_plus = diff + wind_idx * size;
        const double diff_minus = diff - wind_idx * size;

        denominator += exp(-beta_half_k * (diff_plus * diff_plus - shift));
        denominator += exp(-beta_half_k * (diff_minus * diff_minus - shift));
    }

    const double diff_val = diff + winding_number * size;

    return exp(-beta_half_k * (diff_val * diff_val - shift)) / denominator;
}

/**
 * Initializes the winding vectors for periodic boundary conditions.
 * Currently used only for debugging purposes.
 *
 * @param[out] wind_arr The array to store the winding vectors.
 * @param wind_cutoff The cutoff for the winding vectors.
 */
void Simulation::initializeWindingVectors(iVec& wind_arr, int wind_cutoff) {
    // Number of different winding vectors for a given cutoff
    const int num_wind = static_cast<int>(std::pow(2 * wind_cutoff + 1, NDIM));
    wind_arr = iVec(num_wind);

    std::array<int, NDIM> current = {};
    current.fill(-wind_cutoff);  // Initialize with minimum value

    int current_idx = 0; // Index of the current winding vector

    // Generate all the possible winding vectors for the provided cutoff
    while (true) {
        // Add the current vector to the result
        for (int axis = 0; axis < NDIM; ++axis) {
            wind_arr(current_idx, axis) = current[axis];
        }
        current_idx++;

        int index = NDIM - 1;
        while (index >= 0 && current[index] == wind_cutoff) {
            current[index] = -wind_cutoff;
            index--;
        }

        // If all components are at their maximum value, break the loop
        if (index < 0)
            break;

        current[index]++;
    }
}

/**
 * Prints information for debugging purposes.
 * 
 * @param text Text to print.
 * @param target_bead Bead for which the text should be printed.
 */
void Simulation::printDebug(const std::string& text, int target_bead) const {
    if (this_bead == target_bead) {
        std::ofstream debug;
        debug.open(std::format("{}/debug.log", Output::FOLDER_NAME), std::ios::out | std::ios::app);
        debug << text;
        debug.close();
    }
}
