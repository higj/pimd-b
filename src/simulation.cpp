#include <ranges>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <array>
#include <cassert>
#include "mpi.h"
#include "states.h"
#include "observables.h"
#include "propagators.h"
#include "thermostats.h"
#include "normal_modes.h"
#include "simulation.h"

Simulation::Simulation(const int& rank, const int& nproc, Params& param_obj, unsigned int seed) :
    bosonic_exchange(nullptr),
    rand_gen(seed + rank),
    this_bead(rank),
    nproc(nproc) {
    getVariant(param_obj.sim["nbeads"], nbeads);
    getVariant(param_obj.sim["dt"], dt);
    getVariant(param_obj.sim["sfreq"], sfreq);
    getVariant(param_obj.sim["steps"], steps);

    getVariant(param_obj.sim["threshold"], threshold);
    threshold *= steps;  // Threshold in terms of steps

    getVariant(param_obj.sim["gamma"], gamma);
    getVariant(param_obj.sim["nchains"], nchains);
    getVariant(param_obj.sim["bosonic"], bosonic);
    getVariant(param_obj.sim["fixcom"], fixcom);
    getVariant(param_obj.sim["pbc"], pbc);
    getVariant(param_obj.sim["nmthermostat"], nmthermostat);

    getVariant(param_obj.sys["temperature"], temperature);
    getVariant(param_obj.sys["natoms"], natoms);
    getVariant(param_obj.sys["size"], size);

    beta = 1.0 / (Constants::kB * temperature);

    getVariant(param_obj.sys["mass"], mass);

#if IPI_CONVENTION
    // i-Pi convention [J. Chem. Phys. 133, 124104 (2010)]
    thermo_beta = beta / nbeads;
    omega_p = nbeads / (beta * Constants::hbar);
#else
    // Tuckerman convention
    thermo_beta = beta;
    omega_p = sqrt(nbeads) / (beta * Constants::hbar);
#endif

    spring_constant = mass * omega_p * omega_p;
    beta_half_k = thermo_beta * 0.5 * spring_constant;

    // Get the seed from the config file
    getVariant(param_obj.sim["seed"], params_seed);

    // Based on the provided seed, make sure each process has a unique seed
    rand_gen.seed(params_seed + rank);
    mars_gen = std::make_unique<RanMars>(params_seed + rank);

    init_pos_type = std::get<std::string>(param_obj.sim["init_pos_type"]);
    init_vel_type = std::get<std::string>(param_obj.sim["init_vel_type"]);
    propagator_type = std::get<std::string>(param_obj.sim["propagator_type"]);
    thermostat_type = std::get<std::string>(param_obj.sim["thermostat_type"]);
    
    if (propagator_type == "cartesian") {
        propagator = std::make_unique<VelocityVerletPropagator>(*this);
    } else if (propagator_type == "normal_modes") {
        propagator = std::make_unique<NormalModesPropagator>(*this);
    }

    if (thermostat_type == "langevin") {
        thermostat = std::make_unique<LangevinThermostat>(*this, nmthermostat);
    } else if (thermostat_type == "nose_hoover") {
        thermostat = std::make_unique<NoseHooverThermostat>(*this, nmthermostat, nchains);
    } else if (thermostat_type == "nose_hoover_np") {
        thermostat = std::make_unique<NoseHooverNpThermostat>(*this, nmthermostat, nchains);
    } else if (thermostat_type == "nose_hoover_np_dim") {
        thermostat = std::make_unique<NoseHooverNpDimThermostat>(*this, nmthermostat, nchains);
    } else if (thermostat_type == "none") {
        thermostat = std::make_unique<Thermostat>(*this, nmthermostat);
    }
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

    // If the interaction potential is set to "free", then the cutoff distance is meaningless
    int_pot_cutoff = (interaction_potential_name == "free") ? 0.0 : std::get<double>(param_obj.interaction_pot["cutoff"]);

    // For cubic cells with PBC, the cutoff distance must be no greater than L/2 for consistency with
    // the minimum image convention (see 1.6.3 in Allen & Tildesley).
    if (pbc) {
        int_pot_cutoff = std::min(int_pot_cutoff, 0.5 * size);
    }

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
    normal_modes = std::make_unique<NormalModes>(*this);
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
    std::normal_distribution<double> normal(0.0, 1 / sqrt(thermo_beta * mass));

    return normal(rand_gen);
}

/**
 * @brief Perform a molecular dynamics run using the OBABO scheme.
 */
void Simulation::run() {
    printStatus("Running the simulation", this_bead);
    
    MPI_Barrier(MPI_COMM_WORLD);
    const double sim_exec_time_start = MPI_Wtime();
    
    std::filesystem::create_directory(Output::FOLDER_NAME);
    ObservablesLogger obs_logger(Output::MAIN_FILENAME, this_bead, observables);

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
        thermostat->step();

        // If fixcom=true, the center of mass of the ring polymers is fixed during the simulation
        if (fixcom) {
            zeroMomentum();
        }

        propagator->step();
        thermostat->step();

        // Zero momentum after every thermostat step (if needed)
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

        // Calculate the observables (production stage)
        for (const auto& observable : observables) {
            observable->calculate();
        }

        // Save the observables at the specified frequency
        if (step % sfreq == 0) {
            obs_logger.log(step);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    const double sim_exec_time_end = MPI_Wtime();

    const double wall_time = sim_exec_time_end - sim_exec_time_start;

    printStatus(std::format("Simulation finished running successfully (Runtime = {:.3} sec)", wall_time), this_bead);

    printReport(wall_time);
}


/**
 * Receives coordinates from the previous time-slice,
 * and sends the current coordinates to the next time-slice.
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
 * Receives coordinates from the next time-slice,
 * and sends the current coordinates to the previous time-slice.
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
        bosonic_exchange->prepare();
        bosonic_exchange->exteriorSpringForce(spring_force_arr);
        return;
    }

    // If particles are distinguishable, or if the current bead is an interior bead,
    // the force is calculated based on the standard expression for distinguishable particles.
    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double diff_prev = prev_coord(ptcl_idx, axis) - coord(ptcl_idx, axis);
            double diff_next = next_coord(ptcl_idx, axis) - coord(ptcl_idx, axis);

#if MINIM
            if (pbc) {
                applyMinimumImage(diff_prev, size);
                applyMinimumImage(diff_next, size);
            }
#endif

            spring_force_arr(ptcl_idx, axis) = spring_constant * (diff_prev + diff_next);
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
                // Get the vector distance between the two particles.
                // Here "diff" contains just one vector of dimension NDIM.
                dVec diff = getSeparation(ptcl_one, ptcl_two, MINIM);

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

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double diff = prev_coord(ptcl_idx, axis) - coord(ptcl_idx, axis);

#if MINIM
            if (pbc) {
                applyMinimumImage(diff, size);
            }
#endif

            interior_spring_energy += diff * diff;
        }
    }

    interior_spring_energy *= 0.5 * spring_constant;

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
void Simulation::printReport(double wall_time) const {
    if (this_bead != 0)
        return;

    std::ofstream report_file;
    report_file.open(std::format("{}/report.txt", Output::FOLDER_NAME), std::ios::out | std::ios::app);   

    report_file << "---------\nParameters\n---------\n";

    if (bosonic) {
        report_file << formattedReportLine("Statistics", "Bosonic");
        std::string bosonic_alg_name = "Feldman-Hirshberg";

#if OLD_BOSONIC_ALGORITHM
        bosonic_alg_name = "Naive";
#endif

        report_file << formattedReportLine("Bosonic algorithm", bosonic_alg_name);
    } else {
        report_file << formattedReportLine("Statistics", "Boltzmannonic");
    }
    
    report_file << formattedReportLine("Time propagation algorithm", propagator_type);
    report_file << formattedReportLine("Periodic boundary conditions", pbc);
    report_file << formattedReportLine("Dimension", NDIM);
    report_file << formattedReportLine("Seed", params_seed);
    report_file << formattedReportLine("Coordinate initialization method", init_pos_type);
    report_file << formattedReportLine("Number of atoms", natoms);
    report_file << formattedReportLine("Number of beads", nbeads);

    double out_temperature = Units::convertToUser("temperature", "kelvin", temperature);
    report_file << formattedReportLine("Temperature", std::format("{} kelvin", out_temperature));

    double out_sys_size = Units::convertToUser("length", "angstrom", size);
    report_file << formattedReportLine("Linear size of the system", std::format("{} angstroms", out_sys_size));

    double out_mass = Units::convertToUser("mass", "dalton", mass);
    report_file << formattedReportLine("Mass", std::format("{} amu", out_mass));

    report_file << formattedReportLine("Total number of MD steps", steps);
    report_file << formattedReportLine("Interaction potential name", interaction_potential_name);
    report_file << formattedReportLine("External potential name", external_potential_name);

    report_file << "---------\nFeatures\n---------\n";
    report_file << formattedReportLine("Minimum image convention", MINIM);
    report_file << formattedReportLine("Wrapping of coordinates", WRAP);
    report_file << formattedReportLine("Using i-Pi convention", IPI_CONVENTION);

    report_file << "---------\n";
    report_file << formattedReportLine("Wall time", std::format("{:%T}",
        std::chrono::duration<double>(wall_time)
    ));
    report_file << formattedReportLine("Wall time per step (sec)", std::format("{:.5e}", wall_time / steps));

    report_file.close();
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
        loadMomenta(std::format("init/vel_{:02}.dat", this_bead + 1), mass, momentum_arr); // LAMMPS convention
        //loadMomenta(std::format("init/vel_{:02}.dat", this_bead), mass, momenta); // i-Pi convention

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
        if (units == "on" || units == "true" || units == "none") {
            // In the case of "none", it is expected that the State object quantities will not have units,
            // and therefore providing "atomic_unit" as the units in this case is simply a placeholder.
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

    addObservableIfEnabled(sim_params, "gsf", "gsf");
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
