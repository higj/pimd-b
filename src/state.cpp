#include "state.h"
#include "simulation.h"
#include "units.h"

#include <sstream>

/**
 * @brief Generic state class constructor
 */
State::State(const Simulation& _sim, int _freq, const std::string& _out_unit) :
    sim(_sim), freq(_freq), out_unit(_out_unit) {
}

/**
 * @brief Performs cleanup (e.g. closing the output file).
 */
void State::finalize() {
    if (out_file.is_open()) {
        out_file.close();
    }
}

/**
 * @brief Creates a state object based on the state type.
 *
 * @param state_type Type of state to create
 * @param _sim Simulation object
 * @param _freq Frequency at which the state is recorded
 * @param _out_unit Output unit of the state
 * @return std::unique_ptr<State> Pointer to the created state object
 * @throw std::invalid_argument If the state type is unknown
 */
std::unique_ptr<State> StateFactory::createQuantity(const std::string& state_type,
    const Simulation& _sim,
    int _freq,
    const std::string& _out_unit) {
    if (state_type == "position") {
        return std::make_unique<PositionState>(_sim, _freq, _out_unit);
    } else if (state_type == "velocity") {
        return std::make_unique<VelocityState>(_sim, _freq, _out_unit);
    } else if (state_type == "force") {
        return std::make_unique<ForceState>(_sim, _freq, _out_unit);
    } else if (state_type == "wind_prob") {
        return std::make_unique<WindingProbabilityState>(_sim, _freq, _out_unit);
    } else {
        throw std::invalid_argument("Unknown state type.");
    }
}

/**
 * @brief Position state class constructor.
 */
PositionState::PositionState(const Simulation& _sim, int _freq, const std::string& _out_unit) :
    State(_sim, _freq, _out_unit) {
    // Check the correctness of the provided unit ahead of time
    try {
        Units::convertToUser("length", out_unit, 1.0);
    } catch (const std::invalid_argument&) {
        throw std::invalid_argument("Invalid output unit for position state.");
    }
}

/**
 * @brief Initializes the coordinates xyz file.
 */
void PositionState::initialize() {
    out_file.open(std::format("{}/position_{}.xyz", Output::FOLDER_NAME, sim.this_bead), std::ios::out | std::ios::app);
    //out_file << std::format("# Units: {}\n", out_unit);
}

/**
 * Outputs the trajectories.
 *
 * @param step Current step of the simulation.
 */
void PositionState::output(int step) {
    if (step % freq != 0)
        return;

    out_file << std::format("{}\n", sim.natoms);
    //out_file << std::format(" Atoms. MD step: {}\n", step);
    out_file << std::format("Step {}\n", step);

    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        out_file << "1";

        for (int axis = 0; axis < NDIM; ++axis) {
            out_file << std::format(" {:^20.12e}", Units::convertToUser("length", out_unit, sim.coord(ptcl_idx, axis)));
        }
#if NDIM == 1
        out_file << " 0.0 0.0";
#elif NDIM == 2
        out_file << " 0.0";
#endif
        out_file << "\n";
    }
}

/**
 * @brief Velocity state class constructor.
 */
VelocityState::VelocityState(const Simulation& _sim, int _freq, const std::string& _out_unit) :
    State(_sim, _freq, _out_unit) {
    // Check the correctness of the provided unit ahead of time
    try {
        Units::convertToUser("velocity", out_unit, 1.0);
    } catch (const std::invalid_argument&) {
        throw std::invalid_argument("Invalid output unit for velocity state.");
    }
}

/**
 * @brief Initializes the velocities dat file.
 */
void VelocityState::initialize() {
    out_file.open(std::format("{}/velocity_{}.dat", Output::FOLDER_NAME, sim.this_bead), std::ios::out | std::ios::app);
    //out_file << std::format("# Units: {}\n", out_unit);
}

/**
 * Outputs the velocities.
 *
 * @param step Current step of the simulation.
 */
void VelocityState::output(int step) {
    if (step % freq != 0)
        return;

    out_file << std::format("{}\n", sim.natoms);
    out_file << std::format("Step {}\n", step);

    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        out_file << (ptcl_idx + 1) << " 1";

        for (int axis = 0; axis < NDIM; ++axis) {
            out_file << std::format(" {:^20.12e}", Units::convertToUser("velocity", out_unit, sim.momenta(ptcl_idx, axis) / sim.mass));
        }
#if NDIM == 1
        out_file << " 0.0 0.0";
#elif NDIM == 2
        out_file << " 0.0";
#endif
        out_file << "\n";
    }
}

/**
 * @brief Force state class constructor.
 */
ForceState::ForceState(const Simulation& _sim, int _freq, const std::string& _out_unit) :
    State(_sim, _freq, _out_unit) {
    // Check the correctness of the provided unit ahead of time
    try {
        Units::convertToUser("force", out_unit, 1.0);
    } catch (const std::invalid_argument&) {
        throw std::invalid_argument("Invalid output unit for force state.");
    }
}

/**
 * @brief Initializes the forces dat file.
 */
void ForceState::initialize() {
    out_file.open(std::format("{}/force_{}.dat", Output::FOLDER_NAME, sim.this_bead), std::ios::out | std::ios::app);
    //out_file << std::format("# Units: {}\n", out_unit);
}

/**
 * Outputs the forces.
 *
 * @param step Current step of the simulation.
 */
void ForceState::output(int step) {
    if (step % freq != 0)
        return;

    out_file << std::format("{}\n", sim.natoms);
    out_file << std::format("Step {}\n", step);

    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        out_file << (ptcl_idx + 1) << " 1";

        for (int axis = 0; axis < NDIM; ++axis) {
            out_file << std::format(" {:^20.12e}", Units::convertToUser("force", out_unit, sim.forces(ptcl_idx, axis)));
        }
#if NDIM == 1
        out_file << " 0.0 0.0";
#elif NDIM == 2
        out_file << " 0.0";
#endif
        out_file << "\n";
    }
}

/**
 * @brief Winding probability state class constructor.
 */
WindingProbabilityState::WindingProbabilityState(const Simulation& _sim, int _freq, const std::string& _out_unit) :
    State(_sim, _freq, _out_unit) {
}

/**
 * @brief Initializes the winding probabilities log file.
 */
void WindingProbabilityState::initialize() {
    if (!sim.pbc || !sim.apply_wind)
        return;

    out_file.open(std::format("{}/wind-prob-{}.log", Output::FOLDER_NAME, sim.this_bead), std::ios::out | std::ios::app);
    out_file << std::format("{:^16s}", "step");
    out_file << std::format(" {:^8s}", "ptcl");

    for (int wind_idx = 0; wind_idx < sim.wind.len(); ++wind_idx) {
        std::ostringstream wind_ss;
        wind_ss << "(";
        for (int axis = 0; axis < NDIM; ++axis) {
            wind_ss << sim.wind(wind_idx, axis);
            if (axis != NDIM - 1) {
                wind_ss << ",";
            }
        }
        wind_ss << ")";
        std::string wind_str = wind_ss.str();
        out_file << std::format(" {:^16s}", wind_str);
    }

    out_file << '\n';
}

/**
 * Outputs the winding probabilities.
 *
 * @param step Current step of the simulation.
 */
void WindingProbabilityState::output(int step) {
    if (step % freq != 0 || !sim.pbc || !sim.apply_wind)
        return;

    for (int ptcl_idx = 0; ptcl_idx < sim.natoms; ++ptcl_idx) {
        out_file << std::format("{:^16.8e}", static_cast<double>(step));
        out_file << std::format(" {:^8d} ", ptcl_idx + 1);
        for (int wind_idx = 0; wind_idx < sim.wind.len(); ++wind_idx) {
            double prob = 1.0;

            for (int axis = 0; axis < NDIM; ++axis) {
                prob *= sim.getWindingProbability(sim.coord(ptcl_idx, axis) - sim.next_coord(ptcl_idx, axis), sim.wind(wind_idx, axis));
            }

            out_file << std::format(" {:^16.8e}", prob);
        }
        out_file << '\n';
    }
    out_file << '\n';
}