#include "states/position.h"
#include "simulation.h"
#include "units.h"

#include <sstream>

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