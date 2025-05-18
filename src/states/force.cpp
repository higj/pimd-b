#include "states/force.h"
#include "units.h"

#include <sstream>

/**
 * @brief Force state class constructor.
 */
ForceState::ForceState(Params& param_obj, int _freq, const std::string& _out_unit, dVec& _forces) :
    State(param_obj, _freq, _out_unit), forces(_forces) {
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
void ForceState::initialize(int this_bead, std::string& output_dir) {
    out_file.open(std::format("{}/force_{}.dat", output_dir, this_bead), std::ios::out | std::ios::app);
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

    out_file << std::format("{}\n", natoms);
    out_file << std::format("Step {}\n", step);

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        out_file << (ptcl_idx + 1) << " 1";

        for (int axis = 0; axis < NDIM; ++axis) {
            out_file << std::format(" {:^20.12e}", Units::convertToUser("force", out_unit, forces(ptcl_idx, axis)));
        }
#if NDIM == 1
        out_file << " 0.0 0.0";
#elif NDIM == 2
        out_file << " 0.0";
#endif
        out_file << "\n";
    }
}
