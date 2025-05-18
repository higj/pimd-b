#include "dumps/force_dump.h"
#include "core/system_state.h"

ForceDump::ForceDump(const ForceDumpContext& dump_context, int out_freq, const std::string& out_unit) :
    Dump(out_freq, out_unit), m_context(dump_context)
{
    // Check the correctness of the provided unit ahead of time
    try {
        Units::convertToUser("force", out_unit, 1.0);
    } catch (const std::invalid_argument&) {
        throw std::invalid_argument("Invalid output unit for force dump.");
    }
}

void ForceDump::initialize() {
    
    m_out_file.open(std::format("{}/force_{}.dat", Output::FOLDER_NAME, m_context.state->currentBead()), std::ios::out | std::ios::app);
    //m_out_file << std::format("# Units: {}\n", m_out_unit);
}

void ForceDump::output(int step) {
    if (step % m_out_freq != 0)
        return;

    const int natoms = m_context.state->getNumAtoms();

    m_out_file << std::format("{}\n", natoms);
    m_out_file << std::format("Step {}\n", step);

    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        m_out_file << (ptcl_idx + 1) << " 1";

        for (int axis = 0; axis < NDIM; ++axis) {
            m_out_file << std::format(
                " {:^20.12e}", 
                Units::convertToUser(
                    "force", 
                    m_out_unit, 
                    m_context.state->getTotalForce(ptcl_idx, axis)
                )
            );
        }
#if NDIM == 1
        m_out_file << " 0.0 0.0";
#elif NDIM == 2
        m_out_file << " 0.0";
#endif
        m_out_file << "\n";
    }
}
