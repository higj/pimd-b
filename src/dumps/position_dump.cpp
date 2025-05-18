#include "dumps/position_dump.h"

PositionDump::PositionDump(const PositionDumpContext& dump_context, int out_freq, const std::string& out_unit) :
    Dump(out_freq, out_unit), m_context(dump_context) {
    // Check the correctness of the provided unit ahead of time
    try {
        Units::convertToUser("length", out_unit, 1.0);
    } catch (const std::invalid_argument&) {
        throw std::invalid_argument("Invalid output unit for position dump.");
    }
}

void PositionDump::initialize() {
    m_out_file.open(std::format("{}/position_{}.xyz", Output::FOLDER_NAME, m_context.this_bead), std::ios::out | std::ios::app);
    //m_out_file << std::format("# Units: {}\n", m_out_unit);
}

void PositionDump::output(int step) {
    if (step % m_out_freq != 0)
        return;

    m_out_file << std::format("{}\n", m_context.natoms);
    //m_out_file << std::format(" Atoms. MD step: {}\n", step);
    m_out_file << std::format("Step {}\n", step);

    for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
        m_out_file << "1";

        for (int axis = 0; axis < NDIM; ++axis) {
            m_out_file << std::format(
                " {:^20.12e}",
                Units::convertToUser(
                    "length", 
                    m_out_unit, 
                    (*m_context.coord)(ptcl_idx, axis)
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