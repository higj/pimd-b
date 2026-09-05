#include "dumps/position_dump.h"
#include "output_paths.h"

PositionDump::PositionDump(const std::shared_ptr<const VecArray>& coord, int this_bead, int out_freq,
                           const std::string& out_unit) :
    Dump(this_bead, out_freq, out_unit), m_coord(coord), m_natoms(coord->len())
{
    Units::validateUnit("length", out_unit);
}

/*
void PositionDump::initialize()
{
    m_out_file.open(std::format("{}/position_{}.xyz", Output::FOLDER_NAME, m_this_bead), std::ios::out | std::ios::app);
    //m_out_file << std::format("# Units: {}\n", m_out_unit);
}
*/

void PositionDump::output(int step)
{
    if (step % m_out_freq != 0)
        return;

    m_out_file << std::format("{}\n", m_natoms);
    //m_out_file << std::format(" Atoms. MD step: {}\n", step);
    m_out_file << std::format("Step {}\n", step);

    for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx)
    {
        m_out_file << "1";

        for (int axis = 0; axis < NDIM; ++axis)
        {
            m_out_file << std::format(
                " {:^20.12e}",
                Units::convertToUser(
                    "length",
                    m_out_unit,
                    (*m_coord)(ptcl_idx, axis)
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
