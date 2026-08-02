#include "dumps/velocity_dump.h"
#include "output_paths.h"

VelocityDump::VelocityDump(const VelocityContext& dump_context, int this_bead, int out_freq, const std::string& out_unit) :
    Dump(this_bead, out_freq, out_unit), m_context(dump_context), m_natoms(dump_context.momenta->len())
{
    // Check the correctness of the provided unit ahead of time
    try
    {
        Units::convertToUser("velocity", out_unit, 1.0);
    }
    catch (const std::invalid_argument&)
    {
        throw std::invalid_argument("Invalid output unit for velocity dump.");
    }
}

void VelocityDump::initialize()
{
    // Open the output file for appending, creating it if it doesn't exist
    // TODO: Consider using std::filesystem to ensure the output directory exists before opening the file
    // TODO: Note that velocities are dumped to xyz files, which will likely break the tests that expect the old (.dat) format.
    // We should consider changing the file extension to .vel or similar.
    m_out_file.open(std::format("{}/velocity_{}.xyz", Output::FOLDER_NAME, m_this_bead), std::ios::out | std::ios::app);
    //m_out_file << std::format("# Units: {}\n", m_out_unit);
}

void VelocityDump::output(int step)
{
    if (step % m_out_freq != 0)
        return;

    m_out_file << std::format("{}\n", m_natoms);
    m_out_file << std::format("Step {}\n", step);

    for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx)
    {
        //m_out_file << (ptcl_idx + 1) << " 1"; // Old format: particle index and type

        // We want the velocity dump to be in the same format as the position dump,
        // so we use "1" as a placeholder for the atom type. (TODO: Might break tests that expect the old format)
        m_out_file << "1";

        for (int axis = 0; axis < NDIM; ++axis)
        {
            m_out_file << std::format(
                " {:^20.12e}",
                Units::convertToUser(
                    "velocity",
                    m_out_unit,
                    (*m_context.momenta)(ptcl_idx, axis) / m_context.mass
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
