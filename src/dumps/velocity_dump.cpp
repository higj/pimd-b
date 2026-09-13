#include "dumps/velocity_dump.h"
#include "output_paths.h"

VelocityDump::VelocityDump(const VelocityContext& dump_context, int this_bead, int out_freq, const std::string& out_unit) :
    Dump(this_bead, out_freq, out_unit), m_context(dump_context), m_natoms(dump_context.momenta->len())
{
    Units::validateUnit("velocity", out_unit);
}

#ifdef USE_HDF5
void VelocityDump::h5CreateDatasets() {
    m_h5_step_ds = H5Utils::make_1d(m_h5file_id, "step", H5T_NATIVE_INT64);
    m_h5_vel_ds = H5Utils::make_frame_ds(m_h5file_id, "velocities",
        static_cast<hsize_t>(m_natoms),
        static_cast<hsize_t>(NDIM));
}
#endif

void VelocityDump::output(int step)
{
    if (step % m_out_freq != 0)
        return;

#ifdef USE_HDF5
    const hsize_t frame = m_h5_frame_count++;

    H5Utils::append_int64(m_h5_step_ds, frame, step);

    std::vector<double> buf(m_natoms * NDIM);
    for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            buf[ptcl_idx * NDIM + axis] = Units::convertToUser(
                "velocity", 
                m_out_unit,
                (*m_context.momenta)(ptcl_idx, axis) / m_context.mass
            );
        }
    }

    H5Utils::append_frame(
        m_h5_vel_ds, frame,
        buf.data(),
        static_cast<hsize_t>(m_natoms),
        static_cast<hsize_t>(NDIM)
    );
#else
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
#endif
}
