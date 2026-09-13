#include "dumps/force_dump.h"
#include "core/system_state.h"
#include "output_paths.h"

ForceDump::ForceDump(const std::shared_ptr<SystemState>& state, int this_bead, int out_freq, const std::string& out_unit) :
    Dump(this_bead, out_freq, out_unit), m_state(state)
{
    Units::validateUnit("force", out_unit);
}

#ifdef USE_HDF5
void ForceDump::h5CreateDatasets() {
    const int natoms = m_state->getNumAtoms();
    m_h5_step_ds = H5Utils::make_1d(m_h5file_id, "step", H5T_NATIVE_INT64);
    m_h5_frc_ds = H5Utils::make_frame_ds(m_h5file_id, "forces",
        static_cast<hsize_t>(natoms),
        static_cast<hsize_t>(NDIM));
}
#endif

void ForceDump::output(int step) {
    if (step % m_out_freq != 0)
        return;

    const int natoms = m_state->getNumAtoms();

#ifdef USE_HDF5
    const hsize_t frame = m_h5_frame_count++;

    H5Utils::append_int64(m_h5_step_ds, frame, step);

    std::vector<double> buf(natoms * NDIM);
    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx)
    {
        for (int ax = 0; ax < NDIM; ++ax)
        {
            buf[ptcl_idx * NDIM + ax] = Units::convertToUser(
                "force", 
                m_out_unit,
                m_state->getTotalForce(ptcl_idx, ax)
            );
        }
    }
    H5Utils::append_frame(
        m_h5_frc_ds, 
        frame,
        buf.data(),
        static_cast<hsize_t>(natoms),
        static_cast<hsize_t>(NDIM)
    );
#else
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
                    m_state->getTotalForce(ptcl_idx, axis)
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
