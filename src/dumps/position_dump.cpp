#include "dumps/position_dump.h"
#include "output_paths.h"

PositionDump::PositionDump(const std::shared_ptr<const VecArray>& coord, int this_bead, int out_freq,
                           const std::string& out_unit) :
    Dump(this_bead, out_freq, out_unit), m_coord(coord), m_natoms(coord->len())
{
    Units::validateUnit("length", out_unit);
}

#ifdef USE_HDF5
void PositionDump::h5CreateDatasets() {
    // step [N],  positions [N, n_atoms, 3]
    // Third dimension is always 3 regardless of NDIM, consistent with the
    // XYZ text format which always writes three coordinates per atom.
    m_h5_step_ds = H5Utils::make_1d(m_h5file_id, "step", H5T_NATIVE_INT64);
    m_h5_pos_ds = H5Utils::make_frame_ds(
        m_h5file_id, 
        "positions",
        static_cast<hsize_t>(m_natoms),
        3
    );
}
#endif

void PositionDump::output(int step) {
    if (step % m_out_freq != 0)
        return;

#ifdef USE_HDF5
    const hsize_t frame = m_h5_frame_count++;

    H5Utils::append_int64(m_h5_step_ds, frame, step);

#ifdef SINGLE_RPMD_FILE
    if (m_is_multi_run) {
        H5Utils::append_int64(m_h5_run_ds, frame, m_run_idx);
    }
#endif

    std::vector<double> buf(m_natoms * 3, 0.0);  // always 3 columns, zero-padded
    for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            buf[ptcl_idx * 3 + axis] = Units::convertToUser(
                "length",
                m_out_unit,
                (*m_coord)(ptcl_idx, axis)
            );
        }
    }
    H5Utils::append_frame(
        m_h5_pos_ds,
        frame,
        buf.data(),
        static_cast<hsize_t>(m_natoms),
        3
    );
#else
    m_out_file << std::format("{}\n", m_natoms);
#ifdef SINGLE_RPMD_FILE
    if (m_is_multi_run) {
        m_out_file << std::format("Step {} Run {}\n", step, m_run_idx);
    } else {
        //m_out_file << std::format(" Atoms. MD step: {}\n", step);
        m_out_file << std::format("Step {}\n", step);
    }
#else
    //m_out_file << std::format(" Atoms. MD step: {}\n", step);
    m_out_file << std::format("Step {}\n", step);
#endif

    for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx) {
        m_out_file << "1";

        for (int axis = 0; axis < NDIM; ++axis) {
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
#endif
}
