#ifdef USE_HDF5

#include "initializers/h5_data_loader.h"
#include "units.h"

#include <format>
#include <hdf5.h>
#include <stdexcept>
#include <vector>

// RAII guard so HDF5 handles are always closed on exception.
namespace {
    struct H5File {
        hid_t id;
        explicit H5File(const std::string& path) {
            id = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            if (id < 0)
                throw std::runtime_error(
                    std::format("H5DataLoader: cannot open file '{}'", path));
        }
        ~H5File() { if (id >= 0) H5Fclose(id); }
    };

    struct H5Dataset {
        hid_t id;
        H5Dataset(hid_t file, const std::string& name) {
            id = H5Dopen2(file, name.c_str(), H5P_DEFAULT);
            if (id < 0)
                throw std::runtime_error(
                    std::format("H5DataLoader: cannot open dataset '{}'", name));
        }
        ~H5Dataset() { if (id >= 0) H5Dclose(id); }
    };

    struct H5Space {
        hid_t id;
        explicit H5Space(hid_t ds) { id = H5Dget_space(ds); }
        explicit H5Space(int rank, const hsize_t* dims) {
            id = H5Screate_simple(rank, dims, nullptr);
        }
        ~H5Space() { if (id >= 0) H5Sclose(id); }
    };

    // Resolve a step number to an ordinal frame index by scanning the "step" dataset.
    long resolveStep(hid_t file_id, long step_target, const std::string& h5_filename) {
        H5Dataset step_ds(file_id, "step");
        H5Space step_space(step_ds.id);

        hsize_t n_steps = 0;
        H5Sget_simple_extent_dims(step_space.id, &n_steps, nullptr);

        std::vector<int64_t> steps(n_steps);
        H5Dread(step_ds.id, H5T_NATIVE_INT64,
            H5S_ALL, H5S_ALL, H5P_DEFAULT, steps.data());

        for (hsize_t i = 0; i < n_steps; ++i) {
            if (steps[i] == static_cast<int64_t>(step_target))
                return static_cast<long>(i);
        }
        throw std::runtime_error(std::format(
            "H5DataLoader: file '{}' has no frame with step {}",
            h5_filename, step_target));
    }
}

void H5DataLoader::loadFromFile(
    const std::string& h5_filename,
    const std::string& dataset_name,
    const std::string& data_unit,
    const std::string& unit_family,
    long init_frame,
    FrameSelectionMode init_frame_mode,
    VecArray& destination,
    double prefactor
) {
    H5File file(h5_filename);

    // Resolve Step mode by scanning the companion "step" dataset.
    long frame_idx = (init_frame_mode == FrameSelectionMode::Step)
        ? resolveStep(file.id, init_frame, h5_filename)
        : init_frame;

    H5Dataset ds(file.id, dataset_name);
    H5Space fspace(ds.id);

    // Expect a 3D dataset [n_frames, n_atoms, n_coords].
    int ndims = H5Sget_simple_extent_ndims(fspace.id);
    if (ndims != 3) {
        throw std::runtime_error(std::format(
            "H5DataLoader: dataset '{}' in '{}' has {} dimensions; expected 3",
            dataset_name, h5_filename, ndims));
    }

    hsize_t dims[3];
    H5Sget_simple_extent_dims(fspace.id, dims, nullptr);
    // dims[0] = n_frames, dims[1] = n_atoms, dims[2] = n_coords (always 3 in XYZ convention)

    if (frame_idx < 0 || static_cast<hsize_t>(frame_idx) >= dims[0]) {
        throw std::runtime_error(std::format(
            "H5DataLoader: frame index {} is out of range [0, {}) in '{}'",
            frame_idx, dims[0], h5_filename));
    }

    if (dims[2] < static_cast<hsize_t>(NDIM)) {
        throw std::runtime_error(std::format(
            "H5DataLoader: dataset '{}' in '{}' has {} coordinate columns; "
            "at least {} (NDIM) required",
            dataset_name, h5_filename, dims[2], NDIM));
    }

    const int n_atoms = destination.len();
    if (static_cast<hsize_t>(n_atoms) != dims[1]) {
        throw std::runtime_error(std::format(
            "H5DataLoader: file '{}' frame {} has {} atoms; simulation expects {}",
            h5_filename, frame_idx, dims[1], n_atoms));
    }

    // Read one frame: [frame_idx, :, :] -> flat buffer of size n_atoms * dims[2].
    const hsize_t n_cols = dims[2];
    std::vector<double> buffer(static_cast<std::size_t>(n_atoms) * n_cols);

    hsize_t offset[3] = { static_cast<hsize_t>(frame_idx), 0, 0 };
    hsize_t count[3] = { 1, dims[1], dims[2] };
    H5Sselect_hyperslab(fspace.id, H5S_SELECT_SET, offset, nullptr, count, nullptr);

    H5Space mspace(3, count);
    H5Dread(ds.id, H5T_NATIVE_DOUBLE, mspace.id, fspace.id, H5P_DEFAULT, buffer.data());

    // Apply unit conversion and write into destination.
    for (int atom = 0; atom < n_atoms; ++atom) {
        for (int axis = 0; axis < NDIM; ++axis) {
            double val = buffer[static_cast<std::size_t>(atom) * n_cols + axis];
            double converted = Units::convertToInternal(unit_family, data_unit, val);
            if (unit_family == "velocity")
                converted *= prefactor;
            destination(atom, axis) = converted;
        }
    }
}

long H5DataLoader::countFrames(
    const std::string& h5_filename,
    const std::string& dataset_name
) {
    H5File file(h5_filename);
    H5Dataset ds(file.id, dataset_name);
    H5Space space(ds.id);

    int ndims = H5Sget_simple_extent_ndims(space.id);
    if (ndims < 1) {
        throw std::runtime_error(std::format(
            "H5DataLoader: dataset '{}' in '{}' has no dimensions",
            dataset_name, h5_filename));
    }

    std::vector<hsize_t> dims(ndims);
    H5Sget_simple_extent_dims(space.id, dims.data(), nullptr);
    return static_cast<long>(dims[0]);
}

#endif // USE_HDF5