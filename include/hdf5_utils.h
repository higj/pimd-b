#pragma once
#ifdef USE_HDF5

#include <hdf5.h>
#include <stdexcept>
#include <string>

namespace H5Utils {

    /// Create a 1D extensible dataset (chunked + gzip-compressed).
    inline hid_t make_1d(hid_t loc, const char* name, hid_t dtype, hsize_t chunk = 1024) {
        hsize_t init = 0, maxd = H5S_UNLIMITED;
        hid_t space = H5Screate_simple(1, &init, &maxd);
        hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist, 1, &chunk);
        H5Pset_deflate(plist, 6);
        hid_t ds = H5Dcreate2(loc, name, dtype, space, H5P_DEFAULT, plist, H5P_DEFAULT);
        H5Pclose(plist);
        H5Sclose(space);
        if (ds < 0)
            throw std::runtime_error(std::string("H5Utils: failed to create '") + name + "'");
        return ds;
    }

    /// Create a 3D dataset with shape [0, dim1, dim2], extensible on axis 0.
    /// Intended for trajectory frames: [n_frames, n_atoms, NDIM].
    inline hid_t make_frame_ds(hid_t loc, const char* name,
        hsize_t dim1, hsize_t dim2,
        hsize_t chunk_frames = 32) {
        hsize_t init[3] = { 0,              dim1, dim2 };
        hsize_t maxd[3] = { H5S_UNLIMITED,  dim1, dim2 };
        hid_t space = H5Screate_simple(3, init, maxd);
        hsize_t chunk[3] = { chunk_frames, dim1, dim2 };
        hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist, 3, chunk);
        H5Pset_deflate(plist, 6);
        hid_t ds = H5Dcreate2(loc, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, plist, H5P_DEFAULT);
        H5Pclose(plist);
        H5Sclose(space);
        if (ds < 0)
            throw std::runtime_error(std::string("H5Utils: failed to create '") + name + "'");
        return ds;
    }

    /// Append a single double at explicit index idx (extends dataset by 1 if needed).
    inline void append_double(hid_t ds, hsize_t idx, double val) {
        hsize_t newsize = idx + 1;
        H5Dset_extent(ds, &newsize);
        hid_t fspace = H5Dget_space(ds);
        hsize_t one = 1;
        H5Sselect_hyperslab(fspace, H5S_SELECT_SET, &idx, nullptr, &one, nullptr);
        hid_t mspace = H5Screate_simple(1, &one, nullptr);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, mspace, fspace, H5P_DEFAULT, &val);
        H5Sclose(mspace);
        H5Sclose(fspace);
    }

    /// Append a single int32 at explicit index idx.
    inline void append_int32(hid_t ds, hsize_t idx, int val) {
        hsize_t newsize = idx + 1;
        H5Dset_extent(ds, &newsize);
        hid_t fspace = H5Dget_space(ds);
        hsize_t one = 1;
        H5Sselect_hyperslab(fspace, H5S_SELECT_SET, &idx, nullptr, &one, nullptr);
        hid_t mspace = H5Screate_simple(1, &one, nullptr);
        H5Dwrite(ds, H5T_NATIVE_INT, mspace, fspace, H5P_DEFAULT, &val);
        H5Sclose(mspace);
        H5Sclose(fspace);
    }

    /// Append a single int64 at explicit index idx.
    inline void append_int64(hid_t ds, hsize_t idx, int64_t val) {
        hsize_t newsize = idx + 1;
        H5Dset_extent(ds, &newsize);
        hid_t fspace = H5Dget_space(ds);
        hsize_t one = 1;
        H5Sselect_hyperslab(fspace, H5S_SELECT_SET, &idx, nullptr, &one, nullptr);
        hid_t mspace = H5Screate_simple(1, &one, nullptr);
        H5Dwrite(ds, H5T_NATIVE_INT64, mspace, fspace, H5P_DEFAULT, &val);
        H5Sclose(mspace);
        H5Sclose(fspace);
    }

    /// Append a [dim1 x dim2] frame to a 3D dataset at frame_idx.
    inline void append_frame(hid_t ds, hsize_t frame_idx,
        const double* data, hsize_t dim1, hsize_t dim2) {
        hsize_t newsize[3] = { frame_idx + 1, dim1, dim2 };
        H5Dset_extent(ds, newsize);
        hid_t fspace = H5Dget_space(ds);
        hsize_t offset[3] = { frame_idx, 0, 0 };
        hsize_t count[3] = { 1, dim1, dim2 };
        H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        hid_t mspace = H5Screate_simple(3, count, nullptr);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, mspace, fspace, H5P_DEFAULT, data);
        H5Sclose(mspace);
        H5Sclose(fspace);
    }

} // namespace H5Utils
#endif // USE_HDF5