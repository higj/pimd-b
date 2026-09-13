#pragma once
#ifdef USE_HDF5

#include "core/simulation_config.h"  // VecArray, FrameSelectionMode

#include <string>

/**
 * @brief Utility for loading data from HDF5 trajectory files with unit conversion.
 *
 * Reads one frame from a 3D dataset with shape [n_frames, n_atoms, n_coords].
 * This is the HDF5 counterpart of XyzDataLoader; the interface is intentionally
 * parallel so the two paths are interchangeable from the initializer level.
 *
 * Frame selection mirrors XyzDataLoader:
 *   - Index mode  : the ordinal position in the dataset (used by RPMD frame selector)
 *   - Step mode   : the value in the companion "step" dataset that matches init_frame
 */
class H5DataLoader {
public:
    /**
     * @brief Loads one frame from an HDF5 trajectory dataset.
     *
     * @param h5_filename   Path to the HDF5 file
     * @param dataset_name  Name of the 3D dataset ("positions" or "velocities")
     * @param data_unit     Unit declared in config (e.g. "angstrom", "angstrom/fs")
     * @param unit_family   "length" or "velocity"
     * @param init_frame    Frame index (Index mode) or step number (Step mode)
     * @param init_frame_mode  How to interpret init_frame
     * @param destination   Target VecArray to write converted data into
     * @param prefactor     Multiplied after unit conversion (mass, for velocity->momenta)
     */
    static void loadFromFile(
        const std::string& h5_filename,
        const std::string& dataset_name,
        const std::string& data_unit,
        const std::string& unit_family,
        long init_frame,
        FrameSelectionMode init_frame_mode,
        VecArray& destination,
        double prefactor = 1.0
    );

    /**
     * @brief Returns the number of frames in a dataset (extent of axis 0).
     *
     * @param h5_filename   Path to the HDF5 file
     * @param dataset_name  Name of the dataset to query
     * @return Total frame count
     */
    static long countFrames(
        const std::string& h5_filename,
        const std::string& dataset_name
    );
};

#endif // USE_HDF5