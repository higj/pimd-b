#pragma once

#include <string>
#include <fstream>
#include <filesystem>

#ifdef USE_HDF5
#include "hdf5_utils.h"
#endif

class Dump {
public:
    /**
     * @brief Generic state class constructor
     */
    explicit Dump(int this_bead, int out_freq, const std::string& out_unit);

    /**
     * @brief Closes the file upon destruction.
     */
    virtual ~Dump();

    //virtual void initialize() = 0;
    virtual void output(int step) = 0;

    /**
     * @brief Close the current output file and reopen a new one with the specified filename.
     * This is useful for RPMD mode where each run needs a separate output file.
     *
     * @param folder The name of the folder where the new dump file will be created. The actual filename will be determined by the derived class.
     */
    virtual void reopenFile(const std::filesystem::path& folder);

    void setMultiRun(const bool is_multi) { m_is_multi_run = is_multi; }

protected:
    int m_this_bead;           // Index of the current imaginary time slice
    int m_out_freq;            // Frequency at which the dump occurs
    std::string m_out_unit;    // Units of the dump quantities
    std::ofstream m_out_file;  // Output file stream

    bool m_is_multi_run = false; // Flag indicating if this is a multi-run scenario
    int m_run_idx = 0;           // Index of the current run (for multi-run scenarios)

    [[nodiscard]] virtual std::string fileName() const = 0;

#ifdef USE_HDF5
    hid_t   m_h5file_id = H5I_INVALID_HID;
    hsize_t m_h5_frame_count = 0;
#ifdef SINGLE_RPMD_FILE
    hid_t   m_h5_run_ds = H5I_INVALID_HID;
#endif

    /// Derived class must create its datasets inside the already-open m_h5file_id.
    virtual void h5CreateDatasets() = 0;

    /// fileName() with the extension replaced by ".h5".
    [[nodiscard]] std::string h5FileName() const;
#endif
};