#include "dumps/dump.h"

Dump::Dump(int this_bead, int out_freq, const std::string& out_unit) : m_this_bead(this_bead), m_out_freq(out_freq),
                                                                       m_out_unit(out_unit)
{
}

Dump::~Dump()
{
#ifdef USE_HDF5
    if (m_h5file_id != H5I_INVALID_HID) {
        H5Fclose(m_h5file_id);
        m_h5file_id = H5I_INVALID_HID;
    }
#else
    if (m_out_file.is_open()) {
        m_out_file.close();
    }
#endif
}

void Dump::reopenFile(const std::filesystem::path& folder) {
    // Ensure the output directory exists
    std::filesystem::create_directories(folder);

#ifdef USE_HDF5
    // Close any previously open HDF5 file
    if (m_h5file_id != H5I_INVALID_HID) {
        H5Fclose(m_h5file_id);
        m_h5file_id = H5I_INVALID_HID;
    }
    m_h5_frame_count = 0;

    const std::filesystem::path full_path = folder / h5FileName();
    m_h5file_id = H5Fcreate(full_path.string().c_str(),
        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (m_h5file_id < 0) {
        throw std::ios_base::failure("Failed to create HDF5 file: " + full_path.string());
    }

    h5CreateDatasets();   // derived class creates its datasets
#else
    // Close the current file
    if (m_out_file.is_open()) {
        m_out_file.close();
    }

    // Open the new file
    const std::filesystem::path full_path = folder / fileName();

    // Open the file in append mode to avoid overwriting existing data
    m_out_file.open(full_path, std::ios::out | std::ios::app);

    if (!m_out_file.is_open()) {
        throw std::ios_base::failure("Failed to open " + full_path.string());
    }
#endif
}

#ifdef USE_HDF5
std::string Dump::h5FileName() const {
    std::string name = fileName();
    const auto dot = name.rfind('.');
    return (dot != std::string::npos ? name.substr(0, dot) : name) + ".h5";
}
// e.g. "position_0.xyz" -> "position_0.h5", "force_0.dat" -> "force_0.h5"
#endif