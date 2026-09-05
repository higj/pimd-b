#include "dumps/dump.h"

Dump::Dump(int this_bead, int out_freq, const std::string& out_unit) : m_this_bead(this_bead), m_out_freq(out_freq),
                                                                       m_out_unit(out_unit)
{
}

Dump::~Dump()
{
    if (m_out_file.is_open())
    {
        m_out_file.close();
    }
}

void Dump::reopenFile(const std::filesystem::path& folder) {
    // Close the current file
    if (m_out_file.is_open()) {
        m_out_file.close();
    }

    // Ensure the output directory exists
    std::filesystem::create_directories(folder);

    // Open the new file
    const std::filesystem::path full_path = folder / fileName();

    // Open the file in append mode to avoid overwriting existing data
    m_out_file.open(full_path, std::ios::out | std::ios::app);

    if (!m_out_file.is_open()) {
        throw std::ios_base::failure("Failed to open " + full_path.string());
    }
}
