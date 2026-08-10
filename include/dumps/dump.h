#pragma once

#include <string>
#include <fstream>

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

protected:
    int m_this_bead;           // Index of the current imaginary time slice
    int m_out_freq;            // Frequency at which the dump occurs
    std::string m_out_unit;    // Units of the dump quantities
    std::ofstream m_out_file;  // Output file stream

    [[nodiscard]] virtual std::string fileName() const = 0;
};