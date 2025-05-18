#pragma once

#include <string>
#include <fstream>

class Dump {
public:
    /**
     * @brief Generic state class constructor
     */
    explicit Dump(int out_freq, const std::string& out_unit);

    /**
     * @brief Closes the file upon destruction.
     */
    virtual ~Dump();

    virtual void initialize() = 0;
    virtual void output(int step) = 0;   

protected:
    int m_out_freq;            // Frequency at which the dump occurs
    std::string m_out_unit;    // Units of the dump quantities
    std::ofstream m_out_file;  // Output file stream
};