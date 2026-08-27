#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "ordered_map.h"

class Observable {
public:
    /**
     * @brief Generic observable class constructor
     */
    explicit Observable(
        const std::string& name, 
        int out_freq, 
        const std::string& out_unit,
        const std::string& out_filename = ""
    );

    virtual void calculate() = 0;
    virtual ~Observable() = default;

    /**
     * @brief Initializes the output folder for the observable.
     *        By default, all observables are recorded to the same file
     *        in the output folder. However, it might be useful to leave
     *        the option to create different folders/files for different observables.
     */
    static void initializeFolder(std::string_view folder_name);

    /**
     * Initializes observable with the given label.
     *
     * @param label Labels of the quantity to be calculated
     */
    void initializeLabel(const std::string& label);

    /**
     * Initializes observables with the given labels.
     *
     * @param labels Labels of the quantities to be calculated
     */
    void initialize(const std::vector<std::string>& labels);

    /**
     * @brief Resets the values of the observable to zero.
     * Useful for clearing results from the previous molecular dynamics step. /// TODO: Is it actually useful?
     */
    void resetValues();

    /**
     * @brief Get the name of the observable.
     * @return String representing the name of the observable.
     */
    [[nodiscard]] std::string name() const { return m_name; }

    /**
     * @brief Get the output filename for this observable.
     * @return The output filename, relative to the active simulation output directory.
     *         Empty string means use the default main file.
     */
    [[nodiscard]] std::string outputFilename() const { return m_out_filename; }

    /**
     * @brief Select the output filename for this observable.
     * @param filename A filename relative to the active simulation output directory.
     *        An empty filename uses the standard Output::MAIN_FILENAME.
     */
    void setOutputFilename(const std::string_view filename) { m_out_filename = filename; }

    /**
     * @brief Check if this observable uses a custom output file.
     * @return true if a custom output file was specified, false if using default.
     */
    [[nodiscard]] bool usesCustomFile() const { return !m_out_filename.empty(); }

    tsl::ordered_map<std::string, double> quantities;

protected:
    std::string m_name;          // Observable name
    int m_out_freq;              // Frequency at which the observable is recorded
    std::string m_out_unit;      // Units of the output quantities
    std::string m_out_filename;  // Custom output file (optional, empty = use central logger)
};
