#pragma once

#include <string>
#include <vector>

#include "ordered_map.h"

class Observable {
public:
    /**
     * @brief Generic observable class constructor
     */
    explicit Observable(const std::string& name, int out_freq, const std::string& out_unit);

    virtual void calculate() = 0;
    virtual ~Observable() = default;

    /**
     * @brief Initializes the output folder for the observable.
     *        By default, all observables are recorded to the same file
     *        in the output folder. However, it might be useful to leave
     *        the option to create different folders/files for different observables.
     */
    static void initializeFolder(const std::string& folder_name);

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

    tsl::ordered_map<std::string, double> quantities;

protected:
    std::string m_name;      // Observable name
    int m_out_freq;          // Frequency at which the observable is recorded
    std::string m_out_unit;  // Units of the output quantities
};