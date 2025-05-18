#pragma once

#include <string>
#include <vector>

#include "ordered_map.h"

class Observable {
public:
    /**
     * @brief Generic observable class constructor
     */
    explicit Observable(int out_freq, const std::string& out_unit);

    virtual void calculate() = 0;
    virtual ~Observable() = default;

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

    tsl::ordered_map<std::string, double> quantities;

protected:
    int m_out_freq;          // Frequency at which the observable is recorded
    std::string m_out_unit;  // Units of the output quantities
};