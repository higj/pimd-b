#pragma once

#include <string>
#include <memory>
#include <fstream>

class Simulation; // Forward declaration

/* -------------------------------- */

class State {
public:
    explicit State(const Simulation& _sim, int _freq, const std::string& _out_unit);
    ~State();

    virtual void initialize() = 0;
    virtual void output(int step) = 0;   

protected:
    const Simulation& sim;   // Reference to the simulation object
    int freq;                // Frequency at which the state is recorded
    std::ofstream out_file;  // Output file stream
    std::string out_unit;    // Units of the output quantities
};

/* -------------------------------- */

class StateFactory {
public:
    static std::unique_ptr<State> createQuantity(const std::string& state_type, 
        const Simulation& _sim, int _freq, const std::string& _out_unit);
};