#pragma once

#include <string>
#include <memory>
#include <vector>
#include <fstream>
#include "params.h"
#include "ordered_map.h"
#include "mpi.h"

class Params; // Forward declaration
class Simulation; // Forward declaration
/* -------------------------------- */

class Observable {
public:
    explicit Observable(Params& param_obj, int _freq, const std::string& _out_unit, int this_bead);
    virtual void calculate() = 0;
    virtual ~Observable() = default;

    void initialize(const std::vector<std::string>& labels);
    void resetValues();

    tsl::ordered_map<std::string, double> quantities;

protected:
    int freq;              // Frequency at which the observable is recorded
    std::string out_unit;  // Units of the output quantities
    int this_bead;
};

/* -------------------------------- */

class ObservableFactory {
public:
    static std::unique_ptr<Observable> createQuantity(Simulation& sim, Params& param_obj, const std::string& observable_type, int _freq,
        const std::string& _out_unit, int this_bead);
};

/* -------------------------------- */

class ObservablesLogger {
public:
    ObservablesLogger(const std::string& filename, int _bead, const std::vector<std::unique_ptr<Observable>>& _observables);
    ~ObservablesLogger();

    void log(int step, MPI_Comm& walker_comm);

private:
    std::ofstream file;
    int bead;
    const std::vector<std::unique_ptr<Observable>>& observables;
};