#pragma once

#include <string>
#include <memory>
#include <vector>
#include <fstream>

#include "ordered_map.h"

class Simulation; // Forward declaration

/* -------------------------------- */

class Observable {
public:
    explicit Observable(const Simulation& _sim, int _freq, const std::string& _out_unit);
    virtual void calculate() = 0;
    virtual ~Observable() = default;

    void initialize(const std::vector<std::string>& labels);
    void resetValues();

    tsl::ordered_map<std::string, double> quantities;

protected:
    const Simulation& sim; // Reference to the simulation object
    int freq;              // Frequency at which the observable is recorded
    std::string out_unit;  // Units of the output quantities
};

/* -------------------------------- */

class ObservableFactory {
public:
    static std::unique_ptr<Observable> createQuantity(const std::string& observable_type, const Simulation& _sim, int _freq,
        const std::string& _out_unit);
};

/* -------------------------------- */

class ObservablesLogger {
public:
    ObservablesLogger(const std::string& filename, int _bead, const std::vector<std::unique_ptr<Observable>>& _observables);
    ~ObservablesLogger();

    void log(int step);

private:
    std::ofstream file;
    int bead;
    const std::vector<std::unique_ptr<Observable>>& observables;
};

/* -------------------------------- */

class EnergyObservable : public Observable {
public:
    EnergyObservable(const Simulation& _sim, int _freq, const std::string& _out_unit);

    void calculate() override;

private:
    void calculateKinetic();
    void calculatePotential();
};

/* -------------------------------- */

class ClassicalObservable : public Observable {
public:
    ClassicalObservable(const Simulation& _sim, int _freq, const std::string& _out_unit);

    void calculate() override;

private:
    void calculateKineticEnergy();
    void calculateSpringEnergy();
};


/* -------------------------------- */

class BosonicObservable : public Observable {
public:
    BosonicObservable(const Simulation& _sim, int _freq, const std::string& _out_unit);

    void calculate() override;
};
