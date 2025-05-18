#pragma once

#include "common.h"
#include "core/simulation_resources.h"

#include <string>

class Simulation
{
public:
    Simulation(const int& rank, const int& nproc, const std::string& config_filename);
    ~Simulation();

    /**
     * Initializes the simulation box container.
     *
     * @param config Simulation configuration object.
     * @param state System state object.
     * @param rng Random number generator object.
     * @return A shared pointer to the initialized box object.
     */
    /*static std::shared_ptr<Box> initializeBox(
        const std::shared_ptr<SimulationConfig>& config,
        const std::shared_ptr<SystemState>& state,
        const std::shared_ptr<RandomGenerators>& rng);*/

    /**
     * Initializes the bosonic exchange machinery based on the input parameters.
     *
     * @param config Simulation configuration object containing information about bosonic exchange.
     * @param state System state object containing information about the current state of the simulation.
     * @return A shared pointer to the initialized exchange state object.
    */
    static std::shared_ptr<ExchangeState> initializeExchangeState(
        const std::shared_ptr<SimulationConfig>& config, 
        const std::shared_ptr<SystemState>& state
    );

    /**
     * Initializes the normal modes based on the input parameters.
     *
     * @param config Simulation configuration object containing information about normal modes.
     * @param state System state object containing information about the current state of the simulation.
     * @return A shared pointer to the initialized normal modes object.
     */
    static std::shared_ptr<NormalModes> initializeNormalModes(
        const std::shared_ptr<SimulationConfig>& config,
        const std::shared_ptr<SystemState>& state
    );

    /**
     * Initializes the propagator based on the input parameters.
     *
     * @param config Simulation configuration object containing information about the propagator.
     * @param state System state object containing information about the current state of the simulation.
     * @param normal_modes Normal modes object containing information about the normal modes of the system.
     * @param force_mgr Force field manager object containing information about the forces acting on the system.
     * @param exchange_state Exchange state object containing information about the bosonic exchange.
     * @return A shared pointer to the initialized propagator object.
    */
    static std::shared_ptr<Propagator> initializePropagator(
        const std::shared_ptr<SimulationConfig>& config,
        const std::shared_ptr<SystemState>& state,
        const std::shared_ptr<NormalModes>& normal_modes,
        const std::shared_ptr<ForceFieldManager>& force_mgr, 
        const std::shared_ptr<ExchangeState>& exchange_state
    );
    
    /**
     * Initializes the thermostat based on the input parameters.
     *
     * @param config Simulation configuration object containing information about the thermostat.
     * @param state System state object containing information about the current state of the simulation.
     * @param normal_modes Normal modes object containing information about the normal modes of the system.
     * @param rng Random number generator object for generating random numbers.
     * @return A shared pointer to the initialized thermostat object.
    */
    static std::shared_ptr<Thermostat> initializeThermostat(
        const std::shared_ptr<SimulationConfig>& config, 
        const std::shared_ptr<SystemState>& state,
        const std::shared_ptr<NormalModes>& normal_modes,
        const std::shared_ptr<RandomGenerators>& rng
    );

    /**
     * Initializes the observables based on the input parameters.
     *
     * @param config Simulation configuration object.
     * @param state System state object.
     * @param exchange_state Exchange state object.
     * @param force_mgr Force field manager object.
     * @param thermostat Thermostat object.
     * @return A vector of shared pointers to the initialized observable objects.
     */
    static std::vector<std::shared_ptr<Observable>> initializeObservables(
        const std::shared_ptr<SimulationConfig>& config,
        const std::shared_ptr<SystemState>& state,
        const std::shared_ptr<ExchangeState>& exchange_state,
        const std::shared_ptr<ForceFieldManager>& force_mgr,
        const std::shared_ptr<Thermostat>& thermostat
    );

    /**
     * Initializes the positions of the particles based on the input parameters.
     *
     * @param config Simulation configuration object.
     * @param state System state object.
     * @param rng Random number generator object for generating random numbers.
     */
    static void initializePositions(
        const std::shared_ptr<SimulationConfig>& config,
        const std::shared_ptr<SystemState>& state, 
        const std::shared_ptr<RandomGenerators>& rng
    );

    /**
     * Initializes the momenta of the particles based on the input parameters.
     *
     * @param config Simulation configuration object.
     * @param state System state object.
     * @param rng Random number generator object for generating random numbers.
     */
    static void initializeMomenta(
        const std::shared_ptr<SimulationConfig>& config,
        const std::shared_ptr<SystemState>& state, 
        const std::shared_ptr<RandomGenerators>& rng
    );

    /**
     * Initializes the dumps based on the input parameters.
     *
     * @param config Simulation configuration object.
     * @return A vector of shared pointers to the initialized dump objects.
     */
    static std::vector<std::shared_ptr<Dump>> initializeDumps(
        const std::shared_ptr<SimulationConfig>& config,
        const std::shared_ptr<SystemState>& state
    );

    /**
     * @brief Gets the current simulation step.
     */
    [[nodiscard]] int getStep() const;

    /**
     * Sets the current simulation step.
     *
     * @param step The step to set.
     */
    void setStep(int step);

    /**
     * @brief Perform a molecular dynamics run using the OBABO scheme.
     */
    void run();
private:
    int m_step;

    SimulationResources m_context;

    /**
     * @brief Prints a summary of the simulation parameters at the end of the simulation.
     */
    void printReport(double wall_time) const;
    ///void printDebug(const std::string& text, int target_bead = 0) const;
};
