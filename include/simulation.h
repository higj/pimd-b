#pragma once

#include "common.h"
#include "contexts/thermal_context.h"
#include "contexts/bead_context.h"
#include "contexts/spring_context.h"
#include "contexts/box_context.h"
#include "contexts/normal_modes_context.h"
#include "contexts/thermostat_context.h"
#include "contexts/velocity_context.h"

#include <string>
#include <memory>

struct SimulationConfig;
class SystemState;
class RandomGenerators;
class ForceManager;
class Propagator;
class Thermostat;
class NormalModes;
class Observable;
class Dump;
class SimulationReport;
class ObservablesLogger;

class Simulation
{
public:
    Simulation(int rank, int nproc, const std::string& config_filename);
    ~Simulation();

    /**
     * @brief Perform a molecular dynamics run using the OBABO scheme.
     */
    void run();
private:
    long m_steps;
    long m_threshold;
    bool m_is_thermalization_phase;

    /**
     * Sets the current simulation step.
     *
     * @param step The step to set.
     */
    void setStep(long step);

    static double getWallTime();

    /**
     * @brief Reset the observables at the beginning of each step.
     * TODO: Do we need this? What if we want accumulation? Does it account for frequency?
     */
    void resetObservables() const;

    /**
     * @brief Zero momentum after every thermostat step (if needed).
     * TODO: Maybe should be part of the thermostat step?
     */
    void zeroMomentumIfRequired() const;

    /**
     * Dump the desired quantities (e.g., coordinates, forces, etc.) at the specified frequency.
     * @param step The current simulation step.
     */
    void dumpStepInfo(long step) const;

    void performMolecularDynamicsStep() const;

    void calculateObservables() const;

    void calculateAndLogObservables(long step) const;

    /**
     * Calculates the elapsed simulation time and prints a report.
     *
     * @param start_time The time when the simulation started.
     */
    void finalizeSimulation(double start_time) const;

    void executeStep(long step) const;

    //std::shared_ptr<SimulationConfig> m_config;
    std::shared_ptr<SystemState> m_state;
    std::shared_ptr<RandomGenerators> m_rng;
    std::shared_ptr<ForceManager> m_force_mgr;
    std::shared_ptr<NormalModes> m_normal_modes;
    std::shared_ptr<Propagator> m_propagator;
    std::shared_ptr<Thermostat> m_thermostat;

    ThermalContext m_thermal_ctx;
    SpringContext m_spring_ctx;
    BoxContext m_box_ctx;
    BeadContext m_bead_ctx;
    NormalModesContext m_nm_ctx;
    ThermostatContext m_thermostat_ctx;
    VelocityContext m_velocity_ctx;

    std::vector<std::shared_ptr<Observable>> m_observables;
    std::vector<std::shared_ptr<Dump>> m_dumps;

    std::unique_ptr<ObservablesLogger> m_obs_logger;

    std::unique_ptr<SimulationReport> m_report;

    std::vector<std::shared_ptr<Dump>> initializeDumps(
        const std::shared_ptr<SimulationConfig>& config,
        const VelocityContext& vel_ctx
    ) const;

     /**
      * Initializes the quantum statistics for the simulation (bosonic or distinguishable).
      *
      * @param is_bosonic Boolean indicating whether the simulation is bosonic or not.
     */
    void initializeQuantumStatistics(bool is_bosonic) const;

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
     * @param propagator_type Type of the propagator to initialize (e.g., "verlet", "normal modes", etc.).
     * @param mass Mass of the particles in the system.
     * @param dt Time step for the simulation.
     * @param spring_ctx Spring context object containing information about the springs.
     * @param box_ctx Box context object containing information about the simulation box.
     * @param bead_ctx Bead context object containing information about the beads in the system.
     * @param bosonic Boolean indicating whether the simulation is bosonic or not.
     * @param state System state object containing information about the current state of the simulation.
     * @param normal_modes Normal modes object containing information about the normal modes of the system.
     * @param force_mgr Force field manager object containing information about the forces acting on the system.
     * @return A shared pointer to the initialized propagator object.
    */
    std::shared_ptr<Propagator> initializePropagator(
        const std::shared_ptr<SimulationConfig>& config,
        const std::shared_ptr<SystemState>& state,
        const std::shared_ptr<NormalModes>& normal_modes,
        const std::shared_ptr<ForceManager>& force_mgr
    );

    /**
     * Initializes the thermostat based on the input parameters.
     *
     * @param thermal_ctx Thermal context object containing information about the thermal properties of the system.
     * @param nm_ctx Normal modes context object containing information about the normal modes of the system.
     * @param config Simulation configuration object containing information about the thermostat.
     * @param state System state object containing information about the current state of the simulation.
     * @param rng Random number generator object for generating random numbers.
     * @return A shared pointer to the initialized thermostat object.
    */
    std::shared_ptr<Thermostat> initializeThermostat(
        const std::shared_ptr<SimulationConfig>& config,
        const std::shared_ptr<SystemState>& state,
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
     /*static std::vector<std::shared_ptr<Observable>> initializeObservables(
        const std::shared_ptr<SimulationConfig>& config,
        const std::shared_ptr<SystemState>& state,
        const std::shared_ptr<ExchangeState>& exchange_state,
        const std::shared_ptr<ForceManager>& force_mgr,
        const std::shared_ptr<Thermostat>& thermostat
    );*/

    /**
     * Initializes the positions of the particles based on the input parameters.
     *
     * @param config Simulation configuration object.
     */
    void initializePositions(const std::shared_ptr<SimulationConfig>& config);

    /**
     * Initializes the momenta of the particles based on the input parameters.
     *
     * @param config Simulation configuration object.
     */
    void initializeMomenta(const std::shared_ptr<SimulationConfig>& config);

    /**
     * Initializes the dumps based on the input parameters.
     *
     * @param config Simulation configuration object.
     * @return A vector of shared pointers to the initialized dump objects.
     */
    /*static std::vector<std::shared_ptr<Dump>> initializeDumps(
        const std::shared_ptr<SimulationConfig>& config,
        const std::shared_ptr<SystemState>& state
    );*/

    ///void printDebug(const std::string& text, int target_bead = 0) const;
};
