#pragma once

#include "common.h"
#include "contexts.h"

#include <memory>
#include <functional>
#include <optional>
#include <filesystem>

class RpmdFrameSelector;
enum class XyzFrameSelectionMode;
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
    Simulation(
        int rank, 
        int nproc, 
        const std::shared_ptr<SimulationConfig>& config
    );

    ~Simulation();

    /**
     * @brief Perform a molecular dynamics run using the OBABO scheme.
     */
    void run();
private:
    long m_steps;
    long m_threshold;
    long m_sfreq;
    bool m_is_thermalization_phase;

    // Flag indicating whether the simulation is running in RPMD mode
    bool m_is_rpmd;

    // RPMD frame selector, used to determine which frame to use for each run (only initialized if RPMD is enabled)
    std::unique_ptr<RpmdFrameSelector> m_rpmd_frame_selector;

    /**
     * @brief Initializes the basic contexts that depend only on the configuration parameters
     * and not on the state of the simulation.
     *
     * @param config The simulation configuration object containing parameters like temperature, box size, etc.
     */
    void initializeConfigurationDependentContexts(
        const std::shared_ptr<SimulationConfig>& config
    );

    /**
     * Sets the current simulation step.
     *
     * @param step The step to set.
     */
    void setStep(long step);

    static double getWallTime();

    void runPathIntegralMolecularDynamics();

    void runRingPolymerMolecularDynamics();

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

    /**
     * @brief Ensures the dump objects are initialized with the given filename.
     *
     * @param base_folder The path to the folder where the dumps will be saved.
     */
    void initializeDumps(const std::filesystem::path& base_folder);

    //void calculateObservables() const;

    void calculateAndLogObservables(long step) const;

    /**
     * Calculates the elapsed simulation time and prints a report.
     *
     * @param start_time The time when the simulation started.
     */
    void finalizeSimulation(double start_time) const;

    void executeStep(long step) const;

    // ============ Initialization Factories ============
    /**
     * @brief Type alias for position initialization factories.
     * The factory takes an optional (xyz) frame index and performs the initialization.
     */
    using PositionInitFactory = std::function<void(std::optional<long>)>;

    /**
     * @brief Type alias for momentum initialization factories.
     * The factory takes an optional (xyz) frame index and performs the initialization.
     */
    using MomentumInitFactory = std::function<void(std::optional<long>)>;

    /**
     * @brief Stores the factory for position initialization.
     * Bound during construction of the Simulation object based on the configuration parameters.
     * Invoked with optional frame index during RPMD runs.
     */
    PositionInitFactory m_position_init_factory;

    /**
     * @brief Stores the factory for momentum initialization.
     * Bound during construction of the Simulation object based on the configuration parameters.
     * Invoked with optional frame index during RPMD runs.
     */
    MomentumInitFactory m_momentum_init_factory;

    /**
     * @brief Creates and binds the position initialization factory based on the configuration parameters.
     * This captures all configuration details needed for position initialization without storing them in the Simulation object.
     *
     * @param config The simulation configuration object
     */
    void bindPositionInitFactory(
        const std::shared_ptr<SimulationConfig>& config
    );

    /**
     * @brief Creates and binds the momentum initialization factory based on the configuration parameters.
     * This captures all configuration details needed for momentum initialization without storing them in the Simulation object.
     *
     * @param config The simulation configuration object
    */
    void bindMomentumInitFactory(
        const std::shared_ptr<SimulationConfig>& config
    );

    //std::shared_ptr<SimulationConfig> m_config;
    std::shared_ptr<SystemState> m_state;
    std::shared_ptr<RandomGenerators> m_rng;
    std::shared_ptr<ForceManager> m_force_mgr;
    std::shared_ptr<NormalModes> m_normal_modes;
    std::shared_ptr<Propagator> m_propagator;
    std::shared_ptr<Thermostat> m_thermostat;

    RpmdContext m_rpmd_context;
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

    /*std::vector<std::shared_ptr<Dump>> initializeDumps(
        const std::shared_ptr<SimulationConfig>& config,
        const VelocityContext& vel_ctx
    ) const;*/

     /**
      * Initializes the quantum statistics for the simulation (bosonic or distinguishable).
      *
      * @param is_bosonic Boolean indicating whether the simulation is bosonic or not.
      * @param exchange_xi Exchange parameter for fictitious particles ("xi").
     */
    void initializeQuantumStatistics(bool is_bosonic, double exchange_xi) const;

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

    /*
     * Initializes the forces acting on the particles based on the input parameters.
     *
     * @param state System state object containing information about the current state of the simulation.
     * @param force_mgr Force field manager object containing information about the forces acting on the system.
     * @param spring_ctx Spring context object containing information about the springs.
     * @param box_ctx Box context object containing information about the simulation box.
     */
    void initializeForces(
        const std::shared_ptr<SystemState>& state,
        const std::shared_ptr<ForceManager>& force_mgr, 
        const SpringContext& spring_ctx, 
        const BoxContext& box_ctx
    );

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
