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
class BosonicExchangeBase;
class RandomGenerators;
class ForceManager;
class Propagator;
class Thermostat;
class NormalModes;
class Observable;
class Dump;

class Simulation
{
public:
    Simulation(int rank, int nproc, const std::string& config_filename);
    ~Simulation();

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

    std::shared_ptr<SimulationConfig> m_config;
    std::shared_ptr<SystemState> m_state;
    std::shared_ptr<BosonicExchangeBase> m_bosonic_exchange;
    std::shared_ptr<RandomGenerators> m_rng;
    std::shared_ptr<ForceManager> m_force_mgr;
    std::shared_ptr<NormalModes> m_normal_modes;
    std::shared_ptr<Propagator> m_propagator;
    std::shared_ptr<Thermostat> m_thermostat;

    std::vector<std::shared_ptr<Observable>> m_observables;
    std::vector<std::shared_ptr<Dump>> m_dumps;

    std::vector<std::shared_ptr<Dump>> initializeDumps(const VelocityContext& vel_ctx) const;

     /**
      * Initializes the bosonic exchange machinery based on the input parameters.
      *
      * @param bosonic Boolean indicating whether the simulation is bosonic or not.
      * @param thermal_ctx Thermal context object containing information about the thermal properties of the system.
      * @param spring_ctx Spring context object containing information about the springs in the system.
      * @param box_ctx Box context object containing information about the simulation box.
      * @param bead_ctx Bead context object containing information about the beads in the system.
      * @param state System state object containing information about the current state of the simulation.
      * @return A shared pointer to the initialized exchange state object.
     */
    static std::shared_ptr<BosonicExchangeBase> initializeExchange(
        bool bosonic,
        const ThermalContext& thermal_ctx,
        const SpringContext& spring_ctx,
        const BoxContext& box_ctx,
        const BeadContext& bead_ctx,
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
     * @param bosonic_exchange Bosonic exchange object.
     * @return A shared pointer to the initialized propagator object.
    */
    static std::shared_ptr<Propagator> initializePropagator(
        const std::string& propagator_type,
        double mass,
        double dt,
        const SpringContext& spring_ctx,
        const BoxContext& box_ctx,
        const BeadContext& bead_ctx,
        bool bosonic,
        const std::shared_ptr<SystemState>& state,
        const std::shared_ptr<NormalModes>& normal_modes,
        const std::shared_ptr<ForceManager>& force_mgr,
        const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange
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
    static std::shared_ptr<Thermostat> initializeThermostat(
        const ThermalContext& thermal_ctx,
        const NormalModesContext& nm_ctx,
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
     * @param box_ctx Box context object containing information about the simulation box.
     * @param state System state object.
     * @param rng Random number generator object for generating random numbers.
     */
    static void initializePositions(
        const std::shared_ptr<SimulationConfig>& config,
        const BoxContext& box_ctx,
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
    /*static std::vector<std::shared_ptr<Dump>> initializeDumps(
        const std::shared_ptr<SimulationConfig>& config,
        const std::shared_ptr<SystemState>& state
    );*/

    /**
     * @brief Prints a summary of the simulation parameters at the end of the simulation.
     */
    void printReport(double wall_time) const;
    ///void printDebug(const std::string& text, int target_bead = 0) const;
};
