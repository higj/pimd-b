#pragma once

#include <string>
#include <random>
#include <ctime>
#include <memory>

#include "common.h"
#include "params.h"
#include "potentials.h"
#include "bosonic_exchange.h"
#include "mpi.h"

class NormalModes;
class State;
class Observable;
class Propagator;
class Thermostat;
class Coupling;
class WalkersCommunicationBase;

class Simulation
{
public:
    double temperature;
    double beta;        // Thermodynamic beta 1/(kB*T)
    double thermo_beta; // Inverse temperature at which the simulation is actually performed ("thermostat beta")
    double dt;          // Timestep
    double size;        // Linear system size (TODO: Add support for Ly, Lz,...)
    double gamma;       // Friction constant of the Langevin thermostat
    int    nchains;     // Number of Nose-Hoover chains
    double threshold;   // Percentage of steps to throw away (thermalization)

    int natoms;         // Number of atoms in the system
    int nbeads;         // Number of beads

    long wfreq;         // Frequency of walkers communication
    long sfreq;         // Save frequency (how often the observables are recorded)
    long steps;         // Total number of MD steps

    bool bosonic;       // Is the simulation bosonic?
    bool fixcom;        // Fix the center of mass?
    bool pbc;           // Enable periodic boundary conditions?
    bool nmthermostat;  // Couple thermostat to normal modes

    bool out_pos;       // Output trajectories?
    bool out_vel;       // Output velocities?
    bool out_force;     // Output forces?

    bool is_bosonic_bead; // Is the current simulation bosonic and the time-slice is either 1 or P?
    std::unique_ptr<BosonicExchangeBase> bosonic_exchange;

    std::vector<std::unique_ptr<Observable>> observables;
    std::vector<std::unique_ptr<State>> states;

    std::mt19937 rand_gen;

    Simulation(const int& rank, const int& nproc, Params& param_obj, MPI_Comm& walker_world, MPI_Comm& bead_world, int walker_id, 
               unsigned int seed = static_cast<unsigned int>(time(nullptr)));
    ~Simulation();

    dVec coord, momenta, forces, spring_forces, physical_forces;
    dVec prev_coord, next_coord;

    void setStep(int step);

    double mass;
    double spring_constant;  // k=m*omega_p^2 (where omega_p depends on the convention)
    double omega_p;          // Angular frequency of the ring polymer
    double beta_half_k;      // Pre-factor of beta*0.5*k

    void genRandomPositions(dVec& pos_arr);
    void uniformParticleGrid(dVec& pos_arr) const;
    void genMomentum(dVec& momenta_arr);

    void zeroMomentum();
    
    void initializeWalkersCommunication(Params& param_obj);
    void initializePropagator(Params& param_obj);
    void initializeThermostat(Params& param_obj);
    void initializeExchangeAlgorithm(Params& param_obj);
    void initializePositions(dVec& coord_arr, const VariantMap& sim_params);
    void initializeMomenta(dVec& momentum_arr, const VariantMap& sim_params);
    void addStateIfEnabled(Params& param_obj, const std::string& param_key, const std::string& state_name);
    void initializeStates(Params& param_obj);
    void addObservableIfEnabled(Params& param_obj, const std::string& param_key, const std::string& observable_name);
    void initializeObservables(Params& param_obj);
    std::unique_ptr<Potential> initializePotential(const std::string& potential_name,
                                                   const VariantMap& potential_options);
    void initializeWalkersCommunication(const std::string& walkers_communication_name);

    double sampleMaxwellBoltzmann();
    std::unique_ptr<Propagator> propagator;
    std::unique_ptr<Thermostat> thermostat;
    std::unique_ptr<Coupling> thermostat_coupling;
    std::unique_ptr<WalkersCommunicationBase> walker_communication;

    std::unique_ptr<NormalModes> normal_modes;

    void run();

    std::unique_ptr<Potential> ext_potential;
    std::unique_ptr<Potential> int_potential;
    double int_pot_cutoff;

    std::string external_potential_name;
    std::string interaction_potential_name;

    void updateForces();
    void updateSpringForces();
    void updatePhysicalForces();

    void getNextCoords(dVec& next);
    void getPrevCoords(dVec& prev);
    void updateNeighboringCoordinates();

    int this_bead;   // Current process id ("rank" of MPI_Comm_rank)
    int nproc;       // Number of processes ("size" of MPI_Comm_size)
    unsigned int params_seed;
    std::string thermostat_type;

private:
    int md_step;
    int walker_id;
    void printReport(double wall_time, const std::string& filename) const;

    std::string init_pos_type;
    std::string init_vel_type;
    std::string propagator_type;

    void printDebug(const std::string& text, int target_bead = 0) const;

    MPI_Comm& walker_world; // Walker communicator
    MPI_Comm& bead_world;  // Bead communicator
    std::string walker_communication_type;
};
