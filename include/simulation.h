#pragma once

#include <string>
#include <random>
#include <ctime>
#include <memory>

#include "random_mars.h"
#include "common.h"
#include "params.h"
#include "potential.h"
#include "bosonic_exchange_base.h"

class NormalModes;
class State;
class Observable;
class Propagator;
class Thermostat;

class Simulation
{
public:
    double temperature;
    double beta;        // Thermodynamic beta 1/(kB*T)
    double dt;          // Timestep
    double size;        // Linear system size (TODO: Add support for Ly, Lz,...)
    double gamma;       // Friction constant of the Langevin thermostat
    int    nchains;     // Number of Nose-Hoover chains
    double threshold;   // Percentage of steps to throw away (thermalization)

    int natoms;         // Number of atoms in the system
    int nbeads;         // Number of beads

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
    std::unique_ptr<RanMars> mars_gen;

    Simulation(const int& rank, const int& nproc, Params& param_obj, unsigned int seed = static_cast<unsigned int>(time(nullptr)));
    ~Simulation();

    dVec coord, momenta, forces;
    dVec prev_coord, next_coord;

    int getStep() const;
    void setStep(int step);

    double mass;
    double spring_constant;  // k=m*omega_p^2 (where omega_p depends on the convention)
    double omega_p;          // Angular frequency of the ring polymer
    double beta_half_k;      // Pre-factor of beta*0.5*k

    void genRandomPositions(dVec& pos_arr);
    void uniformParticleGrid(dVec& pos_arr) const;
    void genMomentum(dVec& momenta_arr);

    void zeroMomentum();

    void initializePositions(dVec& coord_arr, const VariantMap& sim_params);
    void initializeMomenta(dVec& momentum_arr, const VariantMap& sim_params);
    void addStateIfEnabled(const StringMap& sim_params, const std::string& param_key, const std::string& state_name);
    void initializeStates(const StringMap& sim_params);
    void addObservableIfEnabled(const StringMap& sim_params, const std::string& param_key, const std::string& observable_name);
    void initializeObservables(const StringMap& sim_params);
    std::unique_ptr<Potential> initializePotential(const std::string& potential_name,
                                                   const VariantMap& potential_options);

    double sampleMaxwellBoltzmann();
    
    std::unique_ptr<Propagator> propagator;
    std::unique_ptr<Thermostat> thermostat;

    std::unique_ptr<NormalModes> normal_modes;

    void run();

    std::unique_ptr<Potential> ext_potential;
    std::unique_ptr<Potential> int_potential;
    double int_pot_cutoff;

    std::string external_potential_name;
    std::string interaction_potential_name;

    void updateForces();
    void updateSpringForces(dVec& spring_force_arr) const;
    void updatePhysicalForces(dVec& physical_force_arr) const;

    double classicalSpringEnergy() const;

    void getNextCoords(dVec& next);
    void getPrevCoords(dVec& prev);
    void updateNeighboringCoordinates();

    dVec getSeparation(int first_ptcl, int second_ptcl, bool minimum_image = false) const;

    int this_bead;   // Current process id ("rank" of MPI_Comm_rank)
    int nproc;       // Number of processes ("size" of MPI_Comm_size)
    unsigned int params_seed;
    std::string thermostat_type;

private:
    int md_step;

    void printReport(double wall_time) const;

    std::string init_pos_type;
    std::string init_vel_type;
    std::string propagator_type;

    void printDebug(const std::string& text, int target_bead = 0) const;
};
