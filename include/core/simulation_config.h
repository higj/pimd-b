#pragma once

#include "common.h"

#include <string>
#include "core/potential_config.h"

struct ThermostatConfig {
    std::string type;
    bool couple_to_nm = false;
    double gamma = 0.0;
    bool is_nose_hoover = false;
    int nchains = 0;
};

// RPMD-related parameters
struct RpmdConfig {
    bool enabled;             // Is RPMD mode active?
    int num_runs;             // Number of independent simulations to perform
    double nvt_discard_frac;  // Fraction of the NVT trajectory to discard before sampling starts
};

enum class XyzFrameSelectionMode {
    Index,
    Step
};

// Holds immutable configuration parameters, parsed from the configuration file
struct SimulationConfig {
    double temperature;
    double beta;        // Thermodynamic beta 1/(kB*T)
    double thermo_beta; // Inverse temperature at which the simulation is actually performed ("thermostat beta")
    double dt;          // Timestep
    double box_size;    // Linear system size (TODO: Add support for Ly, Lz,...)
    long threshold;     // Thermalization up to this many steps

    double mass;
    double spring_constant;  // k=m*omega_p^2 (where omega_p depends on the convention)
    double omega_p;          // Angular frequency of the ring polymer
    double beta_half_k;      // Pre-factor of beta*0.5*k

    int natoms;         // Number of atoms in the system
    int nbeads;         // Number of beads

    long sfreq;         // Save frequency (how often the observables are recorded)
    long steps;         // Total number of MD steps

    bool bosonic;        // Is the simulation bosonic?
    double exchange_xi;  // Exchange parameter (only relevant when bosonic=true)

    bool fixcom;        // Fix the center of mass?
    bool pbc;           // Enable periodic boundary conditions?

    // Thermostat related parameters
    ThermostatConfig thermostat;

    bool out_pos;       // Output trajectories?
    bool out_vel;       // Output velocities?
    bool out_force;     // Output forces?

    // Initialization method name for positions and velocities
    std::string init_pos_type, init_vel_type;
    std::string init_pos_unit, init_vel_unit;
    long init_pos_frame, init_vel_frame;
    XyzFrameSelectionMode init_pos_frame_mode = XyzFrameSelectionMode::Index;
    XyzFrameSelectionMode init_vel_frame_mode = XyzFrameSelectionMode::Index;
    std::string init_pos_filename, init_vel_filename;
    int init_pos_index_offset, init_vel_index_offset;

    // Configuration objects that carry parameters and create the concrete Potential
    PotentialConfig ext_potential_cfg;
    PotentialConfig int_potential_cfg;

    // Propagator type
    std::string propagator_type;

    // Map holding the dump parameters
    StringMap dumps_list;
    // Map holding the observable settings
    StringMap observables_list;

    // Map holding the original user-requested units for the various physical quantities
    StringMap units_list;

    // RPMD-related parameters
    RpmdConfig rpmd_config;

    unsigned int seed;  // Seed for random number generation
    int this_bead;      // Current process id ("rank" of MPI_Comm_rank)
};
