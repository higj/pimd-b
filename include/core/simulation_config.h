#pragma once

#include "common.h"
#include <string>

// Holds immutable configuration parameters, parsed from the configuration file
struct SimulationConfig {
    double temperature;
    double beta;        // Thermodynamic beta 1/(kB*T)
    double thermo_beta; // Inverse temperature at which the simulation is actually performed ("thermostat beta")
    double dt;          // Timestep
    double box_size;    // Linear system size (TODO: Add support for Ly, Lz,...)
    double threshold;   // Percentage of steps to throw away (thermalization)

    double mass;
    double spring_constant;  // k=m*omega_p^2 (where omega_p depends on the convention)
    double omega_p;          // Angular frequency of the ring polymer
    double beta_half_k;      // Pre-factor of beta*0.5*k

    int natoms;         // Number of atoms in the system
    int nbeads;         // Number of beads

    long sfreq;         // Save frequency (how often the observables are recorded)
    long steps;         // Total number of MD steps

    bool bosonic;       // Is the simulation bosonic?
    bool fixcom;        // Fix the center of mass?
    bool pbc;           // Enable periodic boundary conditions?

    // Thermostat related parameters
    std::string thermostat_type;
    VariantMap thermostat_params;

    bool out_pos;       // Output trajectories?
    bool out_vel;       // Output velocities?
    bool out_force;     // Output forces?

    // Initialization method name for positions and velocities
    std::string init_pos_type, init_vel_type;
    std::string init_pos_filename, init_vel_filename;
    int init_pos_index_offset, init_vel_index_offset;

    // Map holding the interaction potential parameters
    VariantMap int_pot_params;
    // Map holding the interaction potential parameters
    VariantMap ext_pot_params;

    // Propagator and thermostat types
    std::string propagator_type;

    // External and interaction potential names
    std::string ext_pot_name, int_pot_name;

    // Map holding the dump parameters
    StringMap dumps_list;
    // Map holding the observable settings
    StringMap observables_list;

    unsigned int seed;  // Seed for random number generation
    int this_bead;      // Current process id ("rank" of MPI_Comm_rank)
};
