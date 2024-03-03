#pragma once

#include <string>
#include <unordered_map>
#include <numeric>
#include <random>
#include <cmath>
#include <ctime>
#include <optional>
#include <functional>

#include <iostream>
#include <fstream>

#include <memory>
#include "mpi.h"

#include "common.h"
#include "params.h"
#include "potential.h"

#include "random_mars.h"
#include "random_park.h"

#include "bosonic_exchange_base.h"


class Observable;

class Simulation
{
public:
	double temperature;
	double beta;      // Thermodynamic beta 1/(kB*T)
	double dt;        // Timestep
	double size;      // Linear system size (TODO: Add support for Ly, Lz,...)
	double gamma;     // Friction constant of the Langevin thermostat
	double threshold; // Percentage of steps to throw away (thermalization)
	
	int natoms;       // Number of atoms in the system
	int nbeads;       // Number of beads

	long sfreq;       // Save frequency (how often the observables are recorded)
	long steps;       // Total number of MD steps

	bool enable_t;    // Enable the thermostat?
	bool bosonic;     // Is the simulation bosonic?
	bool fixcom;      // Fix the center of mass?
	bool pbc;         // Enable periodic boundary conditions?

	bool out_pos;     // Output trajectories?
	bool out_vel;     // Output velocities?
	bool out_force;   // Output forces?

	std::unique_ptr<BosonicExchangeBase> bosonic_exchange;

	std::vector<std::unique_ptr<Observable>> observables;

	std::mt19937 rand_gen;
	std::unique_ptr<RanMars> mars_gen;

	Simulation(int &rank, int &nproc, Params &paramObj, unsigned int seed = static_cast<unsigned int>(time(nullptr)));
	~Simulation();

	dVec coord, momenta, forces;
	dVec prev_coord, next_coord;

	double mass;
	double spring_constant;  // k=m*omega_p^2 (where omega_p depends on the convention)
	double omega_p;          // Angular frequency of the ring polymer

	void genRandomPositions(dVec& pos_arr);
	void genMomentum(dVec &momenta_arr);

	void zeroMomentum();
	
	double sampleMB();

	void langevinStep();
	void vvStep();
	void run();
	void forceIniCond(dVec pos_arr, dVec momentum_arr);

	std::unique_ptr<Potential> ext_potential;
	std::unique_ptr<Potential> int_potential;
	double int_pot_cutoff;

	void updateForces();
	void updateSpringForces(dVec &spring_force_arr);
	
	void getNextCoords(dVec& next);
	void getPrevCoords(dVec& prev);
	void updateNeighboringCoordinates();

	void outputTrajectories(int step);
	void outputVelocities(int step);
	void outputForces(int step);

	dVec getSeparation(int first_ptcl, int second_ptcl) const;

	int this_bead;  // Current process id ("rank" of MPI_Comm_rank)
	int nproc;      // Number of processes ("size" of MPI_Comm_size)
	unsigned int params_seed;

private:		
	void printReport(std::ofstream& out_file) const;

	std::string init_pos_type;
	std::string init_vel_type;

	std::string external_potential_name;
	std::string interaction_potential_name;

	void printDebug(const std::string& text);
};
