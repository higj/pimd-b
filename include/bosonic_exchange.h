#pragma once

#include "mpi.h"
#include "common.h"

class BosonicExchange {
public:
    BosonicExchange(int nbosons, int np, int bead_num, double beta, double spring_constant,
        const dVec x, const dVec x_prev, const dVec x_next, bool pbc, double L);
    ~BosonicExchange();

    double get_potential() const;
    double get_Vn(int n) const;
    double get_E_kn_serial_order(int i) const;

    void spring_force(dVec &f);

    double prim_estimator();

private:
    void prepare_with_coordinates();
    void evaluate_cycle_energies();
    void diff_two_beads(const dVec x1, int l1, const dVec x2, int l2, double diff[NDIM]);
    double distance_squared_two_beads(const dVec x1, int l1, const dVec x2, int l2);
    double get_Enk(int m, int k);
    void set_Enk(int m, int k, double val);
    void evaluate_connection_probabilities();
    void spring_force_last_bead(dVec& f);
    void spring_force_first_bead(dVec& f);
    void spring_force_interior_bead(dVec& f);
    void Evaluate_VBn();
    void Evaluate_V_backwards();

    const int nbosons;
    const int np;
    const int bead_num;

    double spring_constant;
    double beta;
    const dVec x;
    const dVec x_prev;
    const dVec x_next;

    bool pbc;
    double L; // Linear size of the system

    double spring_energy;

    std::vector<double> E_kn;
    std::vector<double> V;
    std::vector<double> V_backwards;
    std::vector<double> connection_probabilities;

    std::vector<double> temp_nbosons_array;
    std::vector<double> separate_atom_spring;

    std::vector<double> prim_est;
};