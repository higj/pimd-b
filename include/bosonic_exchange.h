#pragma once

#include "common.h"
#include "bosonic_exchange_base.h"

class BosonicExchange : public BosonicExchangeBase {
public:
    BosonicExchange(int nbosons, int np, int bead_num, double beta, double spring_constant,
        const dVec x, const dVec x_prev, const dVec x_next, bool pbc, double L);
    ~BosonicExchange();

    double get_potential() const;
    double get_Vn(int n) const;
    double get_E_kn_serial_order(int i) const;

    void updateCoordinates(const dVec new_x, const dVec new_x_prev, const dVec new_x_next) override;
    
    double classical_potential();
    double prim_estimator();

protected:
    void spring_force_first_bead(dVec& f) override;
    void spring_force_last_bead(dVec& f) override;

private:
    void prepare_with_coordinates();
    void evaluate_cycle_energies();
    
    double get_Enk(int m, int k);
    void set_Enk(int m, int k, double val);
    void evaluate_connection_probabilities();
    void Evaluate_VBn();
    void Evaluate_V_backwards();

    std::vector<double> E_kn;
    std::vector<double> V;
    std::vector<double> V_backwards;
    std::vector<double> connection_probabilities;

    std::vector<double> temp_nbosons_array;
    std::vector<double> separate_atom_spring;

    std::vector<double> prim_est;
};