#pragma once

#include "common.h"
#include "bosonic_exchange_base.h"

class BosonicExchange : public BosonicExchangeBase {
public:
    BosonicExchange(int nbosons, int np, int bead_num, double beta, double spring_constant,
        const dVec x, const dVec x_prev, const dVec x_next, bool pbc, double L);
    ~BosonicExchange();

    double effectivePotential() override;
    double get_Vn(int n) const;
    double get_E_kn_serial_order(int i) const;

    void updateCoordinates(const dVec new_x, const dVec new_x_prev, const dVec new_x_next) override;

    double primEstimator();

protected:
    void springForceFirstBead(dVec& f) override;
    void springForceLastBead(dVec& f) override;

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