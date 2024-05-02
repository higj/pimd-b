#pragma once

#include "common.h"
#include "bosonic_exchange_base.h"

class BosonicExchange final : public BosonicExchangeBase {
public:
    BosonicExchange(int nbosons_, int np_, int bead_num_, double beta_, double spring_constant_,
                    const dVec& x_, const dVec& x_prev_, const dVec& x_next_, bool pbc_, double size_);
    ~BosonicExchange() override = default;

    double effectivePotential() override;
    double getVn(int n) const;
    double getEknSerialOrder(int i) const;

    void prepare() override;

    double primEstimator() override;

protected:
    void springForceFirstBead(dVec& f) override;
    void springForceLastBead(dVec& f) override;

private:
    void evaluateBosonicEnergies();
    void evaluateCycleEnergies();

    double getEnk(int m, int k) const;
    void setEnk(int m, int k, double val);
    void evaluateConnectionProbabilities();
    void evaluateVBn();
    void evaluateVBackwards();

    std::vector<double> E_kn;
    std::vector<double> V;
    std::vector<double> V_backwards;
    std::vector<double> connection_probabilities;

    std::vector<double> temp_nbosons_array;

    std::vector<double> prim_est;
};
