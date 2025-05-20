#pragma once

#include "common.h"
#include "bosonic_exchange/bosonic_exchange_base.h"

class BosonicExchange final : public BosonicExchangeBase {
public:
    BosonicExchange(Params& param_obj, const dVec& coord, const dVec& prev_coord, const dVec& next_coord, const int this_bead);
    ~BosonicExchange() override = default;

    double effectivePotential() override;
    double getVn(int n) const;
    double getEknSerialOrder(int i) const;

    void prepare() override;
    double primEstimator() override;

    double getDistinctProbability() override;
    double getLongestProbability() override;

    void printBosonicDebug() override;

    double getCumulativeCycleProb() const override;

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

    double log_n_factorial;
};
