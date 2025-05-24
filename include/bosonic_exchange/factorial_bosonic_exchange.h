#pragma once

#include "bosonic_exchange/bosonic_exchange_base.h"

#include <vector>

class FactorialBosonicExchange final : public BosonicExchangeBase {
public:
    explicit FactorialBosonicExchange(const BosonicExchangeContext& context);
    ~FactorialBosonicExchange() override = default;

    double effectivePotential() override;
    void prepare() override;
    double primEstimator() override;

    double getDistinctProbability() override;
    double getLongestProbability() override;

    void printBosonicDebug() override;
protected:
    void springForceFirstBead(dVec& f) override;
    void springForceLastBead(dVec& f) override;

private:
    [[nodiscard]] int firstBeadNeighbor(int ptcl_idx) const;
    [[nodiscard]] int lastBeadNeighbor(int ptcl_idx) const;

    double getMinExteriorSpringEnergy();

    std::vector<int> m_labels;  // Particle labels

    double m_e_shift;
};
