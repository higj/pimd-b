#pragma once

#include "common.h"
#include "bosonic_exchange.h"

class BosonicExchangePIS final : public BosonicExchange {
public:
    BosonicExchangePIS(const Simulation& _sim);
    ~BosonicExchangePIS() override = default;

protected:
    void assignIndirectionCoords() override;
    // double current_measure;
    // double measure_shuffle(std::vector<int>& suggested_indexes) const;
    int get_next_index(int current_index, std::vector<int>& suggested_indexes) const;
};
