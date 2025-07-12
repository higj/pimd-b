#pragma once

#include <memory>

#include "bosonic_exchange/bosonic_exchange_base.h"

struct BosonicExchangeContext;

// Holds data related to quantum exchange
class ExchangeState {
public:
    ExchangeState(
        const std::shared_ptr<const dVec>& coord_first_bead,
        const std::shared_ptr<const dVec>& coord_last_bead,
        const ThermalContext& thermal_ctx,
        const SpringContext& spring_ctx,
        const BoxContext& box_ctx,
        int nbeads,
        int this_bead,
        bool bosonic
    );

    bool is_bosonic;       // Is the current simulation bosonic?
    bool is_bosonic_bead;  // Is the current simulation bosonic and the time-slice is either 1 or P?
    std::unique_ptr<BosonicExchangeBase> bosonic_exchange;
};
