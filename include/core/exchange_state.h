#pragma once

#include <memory>

#include "bosonic_exchange/bosonic_exchange_base.h"

struct BosonicExchangeContext;

// Holds data related to quantum exchange
class ExchangeState {
public:
    ExchangeState(const BosonicExchangeContext& context, bool bosonic);

    bool is_bosonic;       // Is the current simulation bosonic?
    bool is_bosonic_bead;  // Is the current simulation bosonic and the time-slice is either 1 or P?
    std::unique_ptr<BosonicExchangeBase> bosonic_exchange;
};
