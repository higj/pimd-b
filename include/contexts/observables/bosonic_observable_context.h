#pragma once

#include <memory>

class ExchangeState;

/**
 * Parameters unique to the bosonic observable.
 */
struct BosonicObservableContext {
    int this_bead;
    std::shared_ptr<ExchangeState> exchange_state;
};