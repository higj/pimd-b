#include "core/exchange_state.h"

#include "bosonic_exchange.h"
#include "core/simulation_config.h"
#include "contexts/bosonic_exchange_context.h"

void ExchangeState::initialize(const BosonicExchangeContext& context, const bool bosonic) {
    const bool is_bosonic = bosonic && (context.nbeads > 1); // Bosonic exchange is only possible for P>1
    is_bosonic_bead = is_bosonic && (context.this_bead == 0 || context.this_bead == context.nbeads - 1);

    // If the imaginary time-slice is either 1 or P, initialize the bosonic exchange algorithm
    if (is_bosonic_bead) {
#if FACTORIAL_BOSONIC_ALGORITHM
        bosonic_exchange = std::make_unique<FactorialBosonicExchange>(context);
#else
        bosonic_exchange = std::make_unique<BosonicExchange>(context);
#endif
    }
}