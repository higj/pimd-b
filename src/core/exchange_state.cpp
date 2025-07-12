#include "core/exchange_state.h"

#include "bosonic_exchange.h"
#include "core/simulation_config.h"
#include "contexts/bosonic_exchange_context.h"

ExchangeState::ExchangeState(
    const BosonicExchangeContext& context,
    bool bosonic) :
    is_bosonic(bosonic),
    is_bosonic_bead(bosonic && (context.this_bead == 0 || context.this_bead == context.nbeads - 1))
{
    // If the imaginary time-slice is either 1 or P, initialize the bosonic exchange algorithm
    // CR: otherwise it remains uninitialized?
    // JH: Yes
    if (is_bosonic_bead)
    {
#if FACTORIAL_BOSONIC_ALGORITHM
        bosonic_exchange = std::make_unique<FactorialBosonicExchange>(context);
#else
        bosonic_exchange = std::make_unique<BosonicExchange>(context);
#endif
    }
}
