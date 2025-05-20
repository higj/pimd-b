#include "walkers/walkers_communication_base.h"

WalkersCommunicationBase::WalkersCommunicationBase() {}

void WalkersCommunicationBase::communicate(dVec& coord, dVec& momenta) {
    importance_weight = 0.0;
    statistical_weight = 1.0;
}