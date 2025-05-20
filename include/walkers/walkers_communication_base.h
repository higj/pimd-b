#pragma once
#include "common.h"

class WalkersCommunicationBase {
public:
    explicit WalkersCommunicationBase();
    virtual ~WalkersCommunicationBase() = default;
    
    virtual void communicate(dVec& coord, dVec& momenta);
    double importance_weight, statistical_weight;
};