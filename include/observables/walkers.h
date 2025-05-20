#pragma once

#include "observables/observable.h"

class WalkersCommunicationBase; // Forward declaration

/* -------------------------------- */

class WalkersObservable : public Observable {
public:
    WalkersObservable(Params& param_obj, int _freq, const std::string& _out_unit, int this_bead, WalkersCommunicationBase& walker_communication);

    void calculate() override;
private:
    WalkersCommunicationBase& walker_communication;
};