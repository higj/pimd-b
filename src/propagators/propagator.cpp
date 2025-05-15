#include "propagators/propagator.h"
#include "params.h"

Propagator::Propagator(Simulation& _sim, Params& param_obj, dVec& coord, dVec& momenta, dVec& forces) : 
sim(_sim), coord(coord), momenta(momenta), forces(forces) {
    getVariant(param_obj.sys["natoms"], natoms);
    getVariant(param_obj.sim["dt"], dt);
    getVariant(param_obj.sys["mass"], mass); 
}