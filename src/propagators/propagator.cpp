#include "propagators/propagator.h"

Propagator::Propagator(Params& param_obj, dVec& coord, dVec& momenta, dVec& forces) : 
    coord(coord), momenta(momenta), forces(forces) {
    getVariant(param_obj.sys["natoms"], natoms);
    getVariant(param_obj.sim["dt"], dt);
    getVariant(param_obj.sys["mass"], mass);
}