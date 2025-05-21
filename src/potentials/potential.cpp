#include "potentials/potential.h"

Potential::Potential(int start_potential_activation, int finish_potential_activation) 
    : tailV(0.0), start_potential_activation(start_potential_activation), finish_potential_activation(finish_potential_activation) {}

dVec Potential::getGradV(const dVec& x, int md_step) {
    if (md_step < start_potential_activation) {
        return dVec(x.len());
    }

    dVec tempr = gradV(x);
    if (md_step >= finish_potential_activation) {
        return tempr;
    }
    const double prefactor = (md_step - start_potential_activation) / double(finish_potential_activation - start_potential_activation);
    return prefactor * prefactor * tempr;
}

double Potential::getV(const dVec& x, int md_step) {
    if (md_step < start_potential_activation) {
        return 0.;
    }

    double tempr = V(x);
    if (md_step >= finish_potential_activation) {
        return tempr;
    }
    const double prefactor = (md_step - start_potential_activation) / double(finish_potential_activation - start_potential_activation);
    return prefactor * prefactor * tempr;
}
