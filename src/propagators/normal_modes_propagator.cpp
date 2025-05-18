#include "propagators/normal_modes_propagator.h"
#include "normal_modes.h"
#include "params.h"
#include <numbers>

NormalModesPropagator::NormalModesPropagator(Params& param_obj, dVec& coord, dVec& momenta, dVec& forces,
                                            dVec& physical_forces, dVec& spring_forces, dVec& prev_coord, dVec& next_coord,
                          int this_bead, NormalModes& normal_modes) : 
    Propagator(param_obj, coord, momenta, forces),
    physical_forces(physical_forces),
    spring_forces(spring_forces),
    prev_coord(prev_coord),
    next_coord(next_coord),
    this_bead(this_bead),
    normal_modes(normal_modes)
{
    getVariant(param_obj.sim["nbeads"], nbeads);
    getVariant(param_obj.sim["bosonic"], bosonic);

    double temperature;
    getVariant(param_obj.sys["temperature"], temperature);
#if IPI_CONVENTION
    // i-Pi convention [J. Chem. Phys. 133, 124104 (2010)]
    double omega_p = nbeads * Constants::kB * temperature / Constants::hbar;
#else
    // Tuckerman convention
    double omega_p = sqrt(nbeads) * Constants::kB * temperature / Constants::hbar;
#endif

    // Frequencies
    freq = 2 * omega_p * sin(this_bead * std::numbers::pi / nbeads); 
    c = cos(freq * dt);
    s = sin(freq * dt);
    m_omega = mass * freq;
    spring_constant = mass * omega_p * omega_p;
}

void NormalModesPropagator::preForceStep() {
    // Propagate momenta under external forces
    momentaExternalForces();
    normal_modes.shareData();

    MPI_Barrier(MPI_COMM_WORLD);
    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            const int glob_idx = normal_modes.globIndexAtom(axis, ptcl_idx);
            // Cartesian-to-nm transformation
            double coord_nm = normal_modes.coordCarToNM(glob_idx);
            double momentum_nm = normal_modes.momentumCarToNM(glob_idx);
            // Time propagation
            if (freq == 0) {
                normal_modes.arr_coord_nm[glob_idx + this_bead] = coord_nm + dt / mass * momentum_nm;
                normal_modes.arr_momenta_nm[glob_idx + this_bead] = momentum_nm;
            } else {
                normal_modes.arr_coord_nm[glob_idx + this_bead] = c * coord_nm + s / m_omega * momentum_nm;
                normal_modes.arr_momenta_nm[glob_idx + this_bead] = (-1) * m_omega * s * coord_nm + c * momentum_nm;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            int glob_idx = normal_modes.globIndexAtom(axis, ptcl_idx);
            // NM-to-Cartesian transformation
            double coord_cartesian = normal_modes.coordNMToCar(glob_idx); 
            double momentum_cartesian = normal_modes.momentumNMToCar(glob_idx);
            normal_modes.arr_coord_cartesian[glob_idx + this_bead] = coord_cartesian;
            normal_modes.arr_momenta_cartesian[glob_idx + this_bead] = momentum_cartesian;
        }
    }

    // Update forces
    MPI_Barrier(MPI_COMM_WORLD);
    for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            int glob_idx = normal_modes.globIndexAtom(axis, ptcl_idx);
            coord(ptcl_idx, axis) = normal_modes.arr_coord_cartesian[glob_idx + this_bead];
            momenta(ptcl_idx, axis) = normal_modes.arr_momenta_cartesian[glob_idx + this_bead];
        }
    }
}
void NormalModesPropagator::postForceStep() {
    // Propagate momenta under external forces
    momentaExternalForces();
}

void NormalModesPropagator::momentaExternalForces() {
    if (!bosonic) {
        for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx)
            for (int axis = 0; axis < NDIM; ++axis)
#if IPI_CONVENTION
                momenta(ptcl_idx, axis) += 0.5 * dt * physical_forces(ptcl_idx, axis);
#else
                momenta(ptcl_idx, axis) += 0.5 * dt * physical_forces(ptcl_idx, axis) / nbeads;
#endif
    } else if (this_bead == 0 || this_bead == nbeads - 1) {
        double inner_springs;
        for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                inner_springs = -spring_constant * (2 * coord(ptcl_idx, axis) - prev_coord(ptcl_idx, axis) - next_coord(ptcl_idx, axis));
#if IPI_CONVENTION
                momenta(ptcl_idx, axis) += 0.5 * dt * (physical_forces(ptcl_idx, axis) + spring_forces(ptcl_idx, axis) - inner_springs);
#else
                momenta(ptcl_idx, axis) += 0.5 * dt * (physical_forces(ptcl_idx, axis) / nbeads + spring_forces(ptcl_idx, axis) - inner_springs);
#endif
            }
        }
    } else {
        for (int ptcl_idx = 0; ptcl_idx < natoms; ++ptcl_idx)
            for (int axis = 0; axis < NDIM; ++axis)
#if IPI_CONVENTION
                momenta(ptcl_idx, axis) += 0.5 * dt * physical_forces(ptcl_idx, axis);
#else
                momenta(ptcl_idx, axis) += 0.5 * dt * physical_forces(ptcl_idx, axis) / nbeads;
#endif
    }
}
