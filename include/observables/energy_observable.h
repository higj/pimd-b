#pragma once

#include "observables/observable.h"
#include "common.h"
#include "contexts/bead_context.h"
#include "contexts/thermal_context.h"
#include "contexts/spring_context.h"
#include "contexts/box_context.h"

#include <memory>

class ForceManager;
class PrimitiveKineticEnergyStrategy;

class EnergyObservable final : public Observable {
public:
    /**
     * @brief Energy observable class constructor.
     */
    EnergyObservable(
        const std::shared_ptr<const VecArray>& coord,
        const std::shared_ptr<const VecArray>& prev_coord,
        const std::shared_ptr<const VecArray>& physical_forces,
        const std::shared_ptr<const ForceManager>& force_mgr,
        const BeadContext& bead_ctx,
        const ThermalContext& thermal_ctx,
        const SpringContext& spring_ctx,
        const BoxContext& box_ctx,
        int out_freq, 
        const std::string& out_unit
    );

    ~EnergyObservable() override;

    /**
     * @brief Calculates the kinetic and potential energies.
     */
    void calculate() override;

private:
    std::unique_ptr<PrimitiveKineticEnergyStrategy> m_prim_ke_strategy;
    std::shared_ptr<const VecArray> m_coord_this;
    std::shared_ptr<const VecArray> m_coord_prev;
    std::shared_ptr<const VecArray> m_physical_forces;
    std::shared_ptr<const ForceManager> m_force_mgr;

    BeadContext m_bead_ctx;
    ThermalContext m_thermal_ctx;
    SpringContext m_spring_ctx;
    BoxContext m_box_ctx;

    bool m_is_ext_free; // Is the external potential free?
    bool m_is_int_free; // Is the interaction potential free?

    /**
     * @brief Calculates the quantum kinetic energy of the system using the primitive kinetic energy estimator.
     * Works both for distinguishable particles and bosons.
     */
    void calculateKinetic();

    /**
     * @brief Calculates the quantum potential energy of the system, based on the potential energy estimator.
     * The potential energy is the sum of the external potential energy and the interaction potential energy
     * across all time-slices, divided by the number of beads. In addition, the method calculates the virial
     * kinetic energy of the system.
     */
    void calculatePotential();
};