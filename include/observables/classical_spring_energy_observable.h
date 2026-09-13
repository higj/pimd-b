#pragma once

#include "observables/observable.h"
#include "common.h"
#include "contexts/spring_context.h"
#include "contexts/box_context.h"

#include <memory>

class ClassicalSpringEnergyStrategy;

/**
 * @brief Observable for the classical ring-polymer spring energy.
 *
 * Computes the harmonic spring coupling between adjacent beads in the
 * ring polymer, using the appropriate strategy for distinguishable particles
 * or bosons (selected via StatisticsManager at construction).
 */
class ClassicalSpringEnergyObservable final : public Observable {
public:
    ClassicalSpringEnergyObservable(
        const std::shared_ptr<const VecArray>& coord,
        const std::shared_ptr<const VecArray>& prev_coord,
        const SpringContext& spring_ctx,
        const BoxContext& box_ctx,
        int                  out_freq,
        const std::string& out_unit
    );

    ~ClassicalSpringEnergyObservable() override;

    void calculate() override;

private:
    std::unique_ptr<ClassicalSpringEnergyStrategy> m_spring_energy_strategy;
    std::shared_ptr<const VecArray> m_coord;
    std::shared_ptr<const VecArray> m_prev_coord;
    SpringContext m_spring_ctx;
    BoxContext    m_box_ctx;
};