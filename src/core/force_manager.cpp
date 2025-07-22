#include <utility>

#include "core/force_manager.h"
#include "core/statistics_manager.h"
#include "potentials.h"
#include "strategies/forces/spring_force_strategy.h"

ForceManager::ForceManager(
    const SimulationConfig& config,
    const std::shared_ptr<BosonicExchangeBase>& bosonic_exchange,
    const BeadContext& bead_ctx
) : m_spring_force_strategy(StatisticsManager::createSpringForceStrategy(
        bosonic_exchange, bead_ctx, config.bosonic
    ))
{
    ext_potential = initializePotential(config, config.ext_pot_name, config.ext_pot_params);
    int_potential = initializePotential(config, config.int_pot_name, config.int_pot_params);

    // If the interaction potential is set to "free", then the cutoff distance is meaningless
    cutoff = (config.int_pot_name == "free") ? 0.0 : std::get<double>(config.int_pot_params.at("cutoff"));

    // For cubic cells with PBC, the cutoff distance must be no greater than L/2 for consistency with
    // the minimum image convention (see 1.6.3 in Allen & Tildesley).
    if (config.pbc)
    {
        cutoff = std::min(cutoff, 0.5 * config.box_size);
    }
}

ForceManager::~ForceManager() = default;

void ForceManager::updatePhysicalForces(SystemState& state, const BoxContext& box_ctx) const
{
    // Calculate the external forces acting on the particles
    state.physical_forces = (-1.0) * ext_potential->gradV(state.coord);

    if (cutoff == 0.0)
        return;

    const int natoms = state.getNumAtoms();

    for (int ptcl_one = 0; ptcl_one < natoms; ++ptcl_one)
    {
        for (int ptcl_two = ptcl_one + 1; ptcl_two < natoms; ++ptcl_two)
        {
            // Get the vector distance between the two particles.
            // Here "diff" contains just one vector of dimension NDIM.
            dVec diff = state.coord.getSeparation(ptcl_one, ptcl_two);
            /// TODO: MINIM should become a parameter (mic_spring and mic_potential)
            box_ctx.applyMinimumImageIfNeeded(diff);

            // If the distance between the particles exceeds the cutoff length
            // then we assume the interaction is negligible and do not bother
            // calculating the force.
            // We use the convention that when cutoff < 0 then the interaction is
            // calculated for all distances.
            if (const double distance = diff.norm(); distance < cutoff || cutoff < 0.0)
            {
                dVec force_on_one = (-1.0) * int_potential->gradV(diff);

                for (int axis = 0; axis < NDIM; ++axis)
                {
                    state.physical_forces(ptcl_one, axis) += force_on_one(0, axis);
                    state.physical_forces(ptcl_two, axis) -= force_on_one(0, axis);
                }
            }
        }
    }
}

void ForceManager::updateForces(SystemState& state, const SpringContext& spring_ctx, const BoxContext& box_ctx) const
{
    // First, update the spring forces based on the current state of the system.
    m_spring_force_strategy->updateSpringForces(state, spring_ctx, box_ctx);

    // Then, update the physical forces acting on the particles.
    updatePhysicalForces(state, box_ctx);
}

std::unique_ptr<Potential> ForceManager::initializePotential(
    const SimulationConfig& config,
    const std::string& potential_name,
    const VariantMap& potential_options
)
{
    if (potential_name == "free")
    {
        return std::make_unique<Potential>();
    }

    if (potential_name == "harmonic")
    {
        double omega = std::get<double>(potential_options.at("omega"));
        return std::make_unique<HarmonicPotential>(config.mass, omega);
    }

    if (potential_name == "double_well")
    {
        double strength = std::get<double>(potential_options.at("strength"));
        double loc = std::get<double>(potential_options.at("location"));
        return std::make_unique<DoubleWellPotential>(config.mass, strength, loc);
    }

    if (potential_name == "dipole")
    {
        double strength = std::get<double>(potential_options.at("strength"));
        return std::make_unique<DipolePotential>(strength);
    }

    if (potential_name == "cosine")
    {
        double amplitude = std::get<double>(potential_options.at("amplitude"));
        double phase = std::get<double>(potential_options.at("phase"));
        return std::make_unique<CosinePotential>(amplitude, config.box_size, phase);
    }

    if (potential_name == "aziz")
    {
        return std::make_unique<AzizPotential>();
    }

    return std::make_unique<Potential>();
}
