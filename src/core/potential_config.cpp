#include "core/potential_config.h"

// Concrete potential headers are included here and nowhere else.
// This keeps the compile-time dependency on potential types localised to this file.
#include "potentials/potential.h"
#include "potentials/harmonic_potential.h"
#include "potentials/double_well_potential.h"
#include "potentials/cosine_potential.h"
#include "potentials/dipole_potential.h"
#include "potentials/aziz_potential.h"

#include <stdexcept>

PotentialConfig::PotentialConfig(std::string name, PotentialParams params, double cutoff)
    : m_name(std::move(name)), m_params(std::move(params)), m_cutoff(cutoff) {}

std::unique_ptr<Potential> PotentialConfig::createPotential()
{
    if (m_created)
        throw std::logic_error("PotentialConfig::createPotential() called more than once");
    m_created = true;

    return std::visit([]<typename T>(const T& p) -> std::unique_ptr<Potential>
    {
        if constexpr (std::is_same_v<T, FreePotentialParams>)
            return std::make_unique<Potential>(true);

        else if constexpr (std::is_same_v<T, HarmonicPotentialParams>)
            return std::make_unique<HarmonicPotential>(p.mass, p.omega);

        else if constexpr (std::is_same_v<T, DoubleWellPotentialParams>)
            return std::make_unique<DoubleWellPotential>(p.mass, p.strength, p.location);

        else if constexpr (std::is_same_v<T, CosinePotentialParams>)
            return std::make_unique<CosinePotential>(p.amplitude, p.wavelength, p.phase);

        else if constexpr (std::is_same_v<T, DipolePotentialParams>)
            return std::make_unique<DipolePotential>(p.strength);

        else if constexpr (std::is_same_v<T, AzizPotentialParams>)
            return std::make_unique<AzizPotential>();

    }, m_params);
}
