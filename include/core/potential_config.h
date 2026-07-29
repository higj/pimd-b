#pragma once

#include <memory>
#include <string>
#include <variant>

// Forward declaration — concrete headers are only needed in potential_config.cpp
class Potential;

// ============================================================
// Per-potential parameter structs
// Each struct carries only the parameters relevant to that potential.
// Adding a new potential requires only:
//   1. A new struct here
//   2. A new alternative in PotentialParams
//   3. A new branch in PotentialConfig::createPotential()
// ============================================================

struct FreePotentialParams {};

struct HarmonicPotentialParams {
    double mass;
    double omega;
};

struct AnharmonicPotentialParams {
    double mass;
    double omega;
    double cubic_const;
    double quart_const;
};

struct DoubleWellPotentialParams {
    double mass;
    double strength;
    double location;
};

struct CosinePotentialParams {
    double amplitude;
    double wavelength;
    double phase;
};

struct DipolePotentialParams {
    double strength;
};

struct AzizPotentialParams {};

// Type-safe union: holds exactly one alternative at a time.
// Only the active potential's parameters are stored.
using PotentialParams = std::variant<
    FreePotentialParams,
    HarmonicPotentialParams,
    AnharmonicPotentialParams,
    DoubleWellPotentialParams,
    CosinePotentialParams,
    DipolePotentialParams,
    AzizPotentialParams
>;

// ============================================================
// PotentialConfig
// Responsible for carrying the parameters of the chosen potential
// and for instantiating the concrete Potential object on demand.
// ============================================================
class PotentialConfig {
public:
    /// Default-constructs a "free" potential with no cutoff.
    PotentialConfig() = default;

    /**
     * @param name   Human-readable name of the potential (e.g. "harmonic").
     *               Used in the simulation report.
     * @param params Parameters for the specific potential type.
     * @param cutoff Interaction cutoff radius (internal units).
     *               Negative means no cutoff. Relevant only for interaction potentials.
     */
    PotentialConfig(std::string name, PotentialParams params, double cutoff = -1.0);

    /**
     * Construct and return the concrete Potential object.
     * @throws std::logic_error if called more than once.
     */
    [[nodiscard]] std::unique_ptr<Potential> createPotential();

    [[nodiscard]] const std::string& name()   const { return m_name;   }
    [[nodiscard]] double             cutoff() const { return m_cutoff; }

private:
    std::string     m_name    = "free";
    PotentialParams m_params  = FreePotentialParams{};
    double          m_cutoff  = -1.0;
    bool            m_created = false;
};
