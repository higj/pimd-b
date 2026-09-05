#pragma once

#include "initializable_item.h"

#include <cstdint>

enum class ObservableType : std::uint8_t {
    PRIM_KE,
    VIRIAL_KE,
    POT_ENERGY,
    EXT_POT,
    INT_POT,
    CLASSICAL,
    CL_KINETIC_ENERGY,
    CL_SPRING_ENERGY,
    QUANTUM_TEMPERATURE,
    NOSE_HOOVER_ENERGY,
    BOSONIC_PROB_DIST,
    BOSONIC_PROB_ALL,
    BOSONIC_SIGN,
    SC_WEIGHT,             // "sc_weight" or "w_gsf"  - Suzuki-Chin log weight
    SC_POT_ENERGY,         // "sc_pot"                - odd-bead potential estimator
    SC_EVEN_POT_ENERGY,    // "sc_even_pot"           - even-bead potential estimator
    SC_KINETIC_ENERGY,     // "sc_kinetic"            - force-squared kinetic correction
    SC_VIRIAL_ENERGY,      // "sc_virial"             - GSF virial estimator
    CENTER_OF_MASS,
    UNKNOWN
};

class ObservableItem final : public InitializableItem<ObservableType> {
public:
    ObservableItem(const std::string& obs_name, const std::string& obs_unit) {
        name = obs_name;
        unit = obs_unit;
    }

    [[nodiscard]] ObservableType getType() const override {
        if (name == "kinetic") return ObservableType::PRIM_KE;
        if (name == "virial") return ObservableType::VIRIAL_KE;
        if (name == "potential") return ObservableType::POT_ENERGY;
        if (name == "ext_pot") return ObservableType::EXT_POT;
        if (name == "int_pot") return ObservableType::INT_POT;
        if (name == "classical") return ObservableType::CLASSICAL;
        if (name == "cl_kinetic") return ObservableType::CL_KINETIC_ENERGY;
        if (name == "cl_spring") return ObservableType::CL_SPRING_ENERGY;
        if (name == "temperature") return ObservableType::QUANTUM_TEMPERATURE;
        if (name == "nh_energy") return ObservableType::NOSE_HOOVER_ENERGY;
        if (name == "prob_dist") return ObservableType::BOSONIC_PROB_DIST;
        if (name == "prob_all") return ObservableType::BOSONIC_PROB_ALL;
        if (name == "sign") return ObservableType::BOSONIC_SIGN;
        if (name == "w_gsf") return ObservableType::SC_WEIGHT;
        if (name == "pot_gsf") return ObservableType::SC_POT_ENERGY;
        if (name == "even_pot_gsf") return ObservableType::SC_EVEN_POT_ENERGY;
        if (name == "kin_gsf") return ObservableType::SC_KINETIC_ENERGY;
        if (name == "virial_gsf") return ObservableType::SC_VIRIAL_ENERGY;
        if (name == "center_of_mass") return ObservableType::CENTER_OF_MASS;
        return ObservableType::UNKNOWN;
    }
};
