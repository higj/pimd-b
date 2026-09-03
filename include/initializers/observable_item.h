#pragma once

#include "initializable_item.h"

#include <cstdint>

enum class ObservableType : std::uint8_t {
    PRIM_KE,
    VIRIAL_KE,
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
    GSF,
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
        if (name == "gsf") return ObservableType::GSF;
        if (name == "center_of_mass") return ObservableType::CENTER_OF_MASS;
        return ObservableType::UNKNOWN;
    }
};
