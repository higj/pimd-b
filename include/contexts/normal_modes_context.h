#pragma once

#include <memory>

class NormalModes;

struct NormalModesContext {
    std::shared_ptr<NormalModes> normal_modes;
    bool couple_to_nm;
};