#pragma once

#include "dumps/dump.h"

class SystemState;

class ForceDump final : public Dump {
public:
    /**
     * @brief Force dump class constructor.
     */
    ForceDump(const std::shared_ptr<SystemState>& state, int this_bead, int out_freq, const std::string& out_unit);

    /**
     * @brief Initializes the forces dat file.
     */
    void initialize() override;

    /**
     * Outputs the forces.
     *
     * @param step Current step of the simulation.
     */
    void output(int step) override;

private:
    std::shared_ptr<SystemState> m_state;
};