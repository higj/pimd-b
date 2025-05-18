#pragma once

#include "dumps/dump.h"
#include "contexts/dumps/position_dump_context.h"

class PositionDump final : public Dump {
public:
    /**
     * @brief Position dump class constructor.
     */
    PositionDump(const PositionDumpContext& dump_context, int out_freq, const std::string& out_unit);

    /**
     * @brief Initializes the coordinates xyz file.
     */
    void initialize() override;

    /**
     * Outputs the trajectories.
     *
     * @param step Current step of the simulation.
     */
    void output(int step) override;

private:
    PositionDumpContext m_context;
};