#pragma once

#include "dumps/dump.h"
#include "contexts/dumps/velocity_dump_context.h"

class VelocityDump final : public Dump {
public:
    /**
     * @brief Velocity dump class constructor.
     */
    VelocityDump(const VelocityDumpContext& dump_context, int out_freq, const std::string& out_unit);

    /**
     * @brief Initializes the velocities dat file.
     */
    void initialize() override;

    /**
     * Outputs the velocities.
     *
     * @param step Current step of the simulation.
     */
    void output(int step) override;

private:
    VelocityDumpContext m_context;
};