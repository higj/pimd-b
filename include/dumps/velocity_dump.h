#pragma once

#include "dumps/dump.h"
#include "contexts/velocity_context.h"

class VelocityDump final : public Dump {
public:
    /**
     * @brief Velocity dump class constructor.
     */
    VelocityDump(const VelocityContext& dump_context, int this_bead, int out_freq, const std::string& out_unit);

    /**
     * @brief Initializes the velocities dat file.
     */
    //void initialize() override;

    /**
     * Outputs the velocities.
     *
     * @param step Current step of the simulation.
     */
    void output(int step) override;

protected:
    [[nodiscard]] std::string fileName() const override { return std::format("velocity_{}.xyz", m_this_bead); }

private:
    VelocityContext m_context;
    int m_natoms;
};