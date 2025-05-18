#pragma once

#include "dumps/dump.h"
#include "contexts/dumps/force_dump_context.h"

class ForceDump final : public Dump {
public:
    /**
     * @brief Force dump class constructor.
     */
    ForceDump(const ForceDumpContext& dump_context, int out_freq, const std::string& out_unit);

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
    ForceDumpContext m_context;
};