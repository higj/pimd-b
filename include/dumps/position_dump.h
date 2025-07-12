#pragma once

#include "dumps/dump.h"
#include "common.h"

class PositionDump final : public Dump {
public:
    /**
     * @brief Position dump class constructor.
     */
    PositionDump(const std::shared_ptr<const dVec>& coord, int this_bead, int out_freq, const std::string& out_unit);

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
    std::shared_ptr<const dVec> m_coord;  // Pointer to the coordinates array
    int m_natoms;                         // Number of atoms in the quantum system
};