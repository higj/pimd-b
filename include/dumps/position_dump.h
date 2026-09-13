#pragma once

#include "dumps/dump.h"
#include "common.h"

#include <memory>

class PositionDump final : public Dump {
public:
    /**
     * @brief Position dump class constructor.
     */
    PositionDump(
        const std::shared_ptr<const VecArray>& coord, 
        int this_bead, 
        int out_freq, 
        const std::string& out_unit
    );

    /**
     * @brief Initializes the coordinates xyz file.
     */
    //void initialize() override;

    /**
     * Outputs the trajectories.
     *
     * @param step Current step of the simulation.
     */
    void output(int step) override;

protected:
    [[nodiscard]] std::string fileName() const override { return std::format("position_{}.xyz", m_this_bead); }

private:
    std::shared_ptr<const VecArray> m_coord;  // Pointer to the coordinates array
    int m_natoms;                             // Number of atoms in the quantum system

#ifdef USE_HDF5
    hid_t m_h5_step_ds = H5I_INVALID_HID;
    hid_t m_h5_pos_ds = H5I_INVALID_HID;
    void h5CreateDatasets() override;
#endif
};