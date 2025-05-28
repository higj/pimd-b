#pragma once

/**
 * @class TransformationMatrixBuilder
 * @brief Handles the creation of transformation matrices for normal mode calculations.
 *
 * This class encapsulates the logic for building Cartesian-to-NM and NM-to-Cartesian
 * transformation matrix rows based on simulation parameters.
 */
class TransformationMatrixBuilder {
public:
    /**
     * @brief Constructor that initializes the builder with simulation parameters.
     * @param this_bead The index of the current bead in the simulation.
     * @param nbeads Total number of beads in the simulation.
     */
    TransformationMatrixBuilder(int this_bead, int nbeads);

    /**
     * @brief Fills the Cartesian-to-NM transformation matrix row.
     * @param cart_to_nm_mat_row Output array to be filled with matrix row values.
     */
    void buildCartesianToNormalModes(double* cart_to_nm_mat_row) const;

    /**
     * @brief Fills the NM-to-Cartesian transformation matrix row.
     * @param nm_to_cart_mat_row Output array to be filled with matrix row values.
     */
    void buildNormalModesToCartesian(double* nm_to_cart_mat_row) const;
private:
    int m_this_bead;
    int m_nbeads;
    int m_half_nbeads;
    double m_fundamental_frequency;
};