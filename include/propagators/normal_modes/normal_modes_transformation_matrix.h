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
     * @param simulation Reference to the simulation object containing parameters.
     */
    TransformationMatrixBuilder(const int this_bead, const int nbeads);

    /**
     * @brief Fills the Cartesian-to-NM transformation matrix row.
     * @param cart_to_nm_mat_row Output array to be filled with matrix row values.
     */
    void buildCartToNM(double* cart_to_nm_mat_row) const;

    /**
     * @brief Fills the NM-to-Cartesian transformation matrix row.
     * @param nm_to_cart_mat_row Output array to be filled with matrix row values.
     */
    void buildNMToCart(double* nm_to_cart_mat_row) const;
private:
    const int m_this_bead;
    const int m_nbeads;
    const int m_half_nbeads;
    const double m_fundamental_frequency;
};