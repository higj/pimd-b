#include "propagators/normal_modes/normal_modes_transformation_matrix.h"

#include <cmath>
#include <algorithm>
#include <numbers>

TransformationMatrixBuilder::TransformationMatrixBuilder(const int this_bead, const int nbeads)
    : m_this_bead(this_bead), m_nbeads(nbeads), m_half_nbeads(nbeads / 2), 
    m_fundamental_frequency(2 * std::numbers::pi / nbeads * this_bead) {
}

void TransformationMatrixBuilder::buildCartToNM(double* cart_to_nm_mat_row) const {
    double prefactor;

    if (m_this_bead == 0) {
        // Zero mode - constant value
        prefactor = 1 / sqrt(m_nbeads);
        std::ranges::fill(cart_to_nm_mat_row, cart_to_nm_mat_row + m_nbeads, prefactor);
        return;
    }

    if (m_this_bead < m_half_nbeads) {
        // Cosine modes
        prefactor = sqrt(2.0 / m_nbeads);
        for (int i = 0; i < m_nbeads; ++i) {
            cart_to_nm_mat_row[i] = prefactor * cos(m_fundamental_frequency * i);
        }
        return;
    }

    if (m_this_bead == m_half_nbeads) {
        // Nyquist mode (if nbeads is even)
        prefactor = 1 / sqrt(m_nbeads);
        for (int i = 0; i < m_nbeads; ++i) {
            cart_to_nm_mat_row[i] = prefactor * (i % 2 == 0 ? 1.0 : -1.0);
        }
        return;
    }

    // Sine modes
    prefactor = sqrt(2.0 / m_nbeads);
    for (int i = 0; i < m_nbeads; ++i) {
        cart_to_nm_mat_row[i] = -prefactor * sin(m_fundamental_frequency * i);
    }
}

void TransformationMatrixBuilder::buildNMToCart(double* nm_to_cart_mat_row) const {
    const double prefactor = sqrt(2.0 / m_nbeads);

    // Zero mode
    nm_to_cart_mat_row[0] = 1 / sqrt(m_nbeads);

    // Cosine modes
    for (int i = 1; i < m_half_nbeads; ++i) {
        nm_to_cart_mat_row[i] = prefactor * cos(m_fundamental_frequency * i);
    }

    // Nyquist mode (if nbeads is even)
    if (m_nbeads % 2 == 0) {
        nm_to_cart_mat_row[m_nbeads / 2] = (1 / sqrt(m_nbeads)) *
            (m_this_bead % 2 == 0 ? 1.0 : -1.0);
    }

    // Sine modes
    for (int i = std::ceil(m_half_nbeads + 0.5); i < m_nbeads; ++i) {
        nm_to_cart_mat_row[i] = -prefactor * sin(m_fundamental_frequency * i);
    }
}