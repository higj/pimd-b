#include "initializers/grid_position_initializer.h"

#include <cmath>
#include <stdexcept>

GridPositionInitializer::GridPositionInitializer(
    const std::shared_ptr<dVec>& coord,
    double box_size)
    : PositionInitializer(coord, box_size) {
}

/**
 * Places particles in a uniform grid according to the specified box size.
 */
void GridPositionInitializer::initialize() {
    const double volume = std::pow(m_box_size, NDIM);

    // Get the linear size per particle, and the number of particles
    const double init_side = std::pow((1.0 * m_natoms / volume), -1.0 / (1.0 * NDIM));

    // Determine the number and the size of initial grid boxes in each dimension
    int tot_num_grid_boxes = 1;
    iVec num_nn_grid;
    dVec size_nn_grid;

    for (int i = 0; i < NDIM; i++) {
        num_nn_grid(0, i) = static_cast<int>(std::ceil((m_box_size / init_side) - EPS));

        // Make sure we have at least one grid box
        num_nn_grid(0, i) = std::max(num_nn_grid(0, i), 1);

        // Compute the actual size of the grid
        size_nn_grid(0, i) = m_box_size / (1.0 * num_nn_grid(0, i));

        // Determine the total number of grid boxes
        tot_num_grid_boxes *= num_nn_grid(0, i);
    }

    // Place the particles at the middle of each box
    if (tot_num_grid_boxes < m_natoms) {
        throw std::runtime_error("Number of grid boxes is less than the number of particles");
    }

    dVec pos;
    for (int n = 0; n < tot_num_grid_boxes; n++) {
        iVec grid_index;

        for (int i = 0; i < NDIM; i++) {
            int scale = 1;
            for (int j = i + 1; j < NDIM; j++)
                scale *= num_nn_grid(0, j);

            grid_index(0, i) = (n / scale) % num_nn_grid(0, i);
        }

        for (int axis = 0; axis < NDIM; ++axis) {
            pos(0, axis) = (grid_index(0, axis) + 0.5) * size_nn_grid(0, axis) - 0.5 * m_box_size;
            applyMinimumImage(pos(0, axis), m_box_size);
        }

        if (n >= m_natoms) {
            break;
        }

        for (int axis = 0; axis < NDIM; axis++) {
            (*m_coord)(n, axis) = pos(0, axis);
        }
    }

    /// TODO: Remove this
    //state.updateNeighboringCoordinates();
}