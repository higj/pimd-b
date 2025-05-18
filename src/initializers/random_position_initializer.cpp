#include "initializers/random_position_initializer.h"
#include "core/random_generators.h"

RandomPositionInitializer::RandomPositionInitializer(
    const std::shared_ptr<RandomGenerators>& rng, 
    const std::shared_ptr<dVec>& coord, 
    double box_size)
    : PositionInitializer(coord, box_size), m_rng(rng) {
}

void RandomPositionInitializer::initialize() {
    double half_box_size = 0.5 * m_box_size;
    //std::uniform_real_distribution<double> u_dist(-half_box_size, half_box_size);

    for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            //coord(ptcl_idx, axis) = u_dist(rand_gen);
            (*m_coord)(ptcl_idx, axis) = m_rng->uniform(-half_box_size, half_box_size);
        }
    }
}