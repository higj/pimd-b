#include "ring_polymer_utils.h"

namespace RingPolymerUtils
{
    double classicalSpringEnergy(const dVec& coord, const dVec& prev_coord, double spring_constant, bool minimum_image, double box_size) {
        //assert(!m_context.config->bosonic || (m_context.config->bosonic && m_context.config->this_bead != 0));

        double interior_spring_energy = 0.0;

        for (int ptcl_idx = 0; ptcl_idx < coord.len(); ++ptcl_idx) {
            for (int axis = 0; axis < NDIM; ++axis) {
                double diff = prev_coord(ptcl_idx, axis) - coord(ptcl_idx, axis);

                if (minimum_image) {
                    applyMinimumImage(diff, box_size);
                }

                interior_spring_energy += diff * diff;
            }
        }

        interior_spring_energy *= 0.5 * spring_constant;

        return interior_spring_energy;
    }
}