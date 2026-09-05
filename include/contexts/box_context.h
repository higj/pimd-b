#pragma once

#include "common.h"

//struct BoundaryContext {
//    bool pbc;
//    WindingObj winding_obj;
//};

struct BoxContext
{
    double box_size;
    bool pbc;

    // To handle periodic boundary conditions, we employ the Class C storage
    // concept [Z. Phys. Chem. 227 (2013) 345-352], allowing atoms to move
    // outside the primary simulation box. That is, the coordinates are not folded
    // into the simulation box. Instead, we account for the PBC when calculating the
    // distances between particles (or any spatial vector differences).
    // This is done using Algorithm C4. It calculates the remainder of dx
    // on the interval [-L/2, L/2].
    void applyMinimumImage(double& diff) const
    {
#if MINIM
        diff -= box_size * floor(diff / box_size + 0.5);
#endif
    }

    void applyMinimumImageIfNeeded(double& diff) const
    {
        if (pbc)
        {
            applyMinimumImage(diff);
        }
    }

    void applyMinimumImageIfNeeded(Vec& dx_arr) const
    {
        if (pbc)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                applyMinimumImage(dx_arr[axis]);
            }
        }
    }

    void applyMinimumImageIfNeeded(VecArray& dx_arr) const
    {
        if (pbc)
        {
            for (int i = 0; i < dx_arr.len(); ++i)
            {
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    /// TODO: In principle, one could iterate over the elements of the raw (1D) array,
                    ///       but we prefer to keep the nested loops for future generalization of the code
                    ///       to boxes with different side lengths.

                    applyMinimumImage(dx_arr(i, axis));

                    //applyMinimumImageIfNeeded(dx_arr(i, axis));
                    //dx_arr(i, axis) -= box_size * floor(dx_arr(i, axis) / box_size + 0.5);
                }
            }
        }
    }

    /*
    void applyMinimumImage(double& dx, double L) {
        dx -= L * floor(dx / L + 0.5);
    }

    void applyMinimumImage(VecArray& dx_arr, double L)
    {
        for (int i = 0; i < dx_arr.len(); ++i) {
            for (int axis = 0; axis < NDIM; ++axis) {
                applyMinimumImage(dx_arr(i, axis), L);
                //dx_arr(i, axis) -= L * floor(dx_arr(i, axis) / L + 0.5);
            }
        }
    }

    void periodicWrap(double& x, double L) {
        x -= L * std::nearbyint(x / L);
    }
    */
};
