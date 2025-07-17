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

    void applyMinimumImageIfNeeded(double& diff) const
    {
#if MINIM
        if (pbc)
        {
            diff -= box_size * floor(diff / box_size + 0.5);
        }
#endif
    }

    void applyMinimumImageIfNeeded(dVec& dx_arr) const
    {
        for (int i = 0; i < dx_arr.len(); ++i)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                applyMinimumImageIfNeeded(dx_arr(i, axis));
                //dx_arr(i, axis) -= box_size * floor(dx_arr(i, axis) / box_size + 0.5);
            }
        }
    }
};
