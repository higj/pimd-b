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

    void applyMinimumImageIfNeeded(dVec& dx_arr) const
    {
        for (int i = 0; i < dx_arr.len(); ++i)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                /// TODO: In principle, one could iterate over the elements of the raw (1D) array,
                ///       but we prefer to keep the nested loops for future generalization of the code
                ///       to boxes with different side lengths.
                applyMinimumImageIfNeeded(dx_arr(i, axis));
                //dx_arr(i, axis) -= box_size * floor(dx_arr(i, axis) / box_size + 0.5);
            }
        }
    }
};
