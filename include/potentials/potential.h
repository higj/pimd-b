#pragma once

#include "common.h"

/* -------------- Basic potential class -------------- */
class Potential
{
public:
    Potential();
    virtual ~Potential() = default;

    // Potential
    virtual double V(const dVec& x)
    {
        double total = 0.0;
        for (int i = 0; i < x.len(); ++i) {
            SingleVec arr{};
            for (int axis = 0; axis < NDIM; ++axis) arr[axis] = x(i, axis);
            total += V(arr);
        }
        return total;
    }

    // Potential gradient
    virtual dVec gradV(const dVec& x)
    {
        dVec result(x.len());
        for (int i = 0; i < x.len(); ++i) {
            SingleVec arr{};
            for (int axis = 0; axis < NDIM; ++axis) arr[axis] = x(i, axis);
            auto grad = gradV(arr);
            for (int axis = 0; axis < NDIM; ++axis) result(i, axis) = grad[axis];
        }
        return result;
    }

    // Potential laplacian
    virtual double laplacianV(const dVec& x)
    {
        double total = 0.0;
        for (int i = 0; i < x.len(); ++i) {
            SingleVec arr{};
            for (int axis = 0; axis < NDIM; ++axis) arr[axis] = x(i, axis);
            total += laplacianV(arr);
        }
        return total;
    }

    // Single-particle overloads (avoid heap allocations)
    virtual double V(const SingleVec& /* x */) { return 0.0; }
    virtual SingleVec gradV(const SingleVec& /* x */) { return {}; }
    virtual double laplacianV(const SingleVec& /* x */) { return 0.0; }

    [[nodiscard]] virtual bool isFree() const { return true; }

    // Tail correction
    /// @todo Implement the tail correction
    double tail_correction;
};
