#pragma once

#include "common.h"

/* -------------- Basic potential class -------------- */
class Potential
{
public:
    Potential();
    virtual ~Potential() = default;

    // Potential
    virtual double V(const VecArray& x)
    {
        double total = 0.0;
        for (int i = 0; i < x.len(); ++i) {
            Vec arr{};
            for (int axis = 0; axis < NDIM; ++axis) arr[axis] = x(i, axis);
            total += V(arr);
        }
        return total;
    }

    // Potential gradient
    virtual VecArray gradV(const VecArray& x)
    {
        VecArray result(x.len());
        for (int i = 0; i < x.len(); ++i) {
            Vec arr{};
            for (int axis = 0; axis < NDIM; ++axis) arr[axis] = x(i, axis);
            auto grad = gradV(arr);
            for (int axis = 0; axis < NDIM; ++axis) result(i, axis) = grad[axis];
        }
        return result;
    }

    // Potential laplacian
    virtual double laplacianV(const VecArray& x)
    {
        double total = 0.0;
        for (int i = 0; i < x.len(); ++i) {
            Vec arr{};
            for (int axis = 0; axis < NDIM; ++axis) arr[axis] = x(i, axis);
            total += laplacianV(arr);
        }
        return total;
    }

    // Single-particle overloads (avoid heap allocations)
    virtual double V(const Vec& /* x */) { return 0.0; }
    virtual Vec gradV(const Vec& /* x */) { return {}; }
    virtual double laplacianV(const Vec& /* x */) { return 0.0; }

    [[nodiscard]] virtual bool isFree() const { return true; }

    // Tail correction
    /// @todo Implement the tail correction
    double tail_correction;
};
