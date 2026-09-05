#pragma once

#include <vector>
#include <array>
#include <format>
#include <variant>
#include <cmath>
#include "ordered_map.h"

#include "units.h"

#ifndef NDIM
#define NDIM 1                       // Number of spatial dimensions
#endif

#ifndef MINIM
#define MINIM true                  // Apply minimum image convention when PBC are used
#endif

#ifndef WRAP
#define WRAP true                   // Wrap coordinates when PBC are used
#endif

#ifndef FACTORIAL_BOSONIC_ALGORITHM
#define FACTORIAL_BOSONIC_ALGORITHM false  // Enable the old bosonic algorithm that scales as O(N!)?
#endif

// In the "Tau convention" [J. Chem. Phys. 133, 124104 (2010); also J. Chem. Phys. 74, 4078-4095 (1981)], 
// the Boltzmann exponents have the form exp[-(beta/P)*H_tau], where H_tau is the classical Hamiltonian of the 
// ring polymers. This results in a canonical distribution at P times the physical temperature.
// In contrast, "Beta convention" [J. Chem. Phys. 99, 2796-2808 (1993)] uses weights of the form exp(-beta*H_beta),
// such that the temperature of the canonical ensemble coincides with the physical temperature.
// Notably, the classical Hamiltonians of the two conventions differ, with the spring constant
// in the tau convention being P times larger than that in beta convention. Additionally, the tau convention
// lacks a 1/P prefactor in front of the external potential. The Hamiltonians of the two conventions are related through
// H_beta = H_tau / P. Note however that the expressions for the various estimators are unaffected by this choice.
// Setting the following pre-processor directive to false amounts to adopting the beta convention.
#ifndef TAU_CONVENTION
#define TAU_CONVENTION true
#endif

// A small number (but not necessarily the smallest)
constexpr auto EPS = 1.0E-7;

/**
 * A class to store an array of vectors of dimension "dim".
 *
 * @tparam T Type of the elements in the vectors.
 * @tparam dim Dimension of the vectors.
 */
template <typename T, int dim>
class VectorCollection {
public:
    VectorCollection() : m_rows(1), m_arr(m_rows * dim, T()) {}
    explicit VectorCollection(const int rows) : m_rows(rows), m_arr(rows * dim, T()) {}

    /**
     * Retrieve the flattened (absolute) index of the "axis" component of the ith vector.
     * Uses a Structure of Arrays (SoA) memory layout (e.g., [X0, X1... Y0, Y1...]) 
     * instead of an Array of Structures (AoS) layout ([X0, Y0, Z0, X1, Y1, Z1...]).
     * This layout allows the compiler to easily auto-vectorize loops (SIMD) 
     * when computing pairwise forces over coordinates.
     *
     * @param i Vector index.
     * @param axis Axis index.
     * @return Index of the element in the underlying one-dimensional array.
     */
    [[nodiscard]] int index(int i, int axis) const {
        return axis * m_rows + i;
    }

    /**
     * Number of vectors in the array.
     *
     * @return Number of vectors in the array.
     */
    [[nodiscard]] int len() const {
        return m_rows;
    }

    /**
     * Total number of elements in the underlying one-dimensional array.
     *
     * @return Total number of elements in the underlying array.
     */
    [[nodiscard]] int size() const {
        return m_rows * dim;
    }

    /**
     * Calculates the Euclidean norm of the vector at the given index.
     * If no index is provided, the norm of the first vector is calculated.
     *
     * @param vector_idx Location of the vector in the array.
     * @return Norm of the vector.
     */
    [[nodiscard]] double norm(int vector_idx = 0) const {
        double vector_norm = 0.0;

        for (int axis_idx = 0; axis_idx < dim; axis_idx++)
            vector_norm += m_arr[index(vector_idx, axis_idx)] * m_arr[index(vector_idx, axis_idx)];
        return sqrt(vector_norm);
    }

    /**
     * Returns a reference to the underlying one-dimensional array.
     *
     * @return Reference to the underlying array.
     */
    std::vector<T>& raw() {
        return m_arr;
    }

    /**
     * A pointer to the first element in the internal one-dimensional array.
     *
     * @return A pointer to the first element in the underlying array.
     */
    T* data() {
        return m_arr.data();
    }

    /**
     * Access the "axis" component of the ith vector (modifiable).
     *
     * @param i Vector index.
     * @param axis Axis index.
     * @return Reference to the element at the specified location.
     */
    T& operator()(int i, int axis) {
        return m_arr[index(i, axis)];
    }

    /**
     * Access the "axis" component of the ith vector (const version).
     *
     * @param i Vector index.
     * @param axis Axis index.
     * @return Reference to the element at the specified location.
     */
    const T& operator()(int i, int axis) const {
        return m_arr[index(i, axis)];
    }

    /**
     * Retrieve value based on the flattened index (modifiable).
     *
     * @param idx Index in the flattened array.
     * @return Value at the specified index.
     */
    T& operator[](size_t idx) {
        return m_arr[idx];
    }

    /**
     * Retrieve value based on the flattened index (const version).
     *
     * @param idx Index in the flattened array.
     * @return Vale at the specified index.
     */
    const T& operator[](size_t idx) const {
        return m_arr[idx];
    }

    /**
     * Overload the addition operator for vector addition.
     *
     * @param other Vector array to be added.
     * @return New vector resulting from the addition.
     */
    VectorCollection<T, dim> operator+(const VectorCollection<T, dim>& other) const {
        VectorCollection<T, dim> result(*this);  // Create a copy of the current VectorCollection

        // Perform element-wise addition
        for (int i = 0; i < size(); ++i) {
            result.raw()[i] += other[i];
        }

        return result;
    }

    /**
     * Overload the multiplication operator for scalar multiplication (from the right-hand side).
     *
     * @param rhs_scalar Scalar value to multiply the vector by.
     * @return Copy of the vector multiplied by the scalar.
     */
    VectorCollection<T, dim> operator*(const T& rhs_scalar) const {
        VectorCollection<T, dim> result(*this);  // Create a copy of the current VectorCollection

        // Perform scalar multiplication for each element
        for (int i = 0; i < size(); ++i) {
            result.raw()[i] *= rhs_scalar;
        }

        return result;
    }

    /**
     * Overload the *= operator for scalar multiplication of the vector.
     *
     * @param scalar Scalar value to multiply the vector by.
     * @return Reference to the modified vector.
     */
    VectorCollection& operator*=(const T& scalar) {
        for (int i = 0; i < size(); ++i) {
            m_arr[i] *= scalar;
        }
        return *this;
    }

    /**
     * Calculate the difference between two vectors in the array (v_i-v_j).
     *
     * @param i First vector index.
     * @param j Second vector index.
     * @return VectorCollection containing the difference between the two vectors.
     */
    VectorCollection<T, dim> getSeparation(int i, int j) const {
        VectorCollection<T, dim> sep;

        for (int axis = 0; axis < dim; ++axis) {
            sep(0, axis) = (*this)(i, axis) - (*this)(j, axis);
        }

        return sep;
    }

    /**
     * Calculate the difference between two vectors in the array (v_i-v_j).
     * Returns a lightweight std::array to avoid heap allocations in inner loops.
     *
     * @param i First vector index.
     * @param j Second vector index.
     * @return std::array containing the difference between the two vectors.
     */
    std::array<T, dim> getSeparationArray(int i, int j) const {
        std::array<T, dim> sep{};
        for (int axis = 0; axis < dim; ++axis) {
            sep[axis] = (*this)(i, axis) - (*this)(j, axis);
        }
        return sep;
    }

    /**
     * Reset all values in the vector array to zero (or default value of type T).
     */
    void reset() {
        std::fill(m_arr.begin(), m_arr.end(), T());
    }

private:
    int m_rows; // Number of vectors in the array
    std::vector<T> m_arr; // 1D array to store the vectors
};

// Use a non-member operator overload for the right-hand side case
template <typename T, int dim>
VectorCollection<T, dim> operator*(const T& lhs_scalar, VectorCollection<T, dim> rhs_vec) {
    return rhs_vec * lhs_scalar;
}

// Define an array of vectors of doubles of dimension NDIM
using VecArray = VectorCollection<double, NDIM>;

// Define a lightweight array for single-particle / pairwise vectors of dimension NDIM
using Vec = std::array<double, NDIM>;

// Helper function to calculate the norm of an Vec
inline double norm(const Vec& v) {
    double sum = 0.0;
    for (double val : v) sum += val * val;
    return std::sqrt(sum);
}

// Define an array of vectors of integers of dimension NDIM
using IntVecArray = VectorCollection<int, NDIM>;

// Define a map of variant types
using VariantMap = std::unordered_map<std::string, std::variant<int, unsigned int, long, double, bool, std::string>>;

// Define a map of strings
using StringMap = tsl::ordered_map<std::string, std::string>;

// Print a general message on "out_rank" (by default, the root rank is 0)
void printMsg(const std::string& msg, int this_rank, int out_rank = 0);

// Print a warning message on "out_rank" (by default, the root rank is 0)
void printWarning(const std::string& msg, int this_rank, int out_rank = 0);

// Print an info message on "out_rank" (by default, the root rank is 0)
//void printInfo(const std::string& info, bool& info_flag, int this_rank, int out_rank = 0);

// Print a status message on "out_rank" (by default, the root rank is 0)
void printStatus(const std::string& status, int this_rank, int out_rank = 0);

// Print an error message on "out_rank" (by default, the root rank is 0)
void printError(const std::string& msg, int this_rank, const std::string& err_type = std::string(), int out_rank = 0);