#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <stdexcept>
#include <format>
#include <filesystem>
#include <variant>
#include <cmath>
#include <fstream>

#include "units.h"

#ifndef NDIM
#define NDIM 1                       // Number of spatial dimensions
#endif

#ifndef PROGRESS
#define PROGRESS false               // Display progress bar?
#endif

#ifndef MINIM
#define MINIM false                  // Apply minimum image convention when PBC are used
#endif

#ifndef WRAP
#define WRAP false                   // Wrap coordinates when PBC are used
#endif

#ifndef RECENTER
#define RECENTER false               // Recenter polymer when the initial bead moves outside of the fundamental cell?
#endif

#ifndef OLD_BOSONIC_ALGORITHM
#define OLD_BOSONIC_ALGORITHM false  // Enable the old bosonic algorithm that scales as O(N!)?
#endif

// In the "i-Pi convention" [J. Chem. Phys. 133, 124104 (2010); also J. Chem. Phys. 74, 4078–4095 (1981)], 
// the Boltzmann exponents have the form exp[-(beta/P)H], where H is the classical Hamiltonian of the 
// ring polymers. This results in a canonical distribution at P times the physical temperature.
// In constrast, "Tuckerman's convention" [J. Chem. Phys. 99, 2796–2808 (1993)] uses weights of the form exp(-beta*H),
// such that the temperature of the canonical ensemble coincides with the physical temperature.
// Notably, the classical Hamiltonians of the two conventions differ, with the spring constant
// in the i-Pi convention being P times larger than that in Tuckerman's convention. Additionally, the i-Pi convention
// lacks a 1/P prefactor in front of the external potential. The Hamiltonians of the two conventions are related through
// H_tuckerman = H_ipi / P. Note however that the expressions for the various estimators are unaffected by this choice.
// Setting the following pre-processor directive to false amounts to adopting Tuckerman's convention.
#ifndef IPI_CONVENTION
#define IPI_CONVENTION true
#endif

// Progress bar parameters
constexpr int PBWIDTH = 30;
constexpr std::string_view PBSTR = "||||||||||||||||||||||||||||||";

inline const std::string LOGO = R"(
 __       __      __  
|__)||\/||  \ __ |__) 
|   ||  ||__/    |__)
)";

namespace ErrorMessage {
    const std::string GENERAL_ERR = "Error";
    const std::string OVERFLOW_ERR = "Overflow error";
    const std::string INVALID_ARG_ERR = "Invalid argument error";
}

const std::string OUTPUT_FOLDER_NAME = "output";

template<typename T, typename...VariantParams>
void getVariant(const std::variant<VariantParams...>& v, T& value)
{
    value = std::get<T>(v);
}

template <typename T, int dim>
class VectorArray {
public:
    // Initialize an array of vectors.
    // "rows" is the number of vectors contained in the array.
    // Each vector has "dim" components.
    VectorArray() : m_rows(1), m_arr(m_rows * dim, T()) {}
    VectorArray(int rows) : m_rows(rows), m_arr(rows * dim, T()) {}

    // Return the flattened (absolute) index of the "axis" component of the ith vector
    int index(int i, int axis) const { return i * dim + axis; }

    // Number of vectors in the array
    int len() const { return m_rows; }

    // Total number of elements in the underlying raw 1D array
    int size() const { return m_rows * dim; }

    double norm(int vector_idx) const { 
        double vector_norm = 0.0;
        for (int axis_idx = 0; axis_idx < dim; axis_idx++)
            vector_norm += m_arr[index(vector_idx, axis_idx)] * m_arr[index(vector_idx, axis_idx)];
        return sqrt(vector_norm);
    }

    // Return a reference to the underlying 1D array
    std::vector<T>& raw() { return m_arr; }

    // A pointer to the first element in the internal 1D array
    T* data() { return m_arr.data(); }

    // Access the "axis" component of the ith vector
    T& operator()(int i, int axis) { return m_arr[index(i, axis)]; }

    // Access the "axis" component of the ith vector (const version)
    const T& operator()(int i, int axis) const { return m_arr[index(i, axis)]; }

    // Return value based on the flattened index
    T& operator[](size_t idx) { return m_arr[idx]; }

    // Return value based on the flattened index (const version)
    const T& operator[](size_t idx) const { return m_arr[idx]; }

    VectorArray<T, dim> operator*(const T& rhs_scalar) const {
        VectorArray<T, dim> result(*this);  // Create a copy of the current VectorArray

        // Perform scalar multiplication for each element
        for (int i = 0; i < size(); ++i) {
            result.raw()[i] *= rhs_scalar;
        }

        return result;
    }

    VectorArray& operator*=(const T& scalar) {
        for (int i = 0; i < size(); ++i) {
            m_arr[i] *= scalar;
        }
        return *this;
    }

private:
    int m_rows;
    std::vector<T> m_arr;
};

// Use a non-member operator overload for the right-hand side case
template<typename T, int dim>
VectorArray<T, dim> operator*(const T& lhs_scalar, VectorArray<T, dim> rhs_vec) {
    return rhs_vec * lhs_scalar;
}

// Define an array of vectors of doubles of dimension NDIM
using dVec = VectorArray<double, NDIM>;

// Print a general message on "out_rank" (by default, the root rank is 0)
void printMsg(const std::string& msg, int this_rank, int out_rank = 0);

// Print an info message on "out_rank" (by default, the root rank is 0)
void printInfo(const std::string& info, bool& info_flag, int this_rank, int out_rank = 0);

// Print a status message on "out_rank" (by default, the root rank is 0)
void printStatus(const std::string& status, int this_rank, int out_rank = 0);

// Print an error message on "out_rank" (by default, the root rank is 0)
void printError(const std::string& msg, int this_rank, const std::string& err_type = std::string(), int out_rank = 0);

// Print a progress bar
void printProgress(int this_step, int total_steps, int this_rank, int out_rank = 0);

// To handle periodic boundary conditions, we employ the Class C storage
// concept [Z. Phys. Chem. 227 (2013) 345–352], allowing atoms to move
// outside the primary simulation box. That is, the coordinates are not folded
// into the simulation box. Instead, we account for the PBC when calculating the
// distances between particles (or any spatial vector differences).
// This is done using Algorithm C4. It calculates the remainded of dx
// on the interval [-L/2, L/2].
void applyMinimumImage(double& dx, double L);
//void applyMinimumImage(dVec& dx_arr, double L);

void periodicWrap(double& x, double L);
//void periodicWrap(dVec& pos_arr, double L);

// Load positions from an .xyz file to the destination vector
void loadXYZ(const std::string& xyz_filename, dVec& destination);

// Load momenta from a file to the destination vector
void loadMomenta(const std::string& vel_filename, double mass, dVec& destination);

template <typename T>
std::string formattedReportLine(const std::string& property_name, const T& value) {
    return std::format("{:<40}\t:\t{}\n", property_name, value);
}