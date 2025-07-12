#pragma once

#include <vector>
#include "mpi.h"

#include "common.h"

/**
 * @class NormalModes
 * @brief Handles the transformation and management of system coordinates and momenta between Cartesian and normal mode representations, with support for distributed memory via MPI.
 *
 * This class provides methods to convert coordinates and momenta between Cartesian and normal mode representations.
 *
 * The transformation matrices are stored in row-major order in @ref cart_to_nm_mat_row and @ref nm_to_cart_mat_row. The class also manages MPI windows for efficient parallel access to coordinate and momentum arrays.
 *
 * @note The class assumes the context (see @ref NormalModesContext) provides all necessary information about the system size and MPI environment.
 */
class NormalModes
{
public:
    /**
     * @brief Constructs the NormalModes object and allocates necessary resources.
     * @param coord The shared pointer to the Cartesian coordinates of the system.
     * @param momenta The shared pointer to the momenta of the atoms.
     * @param natoms The number of atoms in the system.
     * @param nbeads The number of imaginary time slices.
     * @param this_bead The index of the current MPI bead (process).
     */
    explicit NormalModes(
        const std::shared_ptr<const dVec>& coord, 
        const std::shared_ptr<dVec>& momenta,
        int natoms,
        int nbeads,
        int this_bead
    );

    /**
     * @brief Destructor. Releases allocated resources and MPI windows.
     */
    ~NormalModes();

    /**
     * @brief Shared memory arrays for coordinates and momenta in both representations.
     *
     * - @ref arr_coord_cartesian: Cartesian coordinates of the whole system.
     * - @ref arr_coord_nm: Normal mode coordinates of the whole system.
     * - @ref arr_momenta_cartesian: Cartesian momenta of the whole system.
     * - @ref arr_momenta_nm: Normal mode momenta of the whole system.
     */
    double *arr_coord_cartesian, *arr_coord_nm, *arr_momenta_cartesian, *arr_momenta_nm;

    /**
     * @brief Row-major storage for transformation matrices.
     *
     * - @ref cart_to_nm_mat_row: Cartesian-to-normal mode transformation matrix.
     * - @ref nm_to_cart_mat_row: Normal mode-to-Cartesian transformation matrix.
     */
    std::vector<double> cart_to_nm_mat_row, nm_to_cart_mat_row;

    /**
     * @brief Shares the coordinate and momentum data across MPI processes.
     */
    void shareData() const;

    /**
     * @brief Converts a Cartesian coordinate to its normal mode representation.
     * @param glob_idx Global index of the coordinate.
     * @return The corresponding normal mode coordinate.
     */
    [[nodiscard]] double coordCartesianToNormal(int glob_idx) const;

    /**
     * @brief Converts a Cartesian momentum to its normal mode representation.
     * @param glob_idx Global index of the momentum.
     * @return The corresponding normal mode momentum.
     */
    [[nodiscard]] double momentumCartesianToNormal(int glob_idx) const;

    /**
     * @brief Converts a normal mode coordinate to its Cartesian representation.
     * @param glob_idx Global index of the normal mode coordinate.
     * @return The corresponding Cartesian coordinate.
     */
    [[nodiscard]] double coordNormalToCartesian(int glob_idx) const;

    /**
     * @brief Converts a normal mode momentum to its Cartesian representation.
     * @param glob_idx Global index of the normal mode momentum.
     * @return The corresponding Cartesian momentum.
     */
    [[nodiscard]] double momentumNormalToCartesian(int glob_idx) const;

    /**
     * @brief Computes the global index for a given atom and axis.
     * @param axis The axis index (e.g., 0 for x, 1 for y, 2 for z).
     * @param atom The atom index.
     * @return The global index in the array.
     */
    [[nodiscard]] int globIndexAtom(const int axis, int atom) const
    {
        return axis * m_axis_stride + atom * m_atom_stride;
    }

    /**
     * @brief Updates the Cartesian momenta from the current normal mode momenta.
     */
    void updateCartesianMomenta() const;

private:
    /**
     * @brief Structure for managing shared memory arrays and their associated MPI windows.
     */
    struct SharedMemory
    {
        double* array; ///< Pointer to the shared array.
        MPI_Win window; ///< MPI window for the array.
    };

    /**
     * @brief Allocates shared memory for a given array and bead.
     * @param mem The shared memory structure to initialize.
     * @param size The size of the array.
     * @param this_bead The MPI bead index.
     */
    static void allocateSharedMemory(SharedMemory& mem, size_t size, int this_bead);

    /**
     * @brief Allocates all shared memory arrays for the current bead.
     * @param this_bead The MPI bead index.
     */
    void allocateAllSharedMemory(int this_bead);

    int m_axis_stride; ///< Stride for axis indexing.
    int m_atom_stride; ///< Stride for atom indexing.

    /// MPI windows for the coordinate and momentum arrays.
    MPI_Win m_win_coord_cartesian, m_win_coord_nm, m_win_momenta_cartesian, m_win_momenta_nm;

    /**
     * @brief Computes the global index for a given atom, axis, and bead.
     * @param axis The axis index.
     * @param atom The atom index.
     * @param bead The bead index.
     * @return The global index in the array.
     */
    [[nodiscard]] int globIndexBead(int axis, int atom, int bead) const { return globIndexAtom(axis, atom) + bead; };

    std::shared_ptr<const dVec> m_coord;
    std::shared_ptr<dVec> m_momenta;
    int m_natoms;
    int m_nbeads;
    int m_this_bead;
};
