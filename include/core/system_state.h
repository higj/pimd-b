#pragma once

#include "common.h"

class SystemState {
public:
    dVec coord, momenta;
    dVec spring_forces, physical_forces;
    dVec prev_coord, next_coord;

    /**
     * Initializes the system state (coordinates, momenta and forces) with the given parameters.
     * 
     * @param rank Current process id ("rank" of MPI_Comm_rank), also referred to as "this_bead"
     * @param nproc Number of processes ("size" of MPI_Comm_size)
     * @param natoms Number of atoms in the system
     * @param nbeads Number of beads in the system
     */
    SystemState(int rank, int nproc, int natoms, int nbeads);
    //SystemState(int rank, int nproc, int natoms, int nbeads, double box_size); @param box_size Linear system size (box size)

    /**
     * @brief Updates the neighboring coordinates arrays.
     */
    void updateNeighboringCoordinates();

    /**
     * @brief Zero the linear momentum of a group of atoms by subtracting the velocity
     * of the center of mass from the velocity of each atom.
     * The calculation assumes that all atoms have the same mass, in which case the
     * center of mass momentum is given by p_c=m*v_c=(p_1+...+p_n)/n, where n=N*P
     * is the total number of beads in the system.
     */
    void zeroMomentum();

    /**
     * @brief Get the total number of atoms in the system.
     */
    [[nodiscard]] int getNumAtoms() const { return m_natoms; }

    /**
     * @brief Get the total number of beads in the system.
     */
    [[nodiscard]] int getNumBeads() const { return m_nbeads; }

    /**
     * @brief Get the current bead (imaginary time-slice) index.
     */
    [[nodiscard]] int currentBead() const { return m_rank; }

    /**
     * @brief Get the linear size of the system (box size).
     */
    //[[nodiscard]] double getBoxSize() const { return m_box_size; }

    /**
     * Get the total force acting on a particle along a specific axis.
     *
     * @param ptcl_idx Index of the particle.
     * @param axis Axis along which to get the force.
     * @return Total force acting on the particle along the specified axis.
     */
    [[nodiscard]] double getTotalForce(int ptcl_idx, int axis) const;

    /**
     * Returns the vectorial distance between two particles at the same imaginary time slice.
     *
     * @param first_ptcl Index of the first particle.
     * @param second_ptcl Index of the second particle.
     * @param minimum_image Flag determining whether the minimum image convention should be applied.
     * @return
     */
    //[[nodiscard]] dVec getSeparation(int first_ptcl, int second_ptcl, bool minimum_image) const;

    /**
     * Returns the vectorial distance between two particles at the same imaginary time slice.
     * Mathematically equivalent to r1-r2, where r1 and r2 are the position-vectors
     * of the first and second particles, respectively. In the case of PBC, there is an
     * option to return the minimal distance between the two particles.
     *
     * @param positions Array of particle positions.
     * @param first_ptcl Index of the first particle.
     * @param second_ptcl Index of the second particle.
     * @param minimum_image Flag determining whether the minimum image convention should be applied.
     * @return Vectorial distance between the two particles.
     */
    //[[nodiscard]] dVec getSeparation(const dVec& positions, int first_ptcl, int second_ptcl, bool minimum_image) const;
private:         
    int m_rank; // Current process id ("rank" of MPI_Comm_rank), also referred to as "this_bead"
    int m_nproc; // Number of processes ("size" of MPI_Comm_size)
    /// Useful to have the following constants duplicated here:
    int m_natoms, m_nbeads; // Number of atoms and beads
    //double m_box_size; // Box size (linear system size)

    /**
     * Receives coordinates from the previous time-slice,
     * and sends the current coordinates to the next time-slice.
     *
     * @param prev Vector to store the previous coordinates.
     */
    void updatePreviousCoordinates(dVec& prev);

    /**
     * Receives coordinates from the next time-slice,
     * and sends the current coordinates to the previous time-slice.
     *
     * @param next Vector to store the next coordinates.
     */
    void updateNextCoordinates(dVec& next);
};
