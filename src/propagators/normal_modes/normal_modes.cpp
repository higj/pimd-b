#include "propagators/normal_modes/normal_modes.h"
#include "propagators/normal_modes/normal_modes_transformation_matrix.h"

NormalModes::NormalModes(
    const std::shared_ptr<const dVec>& coord,
    const std::shared_ptr<dVec>& momenta,
    int natoms,
    int nbeads,
    int this_bead
) :
    cart_to_nm_mat_row(nbeads),
    nm_to_cart_mat_row(nbeads),
    m_axis_stride(natoms * nbeads),
    m_atom_stride(nbeads),
    m_coord(coord),
    m_momenta(momenta),
    m_natoms(natoms),
    m_nbeads(nbeads),
    m_this_bead(this_bead)
{
    allocateAllSharedMemory(this_bead);

    const TransformationMatrixBuilder builder(this_bead, nbeads);

    // Build both transformation matrices
    builder.buildCartesianToNormalModes(cart_to_nm_mat_row.data());
    builder.buildNormalModesToCartesian(nm_to_cart_mat_row.data());
}

NormalModes::~NormalModes()
{
    MPI_Win_free(&m_win_coord_cartesian);
    MPI_Win_free(&m_win_coord_nm);
    MPI_Win_free(&m_win_momenta_cartesian);
    MPI_Win_free(&m_win_momenta_nm);
}

// Function to allocate shared memory for a specific array
void NormalModes::allocateSharedMemory(SharedMemory& mem, const size_t size, const int this_bead)
{
    if (this_bead == 0)
    {
        MPI_Win_allocate_shared(size, sizeof(double),
                                MPI_INFO_NULL, MPI_COMM_WORLD,
                                &mem.array, &mem.window);
    }
    else
    {
        MPI_Win_allocate_shared(0, sizeof(double),
                                MPI_INFO_NULL, MPI_COMM_WORLD,
                                &mem.array, &mem.window);

        MPI_Aint win_size;
        int disp_unit;
        MPI_Win_shared_query(mem.window, 0, &win_size, &disp_unit, &mem.array);
    }
}

// Main allocation function that handles all arrays
void NormalModes::allocateAllSharedMemory(int this_bead)
{
    // Calculate total memory size needed for each array
    const size_t array_size = static_cast<size_t>(m_natoms * m_nbeads * NDIM) * sizeof(double);

    // Define all shared memory structures
    SharedMemory coord_cartesian;
    SharedMemory coord_nm;
    SharedMemory momenta_cartesian;
    SharedMemory momenta_nm;

    // Allocate shared memory for all arrays
    allocateSharedMemory(coord_cartesian, array_size, this_bead);
    allocateSharedMemory(coord_nm, array_size, this_bead);
    allocateSharedMemory(momenta_cartesian, array_size, this_bead);
    allocateSharedMemory(momenta_nm, array_size, this_bead);

    arr_coord_cartesian = coord_cartesian.array;
    m_win_coord_cartesian = coord_cartesian.window;

    arr_coord_nm = coord_nm.array;
    m_win_coord_nm = coord_nm.window;

    arr_momenta_cartesian = momenta_cartesian.array;
    m_win_momenta_cartesian = momenta_cartesian.window;

    arr_momenta_nm = momenta_nm.array;
    m_win_momenta_nm = momenta_nm.window;
}

void NormalModes::shareData() const
{
    for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const int glob_idx = globIndexAtom(axis, ptcl_idx);
            arr_coord_cartesian[glob_idx + m_this_bead] = (*m_coord)(ptcl_idx, axis);
            arr_momenta_cartesian[glob_idx + m_this_bead] = (*m_momenta)(ptcl_idx, axis);
        }
    }
}

double NormalModes::coordCartesianToNormal(const int glob_idx) const
{
    double coord_nm = 0;
    for (int bead_idx = 0; bead_idx < m_nbeads; ++bead_idx)
    {
        coord_nm += cart_to_nm_mat_row[bead_idx] * arr_coord_cartesian[glob_idx + bead_idx];
    }
    return coord_nm;
}

double NormalModes::momentumCartesianToNormal(const int glob_idx) const
{
    double momentum_nm = 0;
    for (int bead_idx = 0; bead_idx < m_nbeads; ++bead_idx)
    {
        momentum_nm += cart_to_nm_mat_row[bead_idx] * arr_momenta_cartesian[glob_idx + bead_idx];
    }
    return momentum_nm;
}

double NormalModes::coordNormalToCartesian(const int glob_idx) const
{
    double coord_cartesian = 0;
    for (int bead_idx = 0; bead_idx < m_nbeads; ++bead_idx)
    {
        coord_cartesian += nm_to_cart_mat_row[bead_idx] * arr_coord_nm[glob_idx + bead_idx];
    }
    return coord_cartesian;
}

double NormalModes::momentumNormalToCartesian(const int glob_idx) const
{
    double momentum_cartesian = 0;
    for (int bead_idx = 0; bead_idx < m_nbeads; ++bead_idx)
    {
        momentum_cartesian += nm_to_cart_mat_row[bead_idx] * arr_momenta_nm[glob_idx + bead_idx];
    }
    return momentum_cartesian;
}

void NormalModes::updateCartesianMomenta() const
{
    for (int ptcl_idx = 0; ptcl_idx < m_natoms; ++ptcl_idx)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            int glob_idx = globIndexAtom(axis, ptcl_idx);
            arr_momenta_cartesian[glob_idx + m_this_bead] = momentumNormalToCartesian(glob_idx);
            (*m_momenta)(ptcl_idx, axis) = arr_momenta_cartesian[glob_idx + m_this_bead];
        }
    }
}
