#include "propagators/normal_modes/normal_modes.h"
#include "propagators/normal_modes/normal_modes_transformation_matrix.h"

NormalModes::NormalModes(const NormalModesContext& context) :
    cart_to_nm_mat_row(context.nbeads),
    nm_to_cart_mat_row(context.nbeads),
    axis_stride(context.natoms * context.nbeads),
    atom_stride(context.nbeads),
    m_context(context)
{
    //// Allocate shared memory
    //if (context.this_bead == 0) {
    //    MPI_Win_allocate_shared(sim.natoms*sim.nbeads*NDIM*sizeof(double), sizeof(double),
    //        MPI_INFO_NULL, MPI_COMM_WORLD,
    //        &arr_coord_cartesian, &win_coord_cartesian);
    //    MPI_Win_allocate_shared(sim.natoms*sim.nbeads*NDIM*sizeof(double), sizeof(double),
    //        MPI_INFO_NULL, MPI_COMM_WORLD,
    //        &arr_coord_nm, &win_coord_nm);
    //    MPI_Win_allocate_shared(sim.natoms*sim.nbeads*NDIM*sizeof(double), sizeof(double),
    //        MPI_INFO_NULL, MPI_COMM_WORLD,
    //        &arr_momenta_cartesian, &win_momenta_cartesian);
    //    MPI_Win_allocate_shared(sim.natoms*sim.nbeads*NDIM*sizeof(double), sizeof(double),
    //        MPI_INFO_NULL, MPI_COMM_WORLD,
    //        &arr_momenta_nm, &win_momenta_nm);
    //} else {
    //    int disp_unit;
    //    MPI_Aint size;
    //    MPI_Win_allocate_shared(0, sizeof(double),
    //        MPI_INFO_NULL, MPI_COMM_WORLD,
    //        &arr_coord_cartesian, &win_coord_cartesian);
    //    MPI_Win_allocate_shared(0, sizeof(double),
    //        MPI_INFO_NULL, MPI_COMM_WORLD,
    //        &arr_coord_nm, &win_coord_nm);
    //    MPI_Win_allocate_shared(0, sizeof(double),
    //        MPI_INFO_NULL, MPI_COMM_WORLD,
    //        &arr_momenta_cartesian, &win_momenta_cartesian);
    //    MPI_Win_allocate_shared(0, sizeof(double),
    //        MPI_INFO_NULL, MPI_COMM_WORLD,
    //        &arr_momenta_nm, &win_momenta_nm);
    //    MPI_Win_shared_query(win_coord_cartesian, 0, &size, &disp_unit, &arr_coord_cartesian);
    //    MPI_Win_shared_query(win_coord_nm, 0, &size, &disp_unit, &arr_coord_nm);
    //    MPI_Win_shared_query(win_momenta_cartesian, 0, &size, &disp_unit, &arr_momenta_cartesian);
    //    MPI_Win_shared_query(win_momenta_nm, 0, &size, &disp_unit, &arr_momenta_nm);
    //}
    allocateAllSharedMemory(context.this_bead);
        
    //// Cartesian-to-nm transformation matrix (one row because parallelized)
    //double pref;
    //double fund_freq = 2 * std::numbers::pi / sim.nbeads * sim.this_bead;
    //if (sim.this_bead == 0) {
    //    pref = 1 / sqrt(sim.nbeads);
    //    std::ranges::fill(cart_to_nm_mat_row, pref);
    //} else if (sim.this_bead < 0.5 * sim.nbeads) {
    //    pref = sqrt(2.0 / sim.nbeads);
    //    for (int i = 0; i < sim.nbeads; ++i)
    //        cart_to_nm_mat_row[i] = pref * cos(fund_freq * i);
    //} else if (sim.this_bead == 0.5 * sim.nbeads) {
    //    pref = 1 / sqrt(sim.nbeads);
    //    for (int i = 0; i < sim.nbeads; ++i)
    //        cart_to_nm_mat_row[i] = pref * (i % 2 == 0 ? 1.0 : -1.0);
    //} else {
    //    pref = sqrt(2.0 / sim.nbeads);
    //    for (int i = 0; i < sim.nbeads; ++i)
    //        cart_to_nm_mat_row[i] = -pref * sin(fund_freq * i);
    //}
    //
    //// NM-to-Cartesian transformation matrix row
    //pref = sqrt(2.0 / sim.nbeads);
    //nm_to_cart_mat_row[0] = 1 / sqrt(sim.nbeads);
    //for (int i = 1; i < 0.5 * sim.nbeads; ++i)
    //    nm_to_cart_mat_row[i] = pref * cos(fund_freq * i);
    //if (sim.nbeads % 2 == 0)
    //    nm_to_cart_mat_row[sim.nbeads / 2] = 1 / sqrt(sim.nbeads) * (sim.this_bead % 2 == 0 ? 1.0 : -1.0);
    //for (int i = std::ceil(0.5 * (sim.nbeads + 1)); i < sim.nbeads; ++i)
    //    nm_to_cart_mat_row[i] = -pref * sin(fund_freq * i);

    TransformationMatrixBuilder builder(context.this_bead, context.nbeads);

    // Build both transformation matrices
    builder.buildCartToNM(cart_to_nm_mat_row.data());
    builder.buildNMToCart(nm_to_cart_mat_row.data());
}

NormalModes::~NormalModes() {
    MPI_Win_free(&win_coord_cartesian);
    MPI_Win_free(&win_coord_nm);
    MPI_Win_free(&win_momenta_cartesian);
    MPI_Win_free(&win_momenta_nm);
}

// Function to allocate shared memory for a specific array
void NormalModes::allocateSharedMemory(SharedMemory& mem, const size_t size, const int this_bead) {
    if (this_bead == 0) {
        MPI_Win_allocate_shared(size, sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &mem.array, &mem.window);
    } else {
        MPI_Win_allocate_shared(0, sizeof(double),
            MPI_INFO_NULL, MPI_COMM_WORLD,
            &mem.array, &mem.window);

        MPI_Aint win_size;
        int disp_unit;
        MPI_Win_shared_query(mem.window, 0, &win_size, &disp_unit, &mem.array);
    }
}

// Main allocation function that handles all arrays
void NormalModes::allocateAllSharedMemory(int this_bead) {
    // Calculate total memory size needed for each array
    /// TODO: Fix implicit conversion
    const size_t array_size = m_context.natoms * m_context.nbeads * NDIM * sizeof(double);

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

    // Store the pointers in your global variables (or wherever needed)
    arr_coord_cartesian = coord_cartesian.array;
    win_coord_cartesian = coord_cartesian.window;

    arr_coord_nm = coord_nm.array;
    win_coord_nm = coord_nm.window;

    arr_momenta_cartesian = momenta_cartesian.array;
    win_momenta_cartesian = momenta_cartesian.window;

    arr_momenta_nm = momenta_nm.array;
    win_momenta_nm = momenta_nm.window;
}

void NormalModes::shareData()
{
    for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            int glob_idx = globIndexAtom(axis, ptcl_idx);
            arr_coord_cartesian[glob_idx + m_context.this_bead] = (*m_context.coord)(ptcl_idx, axis);
            arr_momenta_cartesian[glob_idx + m_context.this_bead] = (*m_context.momenta)(ptcl_idx, axis);
        }
    }
}

double NormalModes::coordCartesianToNormal(const int glob_idx) const
{
    double coord_nm = 0;
    for (int bead_idx = 0; bead_idx < m_context.nbeads; ++bead_idx) {
        coord_nm += cart_to_nm_mat_row[bead_idx] * arr_coord_cartesian[glob_idx + bead_idx];
    }
    return coord_nm;
}

double NormalModes::momentumCartesianToNormal(const int glob_idx) const
{
    double momentum_nm = 0;
    for (int bead_idx = 0; bead_idx < m_context.nbeads; ++bead_idx) {
        momentum_nm += cart_to_nm_mat_row[bead_idx] * arr_momenta_cartesian[glob_idx + bead_idx];
    }
    return momentum_nm;
}

double NormalModes::coordNormalToCartesian(const int glob_idx) const
{
    double coord_cartesian = 0;
    for (int bead_idx = 0; bead_idx < m_context.nbeads; ++bead_idx) {
        coord_cartesian += nm_to_cart_mat_row[bead_idx] * arr_coord_nm[glob_idx + bead_idx];
    }
    return coord_cartesian;
}

double NormalModes::momentumNormalToCartesian(const int glob_idx) const
{
    double momentum_cartesian = 0;
    for (int bead_idx = 0; bead_idx < m_context.nbeads; ++bead_idx) {
        momentum_cartesian += nm_to_cart_mat_row[bead_idx] * arr_momenta_nm[glob_idx + bead_idx];
    }
    return momentum_cartesian;
}

void NormalModes::updateCartesianMomenta() {
    for (int ptcl_idx = 0; ptcl_idx < m_context.natoms; ++ptcl_idx) {
        for (int axis = 0; axis < NDIM; ++axis) {
            int glob_idx = globIndexAtom(axis, ptcl_idx);
            arr_momenta_cartesian[glob_idx + m_context.this_bead] = momentumNormalToCartesian(glob_idx);
            (*m_context.momenta)(ptcl_idx, axis) = arr_momenta_cartesian[glob_idx + m_context.this_bead];
        }
    }
}

