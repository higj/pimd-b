#include "bosonic_exchange.h"
#include "mpi.h"

BosonicExchange::BosonicExchange(int nbosons, int np, int bead_num, double _beta, double spring_constant, 
    const dVec x, const dVec x_prev, const dVec x_next, bool pbc, double L) : 
    BosonicExchangeBase(nbosons, np, bead_num, _beta, spring_constant, x, x_prev, x_next, pbc, L),
    temp_nbosons_array(nbosons),
    separate_atom_spring(nbosons),
    E_kn(int(nbosons* (nbosons + 1) / 2)),
    V(nbosons + 1),
    V_backwards(nbosons + 1),
    connection_probabilities(nbosons* nbosons),
    prim_est(nbosons + 1)
{
    prepare_with_coordinates();
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::prepare_with_coordinates()
{
    evaluate_cycle_energies();
    if (bead_num == 0 || bead_num == np - 1) {
        // Exterior beads
        Evaluate_VBn();
        Evaluate_V_backwards();
        evaluate_connection_probabilities();
    }
}

/* ---------------------------------------------------------------------- */

BosonicExchange::~BosonicExchange()
{
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::updateCoordinates(const dVec new_x, const dVec new_x_prev, const dVec new_x_next)
{
    BosonicExchangeBase::updateCoordinates(new_x, new_x_prev, new_x_next);
    prepare_with_coordinates();
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::evaluate_cycle_energies()
{
    for (int i = 0; i < nbosons; i++) {
        temp_nbosons_array[i] = distance_squared_two_beads(x, i, x_next, i);
    }
    
    // Reduce the result and send to bead_num=0
    MPI_Reduce(temp_nbosons_array.data(),
        separate_atom_spring.data(),
        nbosons,
        MPI_DOUBLE,
        MPI_SUM,
        0,
        MPI_COMM_WORLD
    );

    if (bead_num == 0 || bead_num == np - 1) {
        dVec x_first_bead;
        dVec x_last_bead;

        if (bead_num == 0) {
            // Send to bead_num=np-1
            MPI_Send(separate_atom_spring.data(), nbosons, MPI_DOUBLE, np - 1, 0, MPI_COMM_WORLD);

            x_first_bead = x;
            x_last_bead = x_prev;
        }
        else {
            // Receive at bead_num=np-1 from bead_num=0
            MPI_Recv(separate_atom_spring.data(), nbosons, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            x_first_bead = x_next;
            x_last_bead = x;
        }

        for (int v = 0; v < nbosons; v++) {
            set_Enk(v + 1, 1,
                0.5 * spring_constant * separate_atom_spring[v]);

            for (int u = v - 1; u >= 0; u--) {
                double val = get_Enk(v + 1, v - u) +
                    0.5 * spring_constant * (
                        // Eint(u)
                        separate_atom_spring[u] - distance_squared_two_beads(x_first_bead, u, x_last_bead, u)
                        // connect u to u+1
                        + distance_squared_two_beads(x_last_bead, u, x_first_bead, u + 1)
                        // break cycle [u+1,v]
                        - distance_squared_two_beads(x_first_bead, u + 1, x_last_bead, v)
                        // close cycle from v to u
                        + distance_squared_two_beads(x_first_bead, u, x_last_bead, v));

                set_Enk(v + 1, v - u + 1, val);
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

double BosonicExchange::get_Enk(int m, int k)
{
    int end_of_m = m * (m + 1) / 2;
    return E_kn[end_of_m - k];
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::set_Enk(int m, int k, double val)
{
    int end_of_m = m * (m + 1) / 2;
    E_kn[end_of_m - k] = val;
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::Evaluate_VBn()
{
    V[0] = 0.0;

    for (int m = 1; m < nbosons + 1; m++) {
        double Elongest = std::numeric_limits<double>::max();

        for (int k = m; k > 0; k--) {
            double val = get_Enk(m, k) + V[m - k];
            Elongest = std::min(Elongest, val);
            temp_nbosons_array[k - 1] = val;
        }

        double sig_denom = 0.0;
        for (int k = m; k > 0; k--) {
            sig_denom += exp(-beta * (temp_nbosons_array[k - 1] - Elongest));
        }

        V[m] = Elongest - (1.0 / beta) * log(sig_denom / (double)m);

        if (!std::isfinite(V[m])) {
            throw std::overflow_error(
                std::format("Invalid sig_denom {:4.2f} with Elongest {:4.2f} in bosonic exchange potential", sig_denom, Elongest)
            );
        }
    }
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::Evaluate_V_backwards()
{
    V_backwards[nbosons] = 0.0;

    for (int l = nbosons - 1; l > 0; l--) {
        double Elongest = std::numeric_limits<double>::max();
        for (int p = l; p < nbosons; p++) {
            double val = get_Enk(p + 1, p - l + 1) + V_backwards[p + 1];
            Elongest = std::min(Elongest, val);
            temp_nbosons_array[p] = val;
        }

        double sig_denom = 0.0;
        for (int p = l; p < nbosons; p++) {
            sig_denom += 1.0 / (p + 1) * exp(-beta *
                (temp_nbosons_array[p]
                    - Elongest)
            );
        }

        V_backwards[l] = Elongest - log(sig_denom) / beta;

        if (!std::isfinite(V_backwards[l])) {
            throw std::overflow_error(
                std::format("Invalid sig_denom {:4.2f} with Elongest {:4.2f} in bosonic exchange potential", sig_denom, Elongest)
            );
        }
    }

    V_backwards[0] = V[nbosons];
}

/* ---------------------------------------------------------------------- */

double BosonicExchange::get_potential() const
{
    return V[nbosons];
}

/* ---------------------------------------------------------------------- */

double BosonicExchange::get_Vn(int n) const
{
    return V[n];
}

/* ---------------------------------------------------------------------- */

double BosonicExchange::get_E_kn_serial_order(int i) const
{
    return E_kn[i];
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::evaluate_connection_probabilities()
{
    for (int l = 0; l < nbosons - 1; l++) {
        double direct_link_probability = 1.0 - (exp(-beta *
            (V[l + 1] + V_backwards[l + 1] -
                V[nbosons])));
        connection_probabilities[nbosons * l + (l + 1)] = direct_link_probability;
    }
    for (int u = 0; u < nbosons; u++) {
        for (int l = u; l < nbosons; l++) {
            double close_cycle_probability = 1.0 / (l + 1) *
                exp(-beta * (V[u] + get_Enk(l + 1, l - u + 1) + V_backwards[l + 1]
                    - V[nbosons]));
            connection_probabilities[nbosons * l + u] = close_cycle_probability;
        }
    }
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::spring_force_last_bead(dVec& f)
{
    const dVec x_first_bead = x_next;
    const dVec x_last_bead = x;

    for (int l = 0; l < nbosons; l++) {
        std::vector<double> sums(NDIM, 0.0);

        for (int next_l = 0; next_l <= l + 1 && next_l < nbosons; next_l++) {
            double diff_next[NDIM];

            diff_two_beads(x_last_bead, l, x_first_bead, next_l, diff_next);

            double prob = connection_probabilities[nbosons * l + next_l];

            for (int axis = 0; axis < NDIM; ++axis) {
                sums[axis] += prob * diff_next[axis];
            }
        }

        double diff_prev[NDIM];
        diff_two_beads(x_last_bead, l, x_prev, l, diff_prev);

        for (int axis = 0; axis < NDIM; ++axis) {
            sums[axis] += diff_prev[axis];
        }

        for (int axis = 0; axis < NDIM; ++axis) {
            f(l, axis) += sums[axis] * spring_constant;
        }
    }
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::spring_force_first_bead(dVec& f)
{
    const dVec x_first_bead = x;
    const dVec x_last_bead = x_prev;

    for (int l = 0; l < nbosons; l++) {
        std::vector<double> sums(NDIM, 0.0);

        for (int prev_l = std::max(0, l - 1); prev_l < nbosons; prev_l++) {
            double diff_prev[NDIM];

            diff_two_beads(x_first_bead, l, x_last_bead, prev_l, diff_prev);

            double prob = connection_probabilities[nbosons * prev_l + l];

            for (int axis = 0; axis < NDIM; ++axis) {
                sums[axis] += prob * diff_prev[axis];
            }
        }

        double diff_next[NDIM];
        diff_two_beads(x_first_bead, l, x_next, l, diff_next);

        for (int axis = 0; axis < NDIM; ++axis) {
            sums[axis] += diff_next[axis];
        }

        for (int axis = 0; axis < NDIM; ++axis) {
            f(l, axis) += sums[axis] * spring_constant;
        }
    }
}

/* ---------------------------------------------------------------------- */

// Primitive kinetic energy estimator for bosons.
// Corresponds to Eqns. (4)-(5) in SI of pnas.1913365116
double BosonicExchange::prim_estimator()
{
    if (bead_num != 0) return 0.0;

    prim_est[0] = 0.0;

    for (int m = 1; m < nbosons + 1; ++m) {
        double sig = 0.0;

        // Numerical stability (Xiong & Xiong method)
        double Elongest = std::numeric_limits<double>::max();

        for (int k = m; k > 0; k--) {
            Elongest = std::min(Elongest, get_Enk(m, k) + V[m - k]);
        }

        for (int k = m; k > 0; --k) {
            double E_kn_val = get_Enk(m, k);

            sig += (prim_est[m - k] - E_kn_val) * exp(-beta * (E_kn_val + V[m - k] - Elongest));
        }

        double sig_denom_m = m * exp(-beta * (V[m] - Elongest));

        prim_est[m] = sig / sig_denom_m;
    }

#if IPI_CONVENTION
    return prim_est[nbosons] / np;
#else
    return prim_est[nbosons];
#endif
}

/* ---------------------------------------------------------------------- */