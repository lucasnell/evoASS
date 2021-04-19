

#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <random>
#include <vector>
#include <deque>

#include "sim.hpp"
#include "pcg.hpp"
#include "quant_gen.hpp"

#ifdef _OPENMP
#include <omp.h>  // omp
#endif


using namespace Rcpp;



/*
 Note that this is the "instantaneous" version of quantitative genetics.
 It's not the same as what's written in Tobin and Tony's methods section.
 Instead of adding new daughter species, it starts with a number of species
 that should be more than you would ever expect to be produced sequentially.
 */




//' Derivative of fitness with respect to the trait divided by mean fitness.
//'
//' The function below calculates this selection strength for all traits
//' for all species.
//'
//' @noRd
//'
inline void sel_str__(arma::mat& ss_mat,
                      const std::vector<arma::vec>& V,
                      const std::vector<double>& N,
                      const double& f,
                      const double& a0,
                      const arma::mat& C,
                      const double& r0,
                      const arma::mat& D) {

    uint32_t n = V.size();      // # species
    uint32_t q = V[0].n_elem;   // # traits

    if (ss_mat.n_rows != q || ss_mat.n_cols != n) {
        ss_mat.set_size(q, n);
    }

    /*
     For all `i`, calculate `sum(N_j * exp(- transpose(V_j) * D * V_j))`,
     for all `j != i`.
     */
    arma::vec W(n, arma::fill::zeros);
    for (uint32_t j = 0; j < n; j++) {
        // Calculate `N_j * exp(-transpose(V_j) * D * V_j)`:
        double W_ = N[j] * std::exp(-1 * arma::as_scalar(V[j].t() * D * V[j]));
        // Now insert this value at all `i` where `j != i`:
        for (uint32_t i = 0; i < n; i++) {
            if (i == j) continue;
            W(i) += W_;
        }
    }

    // Now go back through and calculate strength of selection:
    arma::rowvec Vtmp(q);
    for (uint32_t i = 0; i < n; i++) {
        const arma::vec& Vi(V[i]);
        Vtmp = 2 * (-f * Vi.t() * C +
            a0 * W[i] * std::exp(arma::as_scalar(-Vi.t() * Vi)) * Vi.t());
        ss_mat.col(i) = Vtmp.t();
    }

    return;

}



//' R-exported version of above, so it can be tested in R for accuracy.
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat sel_str_cpp(const arma::mat& V,
                      const std::vector<double>& N,
                      const double& f,
                      const double& a0,
                      const arma::mat& C,
                      const double& r0,
                      const arma::mat& D) {

    uint32_t n = V.n_cols;

    if (n != N.size()) stop("V.n_cols != N.size()");

    std::vector<arma::vec> VV;
    VV.reserve(n);
    for (uint32_t i = 0; i < n; i++) VV.push_back(V.col(i));

    arma::mat ss_mat;

    sel_str__(ss_mat, VV, N, f, a0, C, r0, D);

    return ss_mat;
}



//' Partial derivative of species i traits at time t+1 with respect to
//' species i traits at time t.
//'
//'
//' @noRd
//'
inline void dVi_dVi_(arma::mat& dVhat,
                     const uint32_t& row_start,
                     const uint32_t& col_start,
                     const arma::vec& Vi,
                     const double& Omega,
                     const arma::mat& C,
                     const double& f,
                     const double& a0,
                     const double& add_var) {
    uint32_t q = Vi.n_elem;
    arma::mat I = arma::eye<arma::mat>(q, q);
    uint32_t row_end = row_start + q - 1;
    uint32_t col_end = col_start + q - 1;
    dVhat(arma::span(row_start, row_end), arma::span(col_start, col_end)) = I +
        2 * add_var * (
                (
                        a0 * Omega * std::exp(-1 *
                            arma::as_scalar(Vi.t() * Vi)) *
                            (I - 2 * Vi * Vi.t())
                ) - (f * C.t())
        );

    return;
}


//' R-exported version of above, to be used in R for testing.
//'
//' NOTE: This does NOT account for step function to keep traits >= 0
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat dVi_dVi_cpp(const uint32_t& i, const arma::mat& V,
                      const double& Omega,
                      const arma::mat& C, const double& f, const double& a0,
                      const double& add_var) {

    if (!C.is_symmetric()) stop("C must be symmetric");

    arma::mat dVhat(V.n_rows, V.n_rows);

    // Fill dVhat:
    dVi_dVi_(dVhat, 0, 0, V.col(i), Omega, C, f, a0, add_var);

    return dVhat;
}

//' Partial derivative of species i traits at time t+1 with respect to
//' species k traits at time t.
//'
//' @noRd
//'
//'
inline void dVi_dVk_(arma::mat& dVhat,
                     const uint32_t& row_start,
                     const uint32_t& col_start,
                     const double& Nk,
                     const arma::vec& Vi,
                     const arma::vec& Vk,
                     const arma::mat& D,
                     const double& a0,
                     const double& add_var) {
    uint32_t row_end = row_start + Vi.n_elem - 1;
    uint32_t col_end = col_start + Vi.n_elem - 1;
    arma::mat M = Vi * arma::exp(-1 * Vk.t() * D * Vk - Vi.t() * Vi) *
        Vk.t() * D;
    M *= (-4 * a0 * add_var * Nk);
    dVhat(arma::span(row_start, row_end), arma::span(col_start, col_end)) = M;
    return;
}


//' R-exported version of above, to be used in R for testing.
//'
//' NOTE: This does NOT account for step function to keep traits >= 0
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat dVi_dVk_cpp(const uint32_t& i,
                      const uint32_t& k,
                      const std::vector<double>& N,
                      const arma::mat& V,
                      const arma::mat& D,
                      const double& a0,
                      const double& add_var) {

    if (!D.is_symmetric()) stop("D must be symmetric");

    arma::mat dVhat(V.n_rows, V.n_rows);
    // Fill dVhat:
    dVi_dVk_(dVhat, 0, 0, N[k], V.col(i), V.col(k), D, a0, add_var);

    return dVhat;
}




//' Partial derivative of species i traits at time t+1 with respect to
//' species i abundance at time t.
//'
//' @noRd
//'
//'
inline void dVi_dNi_(arma::mat& dVhat,
                     const uint32_t& row_start,
                     const uint32_t& col_start,
                     const arma::vec& Vi,
                     const double& a0,
                     const double& add_var) {

    uint32_t row_end = row_start + Vi.n_elem - 1;
    const uint32_t& col_end(col_start);

    dVhat(arma::span(row_start, row_end), arma::span(col_start, col_end)).fill(0);

    return;
}


//' R-exported version of above, to be used in R for testing.
//'
//' NOTE: This does NOT account for step function to keep traits >= 0
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat dVi_dNi_cpp(const uint32_t& i,
                      const arma::mat& V,
                      const double& a0,
                      const double& add_var) {

    arma::mat dVhat(V.n_rows, 1);
    // Fill dVhat:
    dVi_dNi_(dVhat, 0, 0, V.col(i), a0, add_var);

    return dVhat;
}




//' Partial derivative of species i traits at time t+1 with respect to
//' species k abundance at time t.
//'
//' @noRd
//'
//'
inline void dVi_dNk_(arma::mat& dVhat,
                     const uint32_t& row_start,
                     const uint32_t& col_start,
                     const arma::vec& Vi,
                     const arma::vec& Vk,
                     const arma::mat& D,
                     const double& a0,
                     const double& add_var) {

    uint32_t row_end = row_start + Vi.n_elem - 1;
    const uint32_t& col_end(col_start);

    arma::mat M = 2 * add_var * a0 * Vi *
        arma::exp(- Vk.t() * D * Vk - Vi.t() * Vi);

    dVhat(arma::span(row_start, row_end), arma::span(col_start, col_end)) = M;

    return;
}


//' R-exported version of above, to be used in R for testing.
//'
//' NOTE: This does NOT account for step function to keep traits >= 0
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat dVi_dNk_cpp(const uint32_t& i,
                      const uint32_t& k,
                      const arma::mat& V,
                      const arma::mat& D,
                      const double& a0,
                      const double& add_var) {

    if (!D.is_symmetric()) stop("D must be symmetric");

    arma::mat dVhat(V.n_rows, 1);
    // Fill dVhat:
    dVi_dNk_(dVhat, 0, 0, V.col(i), V.col(k), D, a0, add_var);

    return dVhat;
}




//' Partial derivative of species i abundance at time t+1 with respect to
//' species i traits at time t.
//'
//' @noRd
//'
//'
inline void dNi_dVi_(arma::mat& dVhat,
                     const uint32_t& row_start,
                     const uint32_t& col_start,
                     const uint32_t& i,
                     const arma::mat& V,
                     const std::vector<double>& N,
                     const double& f,
                     const double& a0,
                     const arma::mat& C,
                     const double& r0,
                     const arma::mat& D) {

    const uint32_t& row_end(row_start);
    uint32_t col_end = col_start + V.n_rows - 1;

    double F = F_it__(i, V, N, f, a0, C, r0, D);

    double Omega = 0;
    for (uint32_t j = 0; j < N.size(); j++) {
        if (j != i) {
            Omega += (N[j] *
                std::exp(-1 * arma::as_scalar(V.col(j).t() * D * V.col(j))));
        }
    }

    const arma::subview_col<double>& Vi(V.col(i));

    arma::mat M = 2 * F * N[i] * (
        a0 * Omega * arma::exp(-1 * Vi.t() * Vi) * Vi.t() -
        f * Vi.t() * C);

    dVhat(arma::span(row_start, row_end), arma::span(col_start, col_end)) = M;

    return;
}


//' R-exported version of above, to be used in R for testing.
//'
//' NOTE: This does NOT account for step function to keep traits >= 0
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat dNi_dVi_cpp(const uint32_t& i,
                      const arma::mat& V,
                      const std::vector<double>& N,
                      const double& f,
                      const double& a0,
                      const arma::mat& C,
                      const double& r0,
                      const arma::mat& D) {

    if (!C.is_symmetric()) stop("C must be symmetric");
    if (!D.is_symmetric()) stop("D must be symmetric");

    arma::mat dVhat(1, V.n_rows);
    // Fill dVhat:
    dNi_dVi_(dVhat, 0, 0, i, V, N, f, a0, C, r0, D);

    return dVhat;
}




//' Partial derivative of species i abundance at time t+1 with respect to
//' species k traits at time t.
//'
//' @noRd
//'
//'
inline void dNi_dVk_(arma::mat& dVhat,
                     const uint32_t& row_start,
                     const uint32_t& col_start,
                     const uint32_t& i,
                     const uint32_t& k,
                     const arma::mat& V,
                     const std::vector<double>& N,
                     const double& f,
                     const double& a0,
                     const arma::mat& C,
                     const double& r0,
                     const arma::mat& D) {

    const uint32_t& row_end(row_start);
    uint32_t col_end = col_start + V.n_rows - 1;

    double F = F_it__(i, V, N, f, a0, C, r0, D);

    double Omega = 0;
    for (uint32_t j = 0; j < N.size(); j++) {
        if (j != i) {
            Omega += (N[j] *
                std::exp(-1 * arma::as_scalar(V.col(j).t() * D * V.col(j))));
        }
    }

    const arma::subview_col<double>& Vi(V.col(i));
    const arma::subview_col<double>& Vk(V.col(k));

    arma::mat M = 2 * F * N[i] * N[k] * a0 *
        arma::exp(- Vk.t() * D * Vk - Vi.t() * Vi) *
        Vk.t() * D;


    dVhat(arma::span(row_start, row_end), arma::span(col_start, col_end)) = M;

    return;
}


//' R-exported version of above, to be used in R for testing.
//'
//' NOTE: This does NOT account for step function to keep traits >= 0
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat dNi_dVk_cpp(const uint32_t& i,
                      const uint32_t& k,
                      const arma::mat& V,
                      const std::vector<double>& N,
                      const double& f,
                      const double& a0,
                      const arma::mat& C,
                      const double& r0,
                      const arma::mat& D) {

    if (!C.is_symmetric()) stop("C must be symmetric");
    if (!D.is_symmetric()) stop("D must be symmetric");


    arma::mat dVhat(1, V.n_rows);
    // Fill dVhat:
    dNi_dVk_(dVhat, 0, 0, i, k, V, N, f, a0, C, r0, D);

    return dVhat;
}







//' Partial derivative of species i abundance at time t+1 with respect to
//' species i abundance at time t.
//'
//' @noRd
//'
//'
inline void dNi_dNi_(arma::mat& dVhat,
                     const uint32_t& row_start,
                     const uint32_t& col_start,
                     const uint32_t& i,
                     const arma::mat& V,
                     const std::vector<double>& N,
                     const double& f,
                     const double& a0,
                     const arma::mat& C,
                     const double& r0,
                     const arma::mat& D) {

    double F = F_it__(i, V, N, f, a0, C, r0, D);

    double M = F * (1 - N[i] * a0);

    dVhat(row_start, col_start) = M;

    return;
}


//' R-exported version of above, to be used in R for testing.
//'
//' NOTE: This does NOT account for step function to keep traits >= 0
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat dNi_dNi_cpp(const uint32_t& i,
                      const arma::mat& V,
                      const std::vector<double>& N,
                      const double& f,
                      const double& a0,
                      const arma::mat& C,
                      const double& r0,
                      const arma::mat& D) {

    if (!C.is_symmetric()) stop("C must be symmetric");
    if (!D.is_symmetric()) stop("D must be symmetric");

    arma::mat dVhat(1, 1);
    // Fill dVhat:
    dNi_dNi_(dVhat, 0, 0, i, V, N, f, a0, C, r0, D);

    return dVhat;
}





//' Partial derivative of species i abundance at time t+1 with respect to
//' species k abundance at time t.
//'
//' @noRd
//'
//'
inline void dNi_dNk_(arma::mat& dVhat,
                     const uint32_t& row_start,
                     const uint32_t& col_start,
                     const uint32_t& i,
                     const uint32_t& k,
                     const arma::mat& V,
                     const std::vector<double>& N,
                     const double& f,
                     const double& a0,
                     const arma::mat& C,
                     const double& r0,
                     const arma::mat& D) {

    double F = F_it__(i, V, N, f, a0, C, r0, D);

    const arma::subview_col<double>& Vi(V.col(i));
    const arma::subview_col<double>& Vk(V.col(k));

    double M = -1 * F * N[i] * a0 *
        std::exp(arma::as_scalar(-1 * Vi.t() * Vi - Vk.t() * D * Vk));

    dVhat(row_start, col_start) = M;

    return;
}


//' R-exported version of above, to be used in R for testing.
//'
//' NOTE: This does NOT account for step function to keep traits >= 0
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat dNi_dNk_cpp(const uint32_t& i,
                      const uint32_t& k,
                      const arma::mat& V,
                      const std::vector<double>& N,
                      const double& f,
                      const double& a0,
                      const arma::mat& C,
                      const double& r0,
                      const arma::mat& D) {

    if (!C.is_symmetric()) stop("C must be symmetric");
    if (!D.is_symmetric()) stop("D must be symmetric");

    arma::mat dVhat(1, 1);
    // Fill dVhat:
    dNi_dNk_(dVhat, 0, 0, i, k, V, N, f, a0, C, r0, D);

    return dVhat;
}





//
// Account for step function to keep traits >= 0
//
inline void correct_jac(arma::mat& jac,
                        const arma::mat& V,
                        const std::vector<double>& N,
                        const double& f,
                        const double& a0,
                        const double& r0,
                        const arma::mat& D,
                        const arma::mat& C,
                        const arma::vec& add_var,
                        const bool& evo_only) {

    std::vector<arma::vec> VV;
    VV.reserve(V.n_cols);
    for (uint32_t i = 0; i < V.n_cols; i++) VV.push_back(V.col(i));

    arma::mat S = arma::diagmat(add_var);

    arma::mat ss;
    sel_str__(ss, VV, N, f, a0, C, r0, D);
    arma::mat deltaV = ss * S;

    arma::vec newV = arma::vectorise(V + deltaV);
    for (double& d : newV) d = (d > 0) ? 1 : 0;  // <-- Heaviside function

    if (evo_only) {
        jac = arma::diagmat(newV) * jac;
    } else {
        arma::vec newVN(newV.n_elem + V.n_cols);
        newVN.head(newV.n_elem) = newV;
        newVN.tail(V.n_cols).ones(); // No Heaviside for N derivatives!
        jac = arma::diagmat(newVN) * jac;
    }

    return;
}


//' Calculate the Jacobian of first derivatives.
//'
//' Cell [i,j] contains the partial derivative of j with respect to i.
//'
//' NOTE: This DOES account for step function to keep traits >= 0
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat jacobian_cpp(const arma::mat& V,
                       const std::vector<double>& N,
                       const double& f,
                       const double& a0,
                       const double& r0,
                       const arma::mat& D,
                       const arma::mat& C,
                       const arma::vec& add_var,
                       const bool& evo_only) {

    if (!C.is_symmetric()) stop("C must be symmetric");
    if (!D.is_symmetric()) stop("D must be symmetric");

    uint32_t n = N.size();
    uint32_t q = V.n_rows;

    if (V.n_cols != n) stop("V.n_cols != N.size()");
    if (add_var.n_elem != n) stop("add_var.n_elem != N.size()");

    arma::mat jcb_mat;
    if (evo_only) {
        jcb_mat.set_size(n*q, n*q);
    } else jcb_mat.set_size(n*(q+1), n*(q+1));


    /*
     ---------
     1. ∂ Traits / ∂ Traits
     ---------
     */

    arma::vec Omega_vec(n);
    for (uint32_t j = 0; j < n; j++) {
        Omega_vec[j] = (N[j] * std::exp(
            arma::as_scalar(-1 * V.col(j).t() * D * V.col(j))));
    }

    for (uint32_t i = 0; i < n; i++) {

        uint32_t row_start = i * q;

        for (uint32_t k = 0; k < n; k++) {

            uint32_t col_start = k * q;

            if (k == i) {

                double Omega = 0;
                for (uint32_t j = 0; j < n; j++) {
                    if (j != i) Omega += Omega_vec[j];
                }
                // Fill Jacobian:
                dVi_dVi_(jcb_mat, row_start, col_start, V.col(i), Omega,
                         C, f, a0, add_var[i]);

            } else {

                // Fill Jacobian:
                dVi_dVk_(jcb_mat, row_start, col_start, N[k], V.col(i),
                         V.col(k), D, a0, add_var[i]);

            }

        }

    }

    if (evo_only) {
        // account for step function to keep traits >= 0
        correct_jac(jcb_mat, V, N, f, a0, r0, D, C, add_var, evo_only);
        return(jcb_mat);
    }


    /*
     ---------
     2. ∂ Traits / ∂ Abundances
     ---------
     */

    for (uint32_t i = 0; i < n; i++) {

        uint32_t row_start = i * q;

        for (uint32_t k = 0; k < n; k++) {

            uint32_t col_start = n * q + k;

            if (k == i) {

                // Fill Jacobian:
                dVi_dNi_(jcb_mat, row_start, col_start, V.col(i),
                         a0, add_var[i]);

            } else {

                // Fill Jacobian:
                dVi_dNk_(jcb_mat, row_start, col_start, V.col(i), V.col(k),
                         D, a0, add_var[i]);

            }

        }

    }



    /*
     ---------
     3. ∂ Abundances / ∂ Traits
     ---------
     */

    for (uint32_t i = 0; i < n; i++) {

        uint32_t row_start = n * q + i;

        for (uint32_t k = 0; k < n; k++) {

            uint32_t col_start = k * q;

            if (k == i) {

                // Fill Jacobian:
                dNi_dVi_(jcb_mat, row_start, col_start,
                         i, V, N, f, a0, C, r0, D);

            } else {

                // Fill Jacobian:
                dNi_dVk_(jcb_mat, row_start, col_start,
                         i, k, V, N, f, a0, C, r0, D);

            }

        }

    }

    /*
     ---------
     4. ∂ Abundances / ∂ Abundances
     ---------
     */

    for (uint32_t i = 0; i < n; i++) {

        uint32_t row_start = n * q + i;

        for (uint32_t k = 0; k < n; k++) {

            uint32_t col_start = n * q + k;

            if (k == i) {

                // Fill Jacobian:
                dNi_dNi_(jcb_mat, row_start, col_start,
                         i, V, N, f, a0, C, r0, D);

            } else {

                // Fill Jacobian:
                dNi_dNk_(jcb_mat, row_start, col_start,
                         i, k, V, N, f, a0, C, r0, D);

            }

        }

    }

    // account for step function to keep traits >= 0
    correct_jac(jcb_mat, V, N, f, a0, r0, D, C, add_var, evo_only);


    return jcb_mat;
}







//' Search for unique species in a matrix of species trait values.
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::uvec unq_spp_cpp(const std::vector<arma::vec>& V,
                       double precision) {

    precision *= precision;

    uint32_t n = V.size();
    arma::uvec is_unq = arma::ones<arma::uvec>(n);

    for (uint32_t i = 1; i < n; i++) {
        for (uint32_t j = 0; j < i; j++) {

            // don't want to keep looking at non-unique spp:
            if (is_unq(j) == 0) continue;

            arma::vec diff_;
            diff_ = (V[i] - V[j]) % (V[i] - V[j]);
            if (arma::mean(diff_) < precision) {
                is_unq(i) = 0;
                break;
            }
        }
    }

    return is_unq;
}


/*
 Similar to above, except that it groups species in to unique groups.
 */
//[[Rcpp::export]]
IntegerVector group_spp_cpp(const std::vector<arma::vec>& V,
                                    double precision) {


    precision *= precision;

    uint32_t n = V.size();

    IntegerVector spp_groups(n, 0);

    arma::uvec is_unq = arma::ones<arma::uvec>(n);

    uint32_t g = 1;

    for (uint32_t i = 1; i < n; i++) {
        for (uint32_t j = 0; j < i; j++) {

            // don't want to keep looking at non-unique spp:
            if (is_unq(j) == 0) continue;

            arma::vec diff_;
            diff_ = (V[i] - V[j]) % (V[i] - V[j]);
            if (arma::mean(diff_) < precision) {
                is_unq(i) = 0;
                spp_groups[i] = spp_groups[j];
                break;
            }
        }
        if (is_unq(i) == 1) {
            spp_groups[i] = g;
            g++;
        }
    }

    return spp_groups;
}




//' One repetition of quantitative genetics.
//'
//' Higher-up function(s) should handle the info put into `info`.
//'
//'
//' @noRd
//'
int one_quant_gen__(OneRepInfo& info,
                    std::deque<arma::vec> V0,
                    std::deque<arma::vec> Vp0,
                    std::deque<double> N0,
                    const double& f,
                    const double& a0,
                    const arma::mat& C,
                    const double& r0,
                    const arma::mat& D,
                    std::deque<double> add_var,
                    const double& sigma_V0,
                    const double& sigma_N,
                    const std::vector<double>& sigma_V,
                    const uint32_t& spp_gap_t,
                    const uint32_t& final_t,
                    const double& min_N,
                    const bool& adjust_mu_V,
                    const bool& lnorm_V,
                    const uint32_t& save_every,
                    pcg64& eng,
                    Progress& prog_bar) {


    uint32_t n = N0.size();
    uint32_t q = V0.front().n_elem;

    /*
     If using a lognormal distribution for evolution stochasticity
     and if `adjust_mu_V == TRUE`,
     this keeps the mean of the transformed distribution equal to 1.
     `mu_V` isn't used when a lognormal distribution isn't used, so we'll
     adjust it even when a `lnorm_V == TRUE`.
     */
    std::vector<double> mu_V(q, 0);
    if (adjust_mu_V) {
        for (uint32_t j = 0; j < q; j++) {
            mu_V[j] -= (sigma_V[j] * sigma_V[j]) / 2;
        }
    }

    normal_distr distr = normal_distr(0, 1);

    // adding stochasticity to starting genotypes (and phenotypes if desired)
    if (sigma_V0 > 0) {
        Vp0 = V0; // mostly just to resize `Vp0`
        for (uint32_t i = 0; i < n; i++) {
            for (uint32_t j = 0; j < q; j++) {
                V0[i][j] = trunc_rnorm_(V0[i][j], sigma_V0, eng);
                Vp0[i][j] = V0[i][j];
                if (sigma_V[j] > 0) {
                    if (lnorm_V) {
                        Vp0[i][j] *= std::exp(distr(eng) * sigma_V[j] + mu_V[j]);
                    } else {
                        Vp0[i][j] = trunc_rnorm_(Vp0[i][j], sigma_V[j], eng);
                    }
                }
            }
        }
    }

    /*
     adjusting starting phenotypes if `sigma_V0 == 0`
     */
    if (Vp0.size() == 0) {
        Vp0 = V0;
        for (uint32_t j = 0; j < q; j++) {
            if (sigma_V[j] > 0) {
                if (lnorm_V) {
                    for (uint32_t i = 0; i < n; i++) {
                        Vp0[i][j] *= std::exp(distr(eng) * sigma_V[j] + mu_V[j]);
                    }
                } else {
                    for (uint32_t i = 0; i < n; i++) {
                        Vp0[i][j] = trunc_rnorm_(Vp0[i][j], sigma_V[j], eng);
                    }
                }
            }
        }
    }


    if (spp_gap_t == 0) {
        info = OneRepInfo(N0, V0, Vp0, add_var);
        N0.clear();
        V0.clear();
        Vp0.clear();
        add_var.clear();
    } else {
        info = OneRepInfo(N0.front(), V0.front(), Vp0.front(), add_var.front());
        N0.pop_front();
        V0.pop_front();
        Vp0.pop_front();
        add_var.pop_front();
    }


    // Setting size for `info` fields
    if (save_every > 0) {
        uint32_t spp_add_saves = static_cast<uint32_t>(std::ceil(
            static_cast<double>(spp_gap_t) / static_cast<double>(save_every)));
        spp_add_saves += 2U;
        uint32_t final_saves = static_cast<uint32_t>(std::ceil(
            static_cast<double>(final_t) / static_cast<double>(save_every)));
        final_saves += 2U;
        info.reserve(final_saves + (info.n + V0.size() - 1) * spp_add_saves);
    }

    uint32_t t = 0;
    bool all_gone = false;
    uint32_t interrupt_iters = 0;   // checking for user interrupt
    uint32_t n_pb_incr = 0;         // progress bar increments


    // Save starting info:
    if (save_every > 0) info.save_time(t);


    // First iterations with species additions
    bool new_spp = false;
    while (!N0.empty()) {

        n_pb_incr++;

        // Update abundances and traits:
        all_gone = info.iterate(f, a0, C, r0, D, min_N,
                                sigma_N, sigma_V, mu_V, lnorm_V, eng);

        // Add new species if necessary:
        new_spp = (t + 1) == (info.n * spp_gap_t);
        if (new_spp) {
            info.add_species(N0.front(), V0.front(), Vp0.front(),
                             add_var.front());
            N0.pop_front();
            V0.pop_front();
            Vp0.pop_front();
            add_var.pop_front();
        }

        if (save_every > 0 && (t % save_every == 0 || new_spp)) {
            info.save_time(t + 1);
        }

        if (n_pb_incr > 100) {
            prog_bar.increment(n_pb_incr);
            n_pb_incr = 0;
        }

        t++;


        // Check for user interrupt:
        if (interrupt_check(interrupt_iters, prog_bar, 100)) return -1;

    }

    if (final_t == 0) return 0;

    uint32_t total_time = final_t + t;

    // Final iterations with no species additions
    while (!all_gone && t < total_time) {

        n_pb_incr++;

        // Update abundances and traits:
        all_gone = info.iterate(f, a0, C, r0, D, min_N,
                                sigma_N, sigma_V, mu_V, lnorm_V, eng);

        if (save_every > 0 &&
            (t % save_every == 0 || (t+1) == final_t || all_gone)) {
            info.save_time(t + 1);
        }

        if (n_pb_incr > 100) {
            prog_bar.increment(n_pb_incr);
            n_pb_incr = 0;
        }

        t++;

        // Check for user interrupt:
        if (interrupt_check(interrupt_iters, prog_bar, 100)) return -1;

    }

    if (n_pb_incr > 0) prog_bar.increment(n_pb_incr);

    return 0;
}


//' Multiple repetitions of quantitative genetics.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat quant_gen_cpp(const uint32_t& n_reps,
                        const std::deque<arma::vec>& V0,
                        const std::deque<arma::vec>& Vp0,
                        const std::deque<double>& N0,
                        const double& f,
                        const double& a0,
                        const arma::mat& C,
                        const double& r0,
                        const arma::mat& D,
                        const std::deque<double>& add_var,
                        const double& sigma_V0,
                        const double& sigma_N,
                        const std::vector<double>& sigma_V,
                        const uint32_t& spp_gap_t,
                        const uint32_t& final_t,
                        const double& min_N,
                        const bool& adjust_mu_V,
                        const bool& lnorm_V,
                        const uint32_t& save_every,
                        const bool& show_progress,
                        const uint32_t& n_threads) {

    if (!C.is_symmetric()) stop("C must be symmetric");
    if (!D.is_symmetric()) stop("D must be symmetric");

    /*
     Adding stochasticity to starting genotypes invalidates starting
     phenotypes, so I'm not allowing both.
     */
    if (sigma_V0 > 0 && Vp0.size() > 0) {
        stop("\nproviding Vp0 with sigma_V0 > 0 makes no sense");
    }

    const uint32_t n = N0.size();

    if (n == 0) stop("n == 0");

    if (V0.size() != n) stop("V0.size() != n");
    if (add_var.size() != n) stop("add_var.size() != n");
    if (Vp0.size() > 0 && Vp0.size() != n) stop("Vp0.size() != n");


    const uint32_t q = V0[0].n_elem;
    if (C.n_cols != q) stop("C.n_cols != q");
    if (C.n_rows != q) stop("C.n_rows != q");
    if (D.n_cols != q) stop("D.n_cols != q");
    if (D.n_rows != q) stop("D.n_rows != q");

    std::vector<OneRepInfo> rep_infos(n_reps);

    const std::vector<std::vector<uint128_t>> seeds = mc_seeds_rep(n_reps);

    Progress prog_bar(n_reps * (final_t + (n - 1) * spp_gap_t), show_progress);
    bool interrupted = false;

    #ifdef _OPENMP
    #pragma omp parallel default(shared) num_threads(n_threads) if (n_threads > 1)
    {
    #endif


    #ifdef _OPENMP
    uint32_t active_thread = omp_get_thread_num();
    #else
    uint32_t active_thread = 0;
    #endif

    pcg64 eng;

    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (uint32_t i = 0; i < n_reps; i++) {
        eng.seed(seeds[i][0], seeds[i][1]);
        int status = one_quant_gen__(rep_infos[i], V0, Vp0, N0, f, a0, C, r0, D,
                                     add_var, sigma_V0, sigma_N, sigma_V,
                                     spp_gap_t, final_t, min_N, adjust_mu_V,
                                     lnorm_V,
                                     save_every, eng, prog_bar);

        if (active_thread == 0 && status != 0) interrupted = true;
    }
    #ifdef _OPENMP
    }
    #endif

    if (interrupted) {
        throw(Rcpp::exception("\nUser interrupted process.", false));
    }

    /*
     ------------
     Now organize output:
     ------------
     */
    arma::mat nv; // N and V - either through time or just final values

    /*
     Go through one time to calculate the # surviving species for all reps and
     for all time point(s) saved.
     Then go back through and fill in values.
     */
    uint32_t total_n_spp = 0;
    if (save_every > 0) {  //   ---- Saving values through time: ----

        for (uint32_t i = 0; i < n_reps; i++) {
            for (uint32_t j = 0; j < rep_infos[i].t.size(); j++) {
                total_n_spp += rep_infos[i].N_t[j].size();
            }
        }
        nv.set_size(total_n_spp, 4 + 2 * q);

        uint32_t j = 0;
        for (uint32_t i = 0; i < n_reps; i++) {

            Rcpp::checkUserInterrupt();
            const OneRepInfo& info(rep_infos[i]);
            for (uint32_t t = 0; t < info.t.size(); t++) {

                const std::vector<double>& N_t(info.N_t[t]);
                const std::vector<arma::vec>& V_t(info.V_t[t]);
                const std::vector<arma::vec>& Vp_t(info.Vp_t[t]);
                const std::vector<uint32_t>& spp_t(info.spp_t[t]);
                const double& t_(info.t[t]);
                for (uint32_t k = 0; k < N_t.size(); k++) {

                    nv(j+k,0) = i + 1;      // rep
                    nv(j+k,1) = t_;         // time
                    nv(j+k,2) = spp_t[k];   // species
                    nv(j+k,3) = N_t[k];     // N
                    // V and Vp:
                    for (uint32_t l = 0; l < q; l++) {
                        nv(j+k, 4+l) = V_t[k](l);       // V
                        nv(j+k, 4+q+l) = Vp_t[k](l);    //  Vp
                    }

                }

                j += N_t.size();

            }

        }

    } else {  //                ----  Just saving final values:  ----

        for (uint32_t i = 0; i < n_reps; i++) {
            total_n_spp += rep_infos[i].N.size();
            if (rep_infos[i].N.empty()) total_n_spp++;
        }
        nv.set_size(total_n_spp, 3 + 2 * q);
        uint32_t j = 0;
        for (uint32_t i = 0; i < n_reps; i++) {
            Rcpp::checkUserInterrupt();
            const OneRepInfo& info(rep_infos[i]);
            if (!info.N.empty()) {
                for (uint32_t k = 0; k < info.N.size(); k++) {
                    nv(j+k,0) = i + 1;          // rep
                    nv(j+k,1) = info.spp[k];    // species
                    nv(j+k,2) = info.N[k];      // N
                    // V and Vp:
                    for (uint32_t l = 0; l < q; l++) {
                        nv(j+k, 3+l) = info.V[k](l);
                        nv(j+k, 3+q+l) = info.Vp[k](l);
                    }
                }
                j += info.N.size();
            } else {
                nv(j,0) = i + 1;      // rep
                nv(j,1) = 0;          // species
                nv(j,2) = 0;          // N
                // V and Vp:
                for (uint32_t l = 0; l < q; l++) {
                    nv(j, 3+l) = arma::datum::nan;
                    nv(j, 3+q+l) = arma::datum::nan;
                }
                j++;
            }
        }

    }

    return nv;

}

