
// #define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <random>

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
void sel_str__(arma::mat& ss_mat,
               const std::vector<arma::rowvec>& V,
               const std::vector<double>& N,
               const double& f,
               const double& g,
               const arma::mat& C,
               const double& r0,
               const double& d) {

    uint32_t n_spp = V.size();      // # species
    uint32_t n_trt = V[0].size();   // # traits

    if (ss_mat.n_rows != n_spp || ss_mat.n_cols != n_trt) ss_mat.set_size(n_spp, n_trt);

    // For all i, calculate `N_j * exp(-d * V_j * transpose(V_j))`, where `j != i`
    arma::vec W(n_spp, arma::fill::zeros);
    for (uint32_t j = 0; j < n_spp; j++) {
        // Calculate `N_j * exp(-d * V_j * transpose(V_j))`:
        double W_ = N[j] * std::exp(-d * arma::as_scalar(V[j] * V[j].t()));
        // Now insert this value at all `i` where `j != i`:
        for (uint32_t i = 0; i < n_spp; i++) {
            if (i == j) continue;
            W(i) += W_;
        }
    }
    arma::mat C_ = C + C.t();

    // Now go back through and calculate strength of selection:
    for (uint32_t i = 0; i < n_spp; i++) {
        const arma::rowvec& V_i(V[i]);
        const double& N_i(N[i]);
        ss_mat.row(i) = (-f * V_i * C_) +
            2 * g * V_i * std::exp(arma::as_scalar(-V_i * V_i.t())) * (N_i + W(i));
    }

    return;

}


//' R-exported version of above, so it can be tested in R for accuracy.
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat sel_str_cpp(const std::vector<arma::rowvec>& V,
                      const std::vector<double>& N,
                      const double& f,
                      const double& g,
                      const arma::mat& C,
                      const double& r0,
                      const double& d) {

    arma::mat ss_mat;

    sel_str__(ss_mat, V, N, f, g, C, r0, d);

    return ss_mat;
}


//'
//' Calculates for one species' traits at a time.
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::rowvec dF_dVi_cpp(const arma::rowvec& V_i,
                        const std::vector<arma::rowvec>& V_nei,
                        const double& N_i,
                        const std::vector<double>& N_nei,
                        const double& f,
                        const double& g,
                        const arma::mat& C,
                        const double& r0,
                        const double& d) {

    if (N_nei.size() != V_nei.size()) stop("N_nei.size() != V_nei.size()");
    if (V_i.n_elem != V_nei[0].n_elem) stop("V_i.n_elem != V_nei[0].n_elem");
    if (C.n_cols != V_i.n_elem || C.n_rows != V_i.n_elem) {
        stop("C.n_cols != V_i.n_elem || C.n_rows != V_i.n_elem");
    }

    double F = F_t_deriv_cpp(V_i, V_nei, N_i, N_nei, f, g, C, r0, d);

    double tmp_dbl = 0;
    for (uint32_t j = 0; j < V_nei.size(); j++) {
        const arma::rowvec& V_j(V_nei[j]);
        const double& N_j(N_nei[j]);
        tmp_dbl += (N_j / std::exp(d * arma::as_scalar(V_j * V_j.t())));
    }

    arma::rowvec dF_dVi_vec = (-f * V_i * (C + C.t())) +
        2 * V_i * g * std::exp(arma::as_scalar(-V_i * V_i.t())) * (N_i + tmp_dbl);

    dF_dVi_vec *= F;

    return dF_dVi_vec;

}





//' One round of quantitative genetics.
//'
//' Higher-up function(s) should create `C` from `eta` and should handle any outputs.
//'
//'
//' @noRd
//'
void one_quantgen_rep(OneRepInfo& info,
                      const std::vector<arma::rowvec>& V0,
                      const std::vector<double>& N0,
                      const double& f,
                      const double& g,
                      const arma::mat& C,
                      const double& r0,
                      const double& d,
                      const double& add_var,
                      const double& delta,
                      const uint32_t& start_t,
                      const uint32_t& max_t,
                      const double& min_N,
                      const uint32_t& save_every,
                      pcg64& eng) {


    if (save_every > 0) {
        info = OneRepInfo(N0, V0, max_t, save_every);
    } else {
        info = OneRepInfo(N0, V0);
    }

    bool all_gone = false;
    std::lognormal_distribution<double> distr(0.0, delta);

    for (uint32_t t = 0; t < start_t; t++) {

        // Update abundances and traits:
        all_gone = info.iterate(f, g, C, r0, d, min_N, true);
        if (all_gone) return;

    }

    // perturb trait values
    info.perturb(eng, distr);

    for (uint32_t t = 0; t < max_t; t++) {

        // Update abundances and traits:
        all_gone = info.iterate(f, g, C, r0, d, min_N, true);
        if (all_gone) return;

        if (save_every > 0) {
            if ((t+1) % save_every == 0 || (t+1) == max_t) {
                info.save_time(t);
            }
        }

    }



    return;
}

