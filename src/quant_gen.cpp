
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


/*
 ******************
 ******************
 ******************

 Change below to output a vector of rows or a matrix.
 This will be more efficient than doing them individually.

 ******************
 ******************
 ******************
 */
void dF_dVi__(arma::rowvec& dF_dVi_vec,
              const arma::rowvec& V_i,
              const std::vector<arma::rowvec>& V_nei,
              const double& N_i,
              const std::vector<double>& N_nei,
              const double& f,
              const double& g,
              const arma::mat& C,
              const double& r0,
              const double& d) {

    // double F = F_t_deriv_cpp(V_i, V_nei, N_i, N_nei, f, g, C, r0, d);
    //
    // double tmp_dbl = 0;
    // for (uint32_t j = 0; j < V_nei.size(); j++) {
    //     const arma::rowvec& V_j(V_nei[j]);
    //     const double& N_j(N_nei[j]);
    //     tmp_dbl += (N_j / std::exp(d * arma::as_scalar(V_j * V_j.t())));
    // }
    //
    // dF_dVi_vec = (-f * V_i * (C + C.t())) +
    //     2 * V_i * g * std::exp(arma::as_scalar(-V_i * V_i.t())) * (N_i + tmp_dbl);
    //
    // dF_dVi_vec *= F;
    //
    // return;

    return;

}



//' Version of above for use in R for testing.
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

