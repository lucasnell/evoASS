#ifndef __SAURON_SIM_H
#define __SAURON_SIM_H


// #define ARMA_NO_DEBUG



#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "pcg.hpp"

using namespace Rcpp;




// For checking for user interrupts every N iterations:
inline bool interrupt_check(uint32_t& iters,
                            Progress& prog_bar,
                            const uint32_t& N) {
    ++iters;
    if (iters > N) {
        if (prog_bar.is_aborted() || prog_bar.check_abort()) return true;
        iters = 0;
    }
    return false;
}


//' Normal distribution truncated above zero.
//'
//' @noRd
//'
inline double trunc_rnorm_(const double& mu, const double& sigma, pcg64& eng) {

    double a_bar = (0 - mu) / sigma;

    double p = R::pnorm5(a_bar, 0, 1, 1, 0);
    double u = runif_ab(eng, p, 1);

    double x = R::qnorm5(u, 0, 1, 1, 0);
    x = x * sigma + mu;

    return x;
}

//' Same as above, but using R's RNG.
//'
//' Used in `trunc_rnorm_cpp` only, for testing.
//'
//' @noRd
//'
inline double trunc_rnorm_(const double& mu, const double& sigma) {

    double a_bar = (0 - mu) / sigma;

    double p = R::pnorm5(a_bar, 0, 1, 1, 0);
    double u = R::runif(p, 1);

    double x = R::qnorm5(u, 0, 1, 1, 0);
    x = x * sigma + mu;

    return x;
}



//' r (growth rate) based on Vi (Vi = traits for clone i).
//'
//' @noRd
//'
inline double r_V_(const arma::vec& Vi,
                   const double& f,
                   const arma::mat& C,
                   const double& r0) {

    double r = r0 - f * arma::as_scalar(Vi.t() * C * Vi);

    return r;
}




//' A (interspecific + intraspecific density dependence) based on V and N.
//'
//'
//' It's assumed higher-level functions will control the `A` vector size!
//'
//' @noRd
//'
template <typename T>
inline void A_VN_(T& A,
                  const std::vector<arma::vec>& V,
                  const std::vector<double>& N,
                  const double& a0,
                  const arma::mat& D) {

    // Values of sum of squared trait values for each clone
    std::vector<double> W_intra; // `exp(- t(V) %*% V)` for intraspecific component
    std::vector<double> W_inter; // `exp(- t(V) %*% D %*% V)` for interspecific component
    W_intra.reserve(V.size());
    W_inter.reserve(V.size());
    for (uint32_t j = 0; j < V.size(); j++) {
        W_intra.push_back(std::exp(-1 * arma::as_scalar(V[j].t() * V[j])));
        W_inter.push_back(std::exp(-1 * arma::as_scalar(V[j].t() * D * V[j])));
    }

    for (uint32_t i = 0; i < V.size(); i++) {
        // Effects of intra- and inter-specific competition
        double Omega = N[i];
        for (uint32_t j = 0; j < V.size(); j++) {
            if (j == i) continue;
            Omega += N[j] * W_inter[j];
        }
        A[i] = a0 * W_intra[i] * Omega;
    }

    return;
}




//' Same as above, but using a vector of indices `I`.
//'
//' It's assumed higher-level functions will control the `A` vector.
//' Used in `adapt_dyn_cpp` function.
//'
//' @noRd
//'
template <typename T>
inline void A_VNI__(T& A,
                    const std::vector<arma::vec>& V,
                    const std::vector<double>& N,
                    const std::vector<uint32_t>& I,
                    const double& a0,
                    const arma::mat& D) {

    // Values of sum of squared trait values for each clone
    std::vector<double> W_intra; // `exp(- t(V) %*% V)` for intraspecific component
    std::vector<double> W_inter; // `exp(- t(V) %*% D %*% V)` for interspecific component
    W_intra.reserve(I.size());
    W_inter.reserve(I.size());
    for (uint32_t j = 0; j < I.size(); j++) {
        W_intra.push_back(std::exp(-1 * arma::as_scalar(V[I[j]].t() * V[I[j]])));
        W_inter.push_back(std::exp(-1 * arma::as_scalar(V[I[j]].t() * D * V[I[j]])));
    }

    for (uint32_t i = 0; i < I.size(); i++) {
        // Effects of intra- and inter-specific competition
        double Omega = N[i];
        for (uint32_t j = 0; j < V.size(); j++) {
            if (j == i) continue;
            Omega += N[j] * W_inter[j];
        }
        A[i] = a0 * W_intra[i] * Omega;
    }

    return;
}


/*
 Fitness at time t.
*/
template <typename T>
inline void F_t__(T& F,
                  const std::vector<arma::vec>& V,
                  const std::vector<double>& N,
                  const double& f,
                  const double& a0,
                  const arma::mat& C,
                  const double& r0,
                  const arma::mat& D) {

    std::vector<double> A(V.size());
    A_VN_<std::vector<double>>(A, V, N, a0, D);

    for (uint32_t i = 0; i < V.size(); i++) {
        double r = r_V_(V[i], f, C, r0);
        F[i] = std::exp(r - A[i]);
    }

    return;
}



#endif
