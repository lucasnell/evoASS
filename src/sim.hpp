#ifndef __EASS_SIM_H
#define __EASS_SIM_H


#include <RcppArmadillo.h>
#include "pcg.hpp"

using namespace Rcpp;



inline double trunc_rnorm__(const double& mu, const double& sigma) {

    double a_bar = (0 - mu) / sigma;

    double p = R::pnorm5(a_bar, 0, 1, 1, 0);
    double u = R::runif(p, 1);

    double x = R::qnorm5(u, 0, 1, 1, 0);
    x = x * sigma + mu;

    return x;
}


inline double trunc_rnorm__(const double& mu, const double& sigma, pcg32& eng) {

    double a_bar = (0 - mu) / sigma;

    double p = R::pnorm5(a_bar, 0, 1, 1, 0);
    double u = runif_ab(eng, p, 1);

    double x = R::qnorm5(u, 0, 1, 1, 0);
    x = x * sigma + mu;

    return x;
}



//' r (growth rate) based on Vi (Vi = traits for clone i).
//'
//' @noRd
//'
inline double r_V_(const arma::rowvec& Vi,
                   const double& f,
                   const arma::mat& C,
                   const double& r0) {

    double r = r0 - f * arma::as_scalar(Vi * C * Vi.t());

    return r;
}



//' A (interspecific + intraspecific density dependence) for line i based on V and N.
//'
//' This version is only used for `F_t_deriv_cpp` fxn.
//'
//' @noRd
//'
inline double A_i_VN_(const arma::rowvec V_i,
                      const std::vector<arma::rowvec>& V_nei,
                      const double& N_i,
                      const std::vector<double>& N_nei,
                      const double& g,
                      const double& d) {

    // Values of sum of squared trait values for each clone
    std::vector<double> W;
    W.reserve(V_nei.size() + 1);
    W.push_back(arma::as_scalar(V_i * V_i.t()));
    for (uint32_t j = 0; j < V_nei.size(); j++) {
        W.push_back(arma::as_scalar(V_nei[j] * V_nei[j].t()));
    }

    // Effects of intra- and inter-specific competition
    double A = g * std::exp(-1 * W[0]) * N_i; // intra
    for (uint32_t j = 0; j < V_nei.size(); j++) {
        A += (g * std::exp(-1 * (W[0] + d * W[j+1])) * N_nei[j]);
    }

    return A;
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
                  const std::vector<arma::rowvec>& V,
                  const std::vector<double>& N,
                  const double& g,
                  const double& d) {

    // Values of sum of squared trait values for each clone
    std::vector<double> W;
    W.reserve(V.size());
    for (uint32_t j = 0; j < V.size(); j++) {
        W.push_back(arma::as_scalar(V[j] * V[j].t()));
    }

    for (uint32_t i = 0; i < V.size(); i++) {
        // Effects of intra- and inter-specific competition
        double intra = g * std::exp(-1 * W[i]) * N[i];
        double inter = 0;
        for (uint32_t j = 0; j < V.size(); j++) {
            if (j == i) continue;
            inter += (g * std::exp(-1 * (W[i] + d * W[j])) * N[j]);
        }
        A[i] = intra + inter;
    }

    return;
}




//' Same as above, but using a vector of indices `I`.
//'
//' It's assumed higher-level functions will control the `A` vector.
//'
//' @noRd
//'
template <typename T>
inline void A_VNI_(T& A,
                   const std::vector<arma::rowvec>& V,
                   const std::vector<double>& N,
                   const std::vector<uint32_t>& I,
                   const double& g,
                   const double& d) {

    // Values of sum of squared trait values for each clone
    std::vector<double> W;
    W.reserve(I.size());
    for (uint32_t i = 0; i < I.size(); i++) {
        W.push_back(arma::as_scalar(V[I[i]] * V[I[i]].t()));
    }

    for (uint32_t i = 0; i < I.size(); i++) {

        // Effects of intra- and inter-specific competition
        double intra = g * std::exp(-1 * W[i]) * N[i];
        double inter = 0;
        for (uint32_t j = 0; j < I.size(); j++) {
            if (j == i) continue;
            inter += (g * std::exp(-1 * (W[i] + d * W[j])) * N[j]);
        }
        A[i] = intra + inter;
    }

    return;
}



//' Fitness at time t.
//'
//' @inheritParams F_t_cpp
//'
//' @noRd
//'
template <typename T>
inline void F_t__(T& F,
                  const std::vector<arma::rowvec>& V,
                  const std::vector<double>& N,
                  const double& f,
                  const double& g,
                  const arma::mat& C,
                  const double& r0,
                  const double& d) {

    std::vector<double> A(V.size());
    A_VN_<std::vector<double>>(A, V, N, g, d);

    for (uint32_t i = 0; i < V.size(); i++) {
        double r = r_V_(V[i], f, C, r0);
        F[i] = std::exp(r - A[i]);
    }

    return;
}

double F_t_deriv_cpp(const arma::rowvec V_i,
                     const std::vector<arma::rowvec>& V_nei,
                     const double& N_i,
                     const std::vector<double>& N_nei,
                     const double& f,
                     const double& g,
                     const arma::mat& C,
                     const double& r0,
                     const double& d);


#endif
