
#include <RcppArmadillo.h>
#include "sim.hpp"

using namespace Rcpp;




//' Normal distribution truncated above zero.
//'
//' From `http://web.michaelchughes.com/research/sampling-from-truncated-normal`
//'
//' @noRd
//'
//[[Rcpp::export]]
std::vector<double> trunc_rnorm_cpp(const uint32_t& N,
                                    const double& mu,
                                    const double& sigma) {

    std::vector<double> out(N);
    for (double& x : out) x = trunc_rnorm_(mu, sigma);

    return out;
}


/*
 Same as above, but uses a vector of `mu`
 */
//[[Rcpp::export]]
std::vector<double> trunc_rnorm_mu_cpp(const std::vector<double>& mu,
                                       const double& sigma) {

    std::vector<double> out;
    out.reserve(mu.size());
    for (const double& x : mu) out.push_back(trunc_rnorm_(x, sigma));

    return out;
}
/*
 Same as above, but uses a vector of `sigma`
 */
//[[Rcpp::export]]
std::vector<double> trunc_rnorm_sigma_cpp(const double& mu,
                                          const std::vector<double>& sigma) {

    std::vector<double> out;
    out.reserve(sigma.size());
    for (const double& x : sigma) out.push_back(trunc_rnorm_(mu, x));

    return out;
}
/*
 Same as above, but uses a vector of `mu` and `sigma`
 */
//[[Rcpp::export]]
std::vector<double> trunc_rnorm_mu_sigma_cpp(const std::vector<double>& mu,
                                             const std::vector<double>& sigma) {

    uint32_t n = mu.size();
    if (sigma.size() != n) stop("sigma.size() != mu.size()");

    std::vector<double> out;
    out.reserve(n);
    for (uint32_t i = 0; i < n; i++) {
        out.push_back(trunc_rnorm_(mu[i], sigma[i]));
    }

    return out;
}



/*
 Fitness at time t for all species
 */
//[[Rcpp::export]]
arma::vec F_t_cpp(const arma::mat& V,
                  const std::vector<double>& N,
                  const double& f,
                  const double& a0,
                  const arma::mat& C,
                  const double& r0,
                  const arma::mat& D) {

    uint32_t n = V.n_cols;

    if (n != N.size()) stop("V.n_cols != N.size()");

    arma::vec F(n);
    std::vector<arma::vec> VV;
    VV.reserve(n);
    for (uint32_t i = 0; i < n; i++) VV.push_back(V.col(i));

    F_t__<arma::vec>(F, VV, N, f, a0, C, r0, D);

    return F;
}

/*
 Fitness at time t for species i.
 */
//[[Rcpp::export]]
double F_it_cpp(const uint32_t& i,
                const arma::mat& V,
                const std::vector<double>& N,
                const double& f,
                const double& a0,
                const arma::mat& C,
                const double& r0,
                const arma::mat& D) {

    if (V.n_cols != N.size()) stop("V.n_cols != N.size()");

    double F = F_it__(i, V, N, f, a0, C, r0, D);

    return F;
}




//[[Rcpp::export]]
bool using_openmp() {
    bool out = false;
#ifdef _OPENMP
    out = true;
#endif
    return out;
}


