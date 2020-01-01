
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




//' Fitness at time t.
//'
//' @param V List of row vectors, each containing trait values at time t (1x2 vector)
//'     for a particular clone.
//' @param N Row vector of population abundances at time t.
//' @param f Effect of traits on growth rate.
//' @param a0 Base density dependence.
//' @param C Matrix containing non-additive effects of traits on growth rate.
//' @param r0 Starting growth rate.
//' @param d Changes how the focal line is affected by other lines' trait values.
//'     If `d < 0`, then increases in `V_j` (trait that reduces competition
//'     experienced by clone `j`) increases competition experienced by clone `i`,
//'     thereby giving conflicting coevolution.
//'     Conversely, if `d > 0`, then increases in `V_j` decrease competition
//'     experienced by clone `i`, leading to nonconflicting coevolution.
//'
//'
//' @export
//'
//[[Rcpp::export]]
arma::rowvec F_t_cpp(const std::vector<arma::rowvec>& V,
                     const std::vector<double>& N,
                     const double& f,
                     const double& a0,
                     const arma::mat& C,
                     const double& r0,
                     const arma::mat& D) {

    arma::rowvec F(V.size());

    F_t__<arma::rowvec>(F, V, N, f, a0, C, r0, D);

    return F;
}

