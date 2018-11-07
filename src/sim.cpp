
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
std::vector<double> trunc_rnorm_(const uint32_t& N,
                                 const double& mu,
                                 const double& sigma) {

    std::vector<double> out(N);
    for (double& x : out) x = trunc_rnorm__(mu, sigma);

    return out;
}




//' Fitness at time t.
//'
//' @param V List of row vectors, each containing trait values at time t (1x2 vector)
//'     for a particular clone.
//' @param N Row vector of population abundances at time t.
//' @param f Effect of traits on growth rate.
//' @param g Effect of traits on density dependence.
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
arma::rowvec F_t_(const std::vector<arma::rowvec>& V,
                  const std::vector<double>& N,
                  const double& f,
                  const double& g,
                  const arma::mat& C,
                  const double& r0,
                  const double& d) {

    arma::rowvec F(V.size());

    F_t__(F, V, N, f, g, C, r0, d);

    return F;
}


//' Same as above, but for only one clone.
//'
//' This is used for computing partial derivatives.
//'
//' @param V_i Row vector of traits for focal clone.
//' @param V_nei List of row vectors of traits for non-focal clones.
//' @param N_i Abundance for focal clone.
//' @param N_nei Row vector of abundances for non-focal clones.
//' @inheritParams F_t_
//'
//' @noRd
//'
//[[Rcpp::export]]
double F_t_deriv_(const arma::rowvec V_i,
                  const std::vector<arma::rowvec>& V_nei,
                  const double& N_i,
                  const std::vector<double>& N_nei,
                  const double& f,
                  const double& g,
                  const arma::mat& C,
                  const double& r0,
                  const double& d) {

    double A = A_i_VN_(V_i, V_nei, N_i, N_nei, g, d);
    double r = r_V_(V_i, f, C, r0);
    double F = std::exp(r - A);
    return F;
}


