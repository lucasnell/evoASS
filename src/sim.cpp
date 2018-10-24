
#include <RcppArmadillo.h>

using namespace Rcpp;



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


//' A (interspecific + intraspecific density dependence) for line i based on V and N
//'
//' @noRd
//'
inline double A_i_VN_(const arma::rowvec V_i,
                      const std::vector<arma::rowvec>& V_nei,
                      const double& N_i,
                      const arma::rowvec& N_nei,
                      const double& g,
                      const double& d) {

    // Values of sum of squared trait values for each clone
    arma::vec W(V_nei.size() + 1);
    W(0) = arma::as_scalar(V_i * V_i.t());
    for (arma::uword j = 0; j < V_nei.size(); j++) {
        W(j+1) = arma::as_scalar(V_nei[j] * V_nei[j].t());
    }

    // Effects of intra- and inter-specific competition
    double A = g * std::exp(-1 * W(0)) * N_i; // intra
    for (arma::uword j = 0; j < V_nei.size(); j++) {
        A += (g * std::exp(-1 * (W(0) + d * W(j+1))) * N_nei(j));
    }

    return A;
}

//' A (interspecific + intraspecific density dependence) based on V and N
//'
//' @noRd
//'
inline arma::rowvec A_VN_(const std::vector<arma::rowvec>& V,
                          const arma::rowvec& N,
                          const double& g,
                          const double& d) {

    arma::rowvec A(V.size());

    // Values of sum of squared trait values for each clone
    arma::vec W(V.size());
    for (arma::uword j = 0; j < V.size(); j++) {
        W(j) = arma::as_scalar(V[j] * V[j].t());
    }

    for (arma::uword i = 0; i < V.size(); i++) {

        // Effects of intra- and inter-specific competition
        double intra = g * std::exp(-1 * W(i)) * N(i);
        double inter = 0;
        for (arma::uword j = 0; j < V.size(); j++) {
            if (j == i) continue;
            inter += (g * std::exp(-1 * (W(i) + d * W(j))) * N(j));
        }
        A(i) = intra + inter;
    }

    return A;
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
                  const arma::rowvec& N,
                  const double& f,
                  const double& g,
                  const arma::mat& C,
                  const double& r0,
                  const double& d) {

    arma::rowvec A = A_VN_(V, N, g, d);
    arma::rowvec F(V.size());

    for (arma::uword i = 0; i < V.size(); i++) {
        double r = r_V_(V[i], f, C, r0);
        F(i) = std::exp(r - A(i));
    }

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
                  const arma::rowvec& N_nei,
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
