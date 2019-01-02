#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11, openmp)]]

using namespace Rcpp;


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



//' Creat Jacobian matrix for a particular set of species' traits.
//'
//'
//' @noRd
//'
void jacobian_(arma::mat& jacobian,
               const std::vector<arma::rowvec>& V,
               const std::vector<double>& N,
               const double& f,
               const double& g,
               const arma::mat& C,
               const double& r0,
               const double& d,
               const arma::vec& add_var,
               const double& eps) {

    uint32_t n = N.size();
    uint32_t q = V[0].n_elem;
    uint32_t nq = n * q;

    if (jacobian.n_rows != nq || jacobian.n_cols != nq) {
        jacobian.set_size(nq, nq);
    }

    arma::mat ss_mat;
    // Fill strength of selection matrix:
    sel_str__(ss_mat, V, N, f, g, C, r0, d);
    // Trait change based on current trait values:
    arma::rowvec FV(nq);
    for (uint32_t i = 0; i < n; i++) {
        FV(arma::span(i*q, (i+1)*q-1)) = V[i] + add_var(i) * ss_mat.row(i);
    }

    // Version of traits to manipulate:
    std::vector<arma::rowvec> V_m = V;
    // Trait change based on manipulated trait values:
    arma::rowvec FV_m(nq);

    uint32_t k = 0;
    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = 0; j < q; j++) {
            V_m[i](j) = V[i](j) + eps;
            sel_str__(ss_mat, V_m, N, f, g, C, r0, d);
            for (uint32_t ii = 0; ii < n; ii++) {
                FV_m(arma::span(ii*q, (ii+1)*q-1)) = V_m[ii] +
                    add_var(ii) * ss_mat.row(ii);
            }
            jacobian.row(k) = (FV_m - FV) / eps;
            V_m[i](j) = V[i](j);
            k++;
        }
    }

    return;

}

//' R-exported version of above, so it can be tested in R for accuracy.
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat jacobian_cpp(const std::vector<arma::rowvec>& V,
                       const std::vector<double>& N,
                       const double& f,
                       const double& g,
                       const arma::mat& C,
                       const double& r0,
                       const double& d,
                       const arma::vec& add_var,
                       const double& eps) {

    if (V.size() != N.size()) stop("V.size() != N.size()");
    if (add_var.n_elem != N.size()) stop("add_var.n_elem != N.size()");

    arma::mat jacobian;

    jacobian_(jacobian, V, N, f, g, C, r0, d, add_var, eps);

    return jacobian;
}
