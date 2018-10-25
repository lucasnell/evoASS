
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>


using namespace Rcpp;

inline double trunc_rnorm__(const double& mu, const double& sigma) {

    double a_bar = (0 - mu) / sigma;

    double p = R::pnorm5(a_bar, 0, 1, 1, 0);
    double u = R::runif(p, 1);

    double x = R::qnorm5(u, 0, 1, 1, 0);
    x = x * sigma + mu;

    return x;
}

//' Normal distribution truncated above zero.
//'
//' From `http://web.michaelchughes.com/research/sampling-from-truncated-normal`
//'
//' @noRd
//'
//[[Rcpp::export]]
double trunc_rnorm_(const double& mu, const double& sigma) {

    double x = trunc_rnorm__(mu, sigma);

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


//' A (interspecific + intraspecific density dependence) for line i based on V and N
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
    arma::vec W(V_nei.size() + 1);
    W(0) = arma::as_scalar(V_i * V_i.t());
    for (uint32_t j = 0; j < V_nei.size(); j++) {
        W(j+1) = arma::as_scalar(V_nei[j] * V_nei[j].t());
    }

    // Effects of intra- and inter-specific competition
    double A = g * std::exp(-1 * W(0)) * N_i; // intra
    for (uint32_t j = 0; j < V_nei.size(); j++) {
        A += (g * std::exp(-1 * (W(0) + d * W(j+1))) * N_nei[j]);
    }

    return A;
}

//' A (interspecific + intraspecific density dependence) based on V and N
//'
//' @noRd
//'
template <typename T>
inline T A_VN_(const std::vector<arma::rowvec>& V,
               const std::vector<double>& N,
               const double& g,
               const double& d) {

    T A(V.size());

    // Values of sum of squared trait values for each clone
    arma::vec W(V.size());
    for (uint32_t j = 0; j < V.size(); j++) {
        W(j) = arma::as_scalar(V[j] * V[j].t());
    }

    for (uint32_t i = 0; i < V.size(); i++) {

        // Effects of intra- and inter-specific competition
        double intra = g * std::exp(-1 * W(i)) * N[i];
        double inter = 0;
        for (uint32_t j = 0; j < V.size(); j++) {
            if (j == i) continue;
            inter += (g * std::exp(-1 * (W(i) + d * W(j))) * N[j]);
        }
        A[i] = intra + inter;
    }

    return A;
}

// Same as above, but fill an existing vector and using a vector of indices `I`
template <typename T>
inline void A_VNI_(T& A,
                   const std::vector<arma::rowvec>& V,
                   const std::vector<double>& N,
                   const std::vector<uint32_t>& I,
                   const double& g,
                   const double& d) {

    // Values of sum of squared trait values for each clone
    std::vector<double> W(I.size());
    for (uint32_t i = 0; i < I.size(); i++) {
        W[i] = arma::as_scalar(V[I[i]] * V[I[i]].t());
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



//' C++ implementation of the below function
//'
//'
//' @noRd
//'
inline void F_t__(arma::rowvec& F,
                  const std::vector<arma::rowvec>& V,
                  const std::vector<double>& N,
                  const double& f,
                  const double& g,
                  const arma::mat& C,
                  const double& r0,
                  const double& d) {

    arma::rowvec A = A_VN_<arma::rowvec>(V, N, g, d);

    for (uint32_t i = 0; i < V.size(); i++) {
        double r = r_V_(V[i], f, C, r0);
        F(i) = std::exp(r - A(i));
    }

    return;
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

//' Adaptive dynamics.
//'
//'
//' @export
//' @noRd
//'
//[[Rcpp::export]]
List adaptive_dynamics(
        const arma::rowvec& V0,
        const double& N0,
        const double& f = 0.1,
        const double& g = 0.1,
        const double& eta = 0.2,
        const double& r0 = 1,
        const double& d = -0.01,
        const double& max_t = 1e4,
        const double& min_N = 1e-4,
        const double& mut_sd = 0.1,
        const double& mut_prob = 0.01,
        const bool& show_progress = true,
        const uint32_t& max_clones = 1e4) {

    // # traits:
    uint32_t q = 2;

    if (V0.n_cols != q) stop("V0 must have exactly 2 cols");

    arma::mat C(q, q);
    C.fill(eta);
    C.diag().fill(1);



    /*
     -----------------
     Objects that change size over time
     -----------------
     */

    // Current abundances:
    std::vector<double> N;
    N.reserve(max_clones);
    N.push_back(N0);

    // Density dependences
    std::vector<double> A;
    A.reserve(max_clones);
    A.push_back(1); // 1 is simply a placeholder bc this will be recalculated later

    // Indices to trait values for all lines (this will get updated, though):
    std::vector<uint32_t> I;
    uint32_t clone_I = 0;
    I.reserve(max_clones);
    I.push_back(clone_I);
    clone_I++;



    /*
     -----------------
     Objects that store all information (and do not change size over time)
     -----------------
     */
    // Trait values for all clones:
    std::vector<arma::rowvec> all_V;
    all_V.reserve(max_clones);
    all_V.push_back(V0);

    // All abundances from t = 1 to max_t
    std::vector<std::vector<double>> all_N;
    all_N.reserve(max_t + 1);

    // All indices from t = 1 to max_t
    std::vector<std::vector<uint32_t>> all_I;
    all_I.reserve(max_t + 1);

    Progress p(max_t, show_progress);


    for (uint32_t t = 0; t < max_t; t++) {

        Rcpp::checkUserInterrupt();

        // Fill in density dependences:
        A_VNI_<std::vector<double>>(A, all_V, N, I, g, d);

        // Extinct clones (if any):
        std::vector<uint32_t> extinct;
        // Fill in abundances:
        for (uint32_t i = 0; i < A.size(); i++) {
            double r = r_V_(all_V[I[i]], f, C, r0);
            N[i] *= std::exp(r - A[i]);
            // See if it goes extinct:
            if (N[i] < min_N) {
                extinct.push_back(i);
            }
        }
        // If everything is gone, stop simulations:
        if (extinct.size() == N.size()) break;
        // Remove extinct clones (starting at the back):
        for (uint32_t i = 0, j; i < extinct.size(); i++) {
            j = extinct.size() - i - 1;
            N.erase(N.begin() + extinct[j]);
            A.erase(A.begin() + extinct[j]);
            I.erase(I.begin() + extinct[j]);
        }

        // Seeing if I should add new clones:
        uint32_t n_clones = N.size(); // doing this bc N.size() might change
        for (uint32_t i = 0; i < n_clones; i++) {

            double u = R::runif(0, 1);

            if (u < mut_prob) {

                N[i] -= (1.01 * min_N);

                N.push_back(1.01 * min_N);
                A.push_back(0);
                I.push_back(clone_I);
                clone_I++;

                all_V.push_back(arma::rowvec(q));
                arma::rowvec& new_V(all_V.back());
                arma::rowvec& old_V(all_V[I[i]]);
                for (uint32_t j = 0; j < q; j++) {
                    new_V(j) = trunc_rnorm__(old_V(j), mut_sd);
                }

            }

        }

        all_N.push_back(N);
        all_I.push_back(I);
        p.increment();

    }

    List out = List::create(
        _["N"] = all_N,
        _["V"] = all_V,
        _["I"] = all_I
    );

    return out;

}
