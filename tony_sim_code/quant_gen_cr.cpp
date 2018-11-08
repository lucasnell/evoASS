/*
 --------------------------------------------------
 --------------------------------------------------

 Quantitative genetics for consumer-resource model

 --------------------------------------------------
 --------------------------------------------------
 */

#include <RcppArmadillo.h>


using namespace Rcpp;


//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(openmp, cpp11)]]


//' One round of iterations with no mutations to new species.
//'
//' Higher-up function(s) should create `C` from `eta` and should handle any outputs
//' and mutations to new species.
//'
//'
//' @noRd
//'
void one_round(std::vector<double>& N,
               std::vector<double>& P,
               std::vector<arma::rowvec>& V,
               std::vector<arma::rowvec>& U,
               const std::vector<double>& N0,
               const std::vector<double>& P0,
               const std::vector<arma::rowvec>& V0,
               const std::vector<arma::rowvec>& U0,
               const double& r,
               const double& a,
               const double& f,
               const double& b,
               const double& c,
               const double& m,
               const double& g,
               const double& etaN,
               const double& etaP,
               const double& d,
               const double& start_t,
               const double& max_t,
               const double& min_N) {

    N = N0;
    P = P0;
    V = V0;
    U = U0;

    uint32_t n = N.size();      // # resources
    uint32_t p = P.size();      // # consumers
    uint32_t q = V0[0].n_elem;  // # traits

    arma::mat C(q, q);
    C.fill(etaN);
    C.diag().fill(1);
    arma::mat D(q, q);
    D.fill(etaP);
    D.diag().fill(1);

    for (uint32_t t = 0; t < start_t; t++) {

        Rcpp::checkUserInterrupt();

        // for (time in 1:starttime) {
        //
        //     A <- a + f * diag(V %*% C %*% t(V))
        //     B <- b * exp(-(V^2) %*% t(U^2))
        //     M <- m + g/diag(U %*% D %*% t(U))
        //
        //     Nt <- N * exp(r * (1 - A * sum(N) - B %*% P))
        //     Pt <- P * exp(cc * t(B) %*% N - M)
        //     Vt <- V + sig2N * (-2 * r * f * sum(N) * V %*% C + 2 * r * b * V *
        //         exp(-(V^2) %*% t(U^2)) %*% (U^2 * array(P, c(p, q))))
        //     Ut <- U + sig2P * (-2 * b * cc * U * exp(-(U^2) %*% t(V^2)) %*%
        //         (V^2 * array(N, c(n, q))) + 2 * g *
        //         array(1/diag(U %*% D %*% t(U))^2, c(p, q)) * (U %*% D))
        //
        //     N <- Nt
        //     P <- Pt
        //     V <- abs(Vt)
        //     U <- abs(Ut)
        // }

    }



    return;
}

