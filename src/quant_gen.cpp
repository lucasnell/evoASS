#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>

#include "sim.hpp"
#include "quant_gen.hpp"


using namespace Rcpp;



//' One round of iterations with no mutations to new species.
//'
//' Higher-up function(s) should create `C` from `eta` and should handle any outputs
//' and mutations to new species.
//'
//'
//' @noRd
//'
void one_round(std::vector<arma::rowvec>& V,
               std::vector<double>& N,
               const std::vector<arma::rowvec>& V0,
               const std::vector<double>& N0,
               const double& f,
               const double& g,
               const arma::mat& C,
               const double& r0,
               const double& d,
               const double& max_t,
               const double& min_N) {

    V = V0;
    N = N0;
    std::vector<double> A(V.size());

    for (uint32_t t = 0; t < max_t; t++) {

        Rcpp::checkUserInterrupt();

        // Fill in density dependences:
        A_VN_<std::vector<double>>(A, V, N, g, d);

        // Extinct clones (if any):
        std::vector<uint32_t> extinct;
        // Fill in abundances:
        for (uint32_t i = 0; i < A.size(); i++) {
            double r = r_V_(V[i], f, C, r0);
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
            V.erase(V.begin() + extinct[j]);
            A.erase(A.begin() + extinct[j]);
        }

    }



    return;
}

