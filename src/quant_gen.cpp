
// #define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <random>

#include "sim.hpp"
#include "pcg.hpp"
#include "quant_gen.hpp"

#ifdef _OPENMP
#include <omp.h>  // omp
#endif

using namespace Rcpp;

/*
 Note that this is the "instantaneous" version of quantitative genetics.
 It's not the same as what's written in Tobin and Tony's methods section.
 Instead of adding new daughter species, it starts with a number of species
 that should be more than you would ever expect to be produced sequentially.
 */



//' One round of quantitative genetics.
//'
//' Higher-up function(s) should create `C` from `eta` and should handle any outputs.
//'
//'
//' @noRd
//'
void one_quantgen_rep(OneRepInfo& info,
                      const std::vector<arma::rowvec>& V0,
                      const std::vector<double>& N0,
                      const double& f,
                      const double& g,
                      const arma::mat& C,
                      const double& r0,
                      const double& d,
                      const double& add_var,
                      const double& delta,
                      const double& delta2,
                      const uint32_t& start_t,
                      const uint32_t& max_t,
                      const double& min_N,
                      const uint32_t& save_every,
                      pcg32& eng) {


    uint32_t n = N0.size();     // # species
    uint32_t q = V0[0].n_elem;  // # traits

    if (save_every > 0) {
        info = OneRepInfo(N0, V0, max_t, save_every);
    } else {
        info = OneRepInfo(N0, V0);
    }
    std::vector<arma::rowvec>& V(info.V);
    std::vector<double>& N(info.N);
    std::vector<double> A(V.size());

    for (uint32_t t = 0; t < start_t; t++) {

        /*
         LEFT OFF --> TURN THIS INTO YOUR VERSION OF THE QUANTITATIVE GENETICS
         */

        // Fill in density dependences:
        A_VN_<std::vector<double>>(A, V, N, g, d);
        // Fill in abundances:
        for (uint32_t i = 0; i < A.size(); i++) {
            double r = r_V_(V[i], f, C, r0);
            N[i] *= std::exp(r - A[i]);
        }
        // // attack rates of consumer-->resources (<# resources> x <# consumers>):
        // B = b * arma::exp(-(V % V) * (U.t() % U.t()));
        // // consumer mortality rate  (<# consumers> x 1):
        // M = m + g / arma::diagvec(U * D * U.t());
        // // New resources' abundances (<# resources> x 1)
        // Nt = N % arma::exp(r * (1 - A * arma::accu(N) - B * P));
        // // New consumers' abundances (<# consumers> x 1)
        // Pt = P % arma::exp(c * B.t() * N - M);
        //
        // // Traits for each resource (<# resources> x <# resources' traits>):
        // set_Vt_(Vt, V, U, N, P, C, r, f, b, sig2N, p, q);
        // // Traits for each consumer (<# consumers> x <# consumers' traits>)
        // set_Ut_(Ut, V, U, N, D, b, c, g, sig2P, n, p, q);
        //
        // // Setting new values:
        // N = Nt;
        // P = Pt;
        // V = arma::abs(Vt);
        // U = arma::abs(Ut);

    }

    // perturbation
    for (uint32_t j = 0; j < q; j++) {
        for (uint32_t i = 0; i < n; i++) {
            V(i, j) *= distr(eng);
            // V(i, j) *= R::rlnorm(0, delta);
        }
        for (uint32_t i = 0; i < q; i++) {
            U(i, j) *= distr(eng);
            // U(i, j) *= R::rlnorm(0, delta);
        }
    }



    for (uint32_t t = 0; t < max_t; t++) {

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

        if ((t+1) % save_every == 0 || (t+1) == max_t) {  // <-- use this, not quant_gen way
            // all_N.push_back(N);
            // all_I.push_back(I);
            // all_t.push_back(t);
        }

    }



    return;
}

