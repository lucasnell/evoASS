
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>

#include "sim.hpp"


using namespace Rcpp;



//' Adaptive dynamics.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
List adaptive_dynamics_(
        const std::vector<arma::rowvec>& V0,
        const std::vector<double>& N0,
        const double& f,
        const double& g,
        const double& eta,
        const double& r0,
        const double& d,
        const double& max_t,
        const double& min_N,
        const double& mut_sd,
        const double& mut_prob,
        const bool& show_progress,
        const uint32_t& max_clones,
        const uint32_t& save_every) {

    // # traits:
    uint32_t q = 2;

    for (uint32_t i = 0; i < V0.size(); i++) {
        if (V0[i].n_cols != q) stop("all rows in V0 must have exactly 2 cols");
    }
    if (V0.size() != N0.size()) {
        stop("V0 and N0 must be the same size");
    }

    arma::mat C(q, q);
    C.fill(eta);
    C.diag().fill(1);



    /*
     -----------------
     Objects that change size over time
     -----------------
     */

    // Current abundances:
    std::vector<double> N = N0;
    N.reserve(max_clones);

    // Density dependences
    std::vector<double> A(N0.size());
    A.reserve(max_clones);

    // Indices to trait values for all lines (this will get updated, though):
    std::vector<uint32_t> I;
    uint32_t clone_I = 0;
    I.reserve(max_clones);
    while (clone_I < N0.size()) {
        I.push_back(clone_I);
        clone_I++;
    }



    /*
     -----------------
     Objects that store all information (and do not change size over time)
     -----------------
     */
    // Trait values for all clones:
    std::vector<arma::rowvec> all_V = V0;
    all_V.reserve(V0.size() + max_clones);

    // All abundances from t = 1 to max_t
    std::vector<std::vector<double>> all_N;
    all_N.reserve((max_t / save_every) + 2);
    all_N.push_back(N);

    // All indices from t = 1 to max_t
    std::vector<std::vector<uint32_t>> all_I;
    all_I.reserve((max_t / save_every) + 2);
    all_I.push_back(I);

    // All time points sampled
    std::vector<uint32_t> all_t;
    all_t.reserve((max_t / save_every) + 2);
    all_t.push_back(0);


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
        // If everything is gone, output results and stop simulations:
        if (extinct.size() == N.size()) {
            all_N.push_back(N);
            all_I.push_back(I);
            all_t.push_back(t);
            break;
        }
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
                    // new_V(j) = trunc_rnorm__(old_V(j), mut_sd);
                    new_V(j) = R::rnorm(old_V(j), mut_sd);
                }

            }

        }

        if ((t+1) % save_every == 0 || (t+1) == max_t) {
            all_N.push_back(N);
            all_I.push_back(I);
            all_t.push_back(t);
        }

        p.increment();

    }

    List out = List::create(
        _["N"] = all_N,
        _["V"] = all_V,
        _["I"] = all_I,
        _["T"] = all_t
    );

    return out;

}
