
#include <RcppArmadillo.h>
#include <random>
#include <progress.hpp>
#include <progress_bar.hpp>

#include "sim.hpp"
#include "adapt_dyn.hpp"


using namespace Rcpp;




//' Multiple repetitions of adaptive dynamics.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat adapt_dyn_cpp(const uint32_t& n_reps,
                        const std::vector<arma::rowvec>& V0,
                        const std::vector<double>& N0,
                        const double& f,
                        const double& a0,
                        const double& eta,
                        const double& r0,
                        const double& d,
                        const double& max_t,
                        const double& min_N,
                        const double& mut_sd,
                        const double& mut_prob,
                        const bool& show_progress,
                        const uint32_t& max_clones,
                        const uint32_t& save_every,
                        const uint32_t& n_cores) {

    if (V0.size() == 0) stop("empty V0 vector");
    if (V0[0].n_elem == 0) stop("empty V0[0] vector");
    if (N0.size() == 0) stop("empty N0 vector");
    if (V0.size() != N0.size()) stop("V0 and N0 must be the same size");

    // # traits:
    uint32_t q = V0[0].n_elem;

    for (uint32_t i = 0; i < V0.size(); i++) {
        if (V0[i].n_cols != q) stop("all items in V0 must have the same length / # cols");
    }

    arma::mat C(q, q);
    C.fill(eta);
    C.diag().fill(1);

    std::vector<OneRepInfoAD> rep_infos(n_reps);

    const std::vector<std::vector<uint64_t>> seeds = mc_seeds(n_cores);
    Progress prog_bar(n_reps, show_progress);
    bool interrupted = false;

    #ifdef _OPENMP
    #pragma omp parallel default(shared) num_threads(n_cores) if (n_cores > 1)
    {
    #endif

    std::vector<uint64_t> active_seeds;

    // Write the active seed per core or just write one of the seeds.
    #ifdef _OPENMP
    uint32_t active_thread = omp_get_thread_num();
    #else
    uint32_t active_thread = 0;
    #endif
    active_seeds = seeds[active_thread];

    pcg64 eng = seeded_pcg(active_seeds);

    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (uint32_t i = 0; i < n_reps; i++) {

        if (!Progress::check_abort()) {

            rep_infos[i] = OneRepInfoAD(V0, N0, max_clones, max_t, save_every, mut_sd);

            for (uint32_t t = 0; t < max_t; t++) {

                rep_infos[i].iterate(t, f, a0, C, r0, d, max_t, min_N, mut_sd, mut_prob,
                                     save_every, eng);

            }

            prog_bar.increment();

        } else if (active_thread == 0) interrupted = true;
    }
    #ifdef _OPENMP
    }
    #endif

    if (interrupted) {
        throw(Rcpp::exception("\nUser interrupted process.", false));
    }

    /*
     Go through one time to calculate the # surviving species for all reps and
     for all time point(s) saved.
     Then go back through and fill in values.
    */
    uint32_t n_rows = rep_infos[0].n_rows();
    std::vector<uint32_t> cum_rows(n_reps, 0);
    cum_rows[0] = n_rows;
    for (uint32_t i = 1; i < n_reps; i++) {
        uint32_t nr = rep_infos[i].n_rows();
        n_rows += nr;
        cum_rows[i] = nr + cum_rows[i-1];
    }
    // Make output matrix:
    arma::mat output(n_rows, 4 + q);
    // Fill output matrix:
    #ifdef _OPENMP
    #pragma omp parallel for default(shared) num_threads(n_cores) schedule(static)
    #endif
    for (uint32_t i = 0; i < n_reps; i++) {
        uint32_t start = 0;
        if (i > 0) start = cum_rows[i-1];
        rep_infos[i].fill_matrix(output, i, start);
    }




    return output;
}
