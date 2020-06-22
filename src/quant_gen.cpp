

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




//' Derivative of fitness with respect to the trait divided by mean fitness.
//'
//' The function below calculates this selection strength for all traits
//' for all species.
//'
//' @noRd
//'
inline void sel_str__(arma::mat& ss_mat,
                      const std::vector<arma::rowvec>& V,
                      const std::vector<double>& N,
                      const double& f,
                      const double& a0,
                      const arma::mat& C,
                      const double& r0,
                      const arma::mat& D) {

    uint32_t n_spp = V.size();      // # species
    uint32_t n_trt = V[0].size();   // # traits

    if (ss_mat.n_rows != n_spp || ss_mat.n_cols != n_trt) ss_mat.set_size(n_spp, n_trt);

    // For all i, calculate `N_j * exp(- V_j * D * transpose(V_j))`, where `j != i`
    arma::vec W(n_spp, arma::fill::zeros);
    for (uint32_t j = 0; j < n_spp; j++) {
        // Calculate `N_j * exp(-V_j * D * transpose(V_j))`:
        double W_ = N[j] * std::exp(-1 * arma::as_scalar(V[j] * D * V[j].t()));
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
            2 * a0 * V_i * std::exp(arma::as_scalar(-V_i * V_i.t())) * (N_i + W(i));
    }

    return;

}



//' R-exported version of above, so it can be tested in R for accuracy.
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat sel_str_cpp(const std::vector<arma::rowvec>& V,
                      const std::vector<double>& N,
                      const double& f,
                      const double& a0,
                      const arma::mat& C,
                      const double& r0,
                      const arma::mat& D) {

    arma::mat ss_mat;

    sel_str__(ss_mat, V, N, f, a0, C, r0, D);

    return ss_mat;
}



//' Partial derivative of species i traits at time t+1 with respect to species i traits
//' at time t.
//'
//'
//' @noRd
//'
inline void dVi_dVi_(arma::mat& dVhat,
                     const uint32_t& row_start,
                     const uint32_t& col_start,
                     const arma::rowvec& Vi,
                     const double& Z,
                     const arma::mat& C,
                     const double& f,
                     const double& a0,
                     const double& add_var) {
    uint32_t q = Vi.n_elem;
    arma::mat I = arma::eye<arma::mat>(q, q);
    uint32_t row_end = row_start + q - 1;
    uint32_t col_end = col_start + q - 1;
    dVhat(arma::span(row_start, row_end), arma::span(col_start, col_end)) = I +
        2 * add_var * (
                (
                        a0 * Z * std::exp(-1 * arma::as_scalar(Vi * Vi.t())) *
                            (I - 2 * (Vi.t() * Vi))
                ) - (f * C)
        );
    return;
}


//' R-exported version of above, to be used in R for testing.
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat dVi_dVi_cpp(const uint32_t& i, const arma::mat& V, const double& Z,
                      const arma::mat& C, const double& f, const double& a0,
                      const double& add_var) {
    uint32_t q = V.n_cols;
    arma::mat dVhat(q, q);

    // Fill dVhat:
    dVi_dVi_(dVhat, 0, 0, V.row(i), Z, C, f, a0, add_var);

    dVhat = dVhat.t();

    return dVhat;
}


//' Partial derivative of species i traits at time t+1 with respect to species k traits
//' at time t.
//'
//' @noRd
//'
//'
inline void dVi_dVk_(arma::mat& dVhat,
                     const uint32_t& row_start,
                     const uint32_t& col_start,
                     const double& Nk,
                     const arma::rowvec& Vi,
                     const arma::rowvec& Vk,
                     const arma::mat& D,
                     const double& a0,
                     const double& add_var) {
    uint32_t row_end = row_start + Vi.n_elem - 1;
    uint32_t col_end = col_start + Vi.n_elem - 1;
    arma::mat M = D * Vk.t() * arma::exp(-1 * Vk * D * Vk.t() - Vi * Vi.t()) * Vi;
    M *= (-4 * a0 * add_var * Nk);
    dVhat(arma::span(row_start, row_end), arma::span(col_start, col_end)) = M;
    return;
}


//' R-exported version of above, to be used in R for testing.
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat dVi_dVk_cpp(const uint32_t& i,
                      const uint32_t& k,
                      const std::vector<double>& N,
                      const arma::mat& V,
                      const arma::mat& D,
                      const double& a0,
                      const double& add_var) {
    uint32_t q = V.n_cols;
    arma::mat dVhat(q, q);
    // Fill dVhat:
    dVi_dVk_(dVhat, 0, 0, N[k], V.row(i), V.row(k), D, a0, add_var);

    dVhat = dVhat.t();

    return dVhat;
}







//' Calculate the Jacobian of first derivatives.
//'
//' Cell [i,j] contains the partial derivative of j with respect to i.
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::mat jacobian_cpp(const std::vector<arma::rowvec>& V,
                        const std::vector<double>& N,
                        const double& f,
                        const double& a0,
                        const arma::mat& D,
                        const arma::mat& C,
                        const arma::vec& add_var) {

    uint32_t n = N.size();
    uint32_t q = V[0].n_cols;

    if (V.size() != n) stop("V.size() != N.size()");
    if (add_var.n_elem != n) stop("add_var.n_elem != N.size()");
    for (uint32_t i = 0; i < n; i++) {
        if (V[i].n_elem != q) stop("V[i].n_cols != q");
    }

    arma::mat jcb_mat(n*q, n*q);

    arma::vec Z_vec(n);
    for (uint32_t j = 0; j < n; j++) {
        Z_vec[j] = (N[j] * std::exp(arma::as_scalar(-1 * V[j] * D * V[j].t())));
    }

    for (uint32_t i = 0; i < n; i++) {

        const arma::rowvec& Vi(V[i]);
        const double& add_var_i(add_var[i]);
        uint32_t col_start = i * q;

        for (uint32_t k = 0; k < n; k++) {

            uint32_t row_start = k * q;

            if (k == i) {

                double Z = N[i];
                for (uint32_t j = 0; j < n; j++) {
                    if (j != i) Z += Z_vec[j];
                }
                // Fill Jacobian:
                dVi_dVi_(jcb_mat, row_start, col_start, Vi, Z, C, f, a0,
                         add_var_i);

            } else {

                // Fill Jacobian:
                dVi_dVk_(jcb_mat, row_start, col_start, N[k], Vi, V[k], D, a0,
                         add_var_i);

            }

        }

    }

    jcb_mat = jcb_mat.t();


    return jcb_mat;
}





//' Search for unique species in a matrix of species trait values.
//'
//' @noRd
//'
inline void unq_spp_(arma::uvec& is_unq,
                     const std::vector<arma::rowvec>& V,
                     const double& precision) {

    double precision_ = precision * precision;

    uint32_t n = V.size();
    is_unq = arma::ones<arma::uvec>(n);

    for (uint32_t i = 1; i < n; i++) {
        for (uint32_t j = 0; j < i; j++) {
            if (is_unq(j) == 0) continue; // don't want to keep looking at non-unique spp
            arma::rowvec diff_;
            diff_ = (V[i] - V[j]) % (V[i] - V[j]);
            if (arma::mean(diff_) < precision_) {
                is_unq(i) = 0;
                break;
            }
        }
    }

    return;
}

//' Same as above, but exported for use in R
//'
//' @noRd
//'
//[[Rcpp::export]]
arma::uvec unq_spp_cpp(const std::vector<arma::rowvec>& V,
                       const double& precision) {

    arma::uvec is_unq;

    unq_spp_(is_unq, V, precision);

    return is_unq;
}


/*
 Similar to above, except that it groups species in to unique groups.
 */
//[[Rcpp::export]]
IntegerVector group_spp_cpp(const std::vector<arma::rowvec>& V,
                                    double precision) {


    precision *= precision;

    uint32_t n = V.size();

    IntegerVector spp_groups(n, 0);

    arma::uvec is_unq = arma::ones<arma::uvec>(n);

    uint32_t g = 1;

    for (uint32_t i = 1; i < n; i++) {
        for (uint32_t j = 0; j < i; j++) {
            if (is_unq(j) == 0) continue; // don't want to keep looking at non-unique spp
            arma::rowvec diff_;
            diff_ = (V[i] - V[j]) % (V[i] - V[j]);
            if (arma::mean(diff_) < precision) {
                is_unq(i) = 0;
                spp_groups[i] = spp_groups[j];
                break;
            }
        }
        if (is_unq(i) == 1) {
            spp_groups[i] = g;
            g++;
        }
    }

    return spp_groups;
}




//' One repetition of quantitative genetics.
//'
//' Higher-up function(s) should handle the info put into `info`.
//'
//'
//' @noRd
//'
int one_quant_gen__(OneRepInfo& info,
                     const std::vector<arma::rowvec>& V0,
                     const std::vector<double>& N0,
                     const double& f,
                     const double& a0,
                     const arma::mat& C,
                     const double& r0,
                     const arma::mat& D,
                     const arma::vec& add_var,
                     const double& perturb_sd,
                     const uint32_t& start_t,
                     const uint32_t& max_t,
                     const double& min_N,
                     const uint32_t& save_every,
                     pcg64& eng,
                     Progress& prog_bar) {


    info = OneRepInfo(N0, V0, max_t, save_every, perturb_sd);

    uint32_t t = 0;
    bool all_gone = false;

    uint32_t iters = 0;

    while (!all_gone && t < start_t) {

        // Update abundances and traits:
        all_gone = info.iterate(f, a0, C, r0, D, add_var, min_N);
        t++;
        prog_bar.increment();

        // Check for user interrupt:
        if (interrupt_check(iters, prog_bar, 1000)) return -1;

    }


    // perturb trait values
    if (perturb_sd > 0) info.perturb(eng);

    t = 0;
    // Save starting info:
    if (save_every > 0) info.save_time(t);

    uint32_t n_incr = 0;
    while (!all_gone && t < max_t) {

        n_incr++;

        // Update abundances and traits:
        all_gone = info.iterate(f, a0, C, r0, D, add_var, min_N);

        if (save_every > 0 && (t % save_every == 0 || (t+1) == max_t || all_gone)) {
            info.save_time(t + 1);
            prog_bar.increment(n_incr);
            n_incr = 0;
        }

        t++;

        // Check for user interrupt:
        if (interrupt_check(iters, prog_bar, 1000)) return -1;

    }

    // Calculate final fitnesses and selection pressure to see if we're at equilibrium
    info.fitness_selection(f, a0, C, r0, D);


    return 0;
}


//' Multiple repetitions of quantitative genetics.
//'
//'
//' @noRd
//'
//[[Rcpp::export]]
List quant_gen_cpp(const uint32_t& n_reps,
                  const std::vector<arma::rowvec>& V0,
                  const std::vector<double>& N0,
                  const double& f,
                  const double& a0,
                  const arma::mat& C,
                  const double& r0,
                  const arma::mat& D,
                  const arma::vec& add_var,
                  const double& perturb_sd,
                  const uint32_t& start_t,
                  const uint32_t& max_t,
                  const double& min_N,
                  const uint32_t& save_every,
                  const bool& show_progress,
                  const uint32_t& n_threads) {

    if (N0.size() != V0.size()) stop("N0.size() != V0.size()");
    if (add_var.n_elem != V0.size()) stop("add_var.n_elem != V0.size()");


    if (C.n_cols != V0[0].n_elem) stop("C.n_cols != q");
    if (C.n_rows != V0[0].n_elem) stop("C.n_rows != q");
    if (D.n_cols != V0[0].n_elem) stop("D.n_cols != q");
    if (D.n_rows != V0[0].n_elem) stop("D.n_rows != q");

    std::vector<OneRepInfo> rep_infos(n_reps);

    const std::vector<std::vector<uint128_t>> seeds = mc_seeds_rep(n_reps);

    Progress prog_bar(n_reps * (max_t + start_t), show_progress);
    bool interrupted = false;

    #ifdef _OPENMP
    #pragma omp parallel default(shared) num_threads(n_threads) if (n_threads > 1)
    {
    #endif


    #ifdef _OPENMP
    uint32_t active_thread = omp_get_thread_num();
    #else
    uint32_t active_thread = 0;
    #endif

    pcg64 eng;

    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (uint32_t i = 0; i < n_reps; i++) {
        eng.seed(seeds[i][0], seeds[i][1]);
        int status = one_quant_gen__(rep_infos[i], V0, N0, f, a0, C, r0, D, add_var,
                                     perturb_sd, start_t, max_t, min_N, save_every,
                                     eng, prog_bar);
        if (active_thread == 0 && status != 0) interrupted = true;
    }
    #ifdef _OPENMP
    }
    #endif

    if (interrupted) {
        throw(Rcpp::exception("\nUser interrupted process.", false));
    }

    /*
     ------------
     Now organize output:
     ------------
     */
    arma::mat nv; // N and V - either through time or just final values
    arma::mat fs(n_reps, 2); // fitness and selection
    uint32_t q = V0[0].n_cols;

    /*
     Go through one time to calculate the # surviving species for all reps and
     for all time point(s) saved.
     Then go back through and fill in values.
     */
    uint32_t total_n_spp = 0;
    if (save_every > 0) {  //   ---- Saving values through time: ----

        for (uint32_t i = 0; i < n_reps; i++) {
            for (uint32_t j = 0; j < rep_infos[i].t.size(); j++) {
                total_n_spp += rep_infos[i].N_t[j].size();
            }
        }
        nv.set_size(total_n_spp, 4 + q);

        uint32_t j = 0;
        for (uint32_t i = 0; i < n_reps; i++) {

            Rcpp::checkUserInterrupt();
            const OneRepInfo& info(rep_infos[i]);
            for (uint32_t t = 0; t < info.t.size(); t++) {

                const std::vector<double>& N_t(info.N_t[t]);
                const std::vector<arma::rowvec>& V_t(info.V_t[t]);
                const std::vector<uint32_t>& spp_t(info.spp_t[t]);
                const double& t_(info.t[t]);
                for (uint32_t k = 0; k < N_t.size(); k++) {

                    nv(j+k,0) = i + 1;      // rep
                    nv(j+k,1) = t_;         // time
                    nv(j+k,2) = spp_t[k];   // species
                    nv(j+k,3) = N_t[k];     // N
                    // V:
                    for (uint32_t l = 0; l < q; l++) nv(j+k, 4+l) = V_t[k](l);

                }

                j += N_t.size();

            }

            fs(i,0) = info.fitness;
            fs(i,1) = info.selection;

        }

    } else {  //                ----  Just saving final values:  ----

        for (uint32_t i = 0; i < n_reps; i++) {
            total_n_spp += rep_infos[i].N.size();
        }
        nv.set_size(total_n_spp, 3 + q);
        uint32_t j = 0;
        for (uint32_t i = 0; i < n_reps; i++) {
            Rcpp::checkUserInterrupt();
            const OneRepInfo& info(rep_infos[i]);
            for (uint32_t k = 0; k < info.N.size(); k++) {
                nv(j+k,0) = i + 1;      // rep
                nv(j+k,1) = k + 1;      // species
                nv(j+k,2) = info.N[k];  // N
                // V:
                for (uint32_t l = 0; l < q; l++) nv(j+k, 3+l) = info.V[k](l);
            }
            j += info.N.size();
            fs(i,0) = info.fitness;
            fs(i,1) = info.selection;
        }

    }

    List out = List::create(_["NV"] = nv, _["FS"] = fs);

    return out;

}

