#ifndef __SAURON_QUANT_GEN_H
#define __SAURON_QUANT_GEN_H


#include <RcppArmadillo.h>
#include <random>
#include "sim.hpp"

using namespace Rcpp;


typedef std::normal_distribution<double> normal_distr;


void sel_str__(arma::mat& ss_mat,
               const std::vector<arma::vec>& V,
               const std::vector<double>& N,
               const double& f,
               const double& a0,
               const arma::mat& C,
               const double& r0,
               const arma::mat& D);

/*
 Output info for one repetition:
 */
class OneRepInfo {
public:

    std::vector<double> N;      // abundances
    std::vector<arma::vec> V;   // traits - genotypes
    std::vector<arma::vec> Vp;  // traits - phenotypes
    std::vector<uint32_t> spp;  // species indexes (based on N0 and V0)
    // Info for output if tracking through time:
    std::vector<double> t;
    std::vector<std::vector<double>> N_t;
    std::vector<std::vector<arma::vec>> V_t;
    std::vector<std::vector<arma::vec>> Vp_t;
    std::vector<std::vector<uint32_t>> spp_t;

    OneRepInfo () {};
    OneRepInfo(const std::vector<double>& N_,
               const std::vector<arma::vec>& V_,
               const uint32_t& max_t,
               const uint32_t& save_every,
               const double& perturb_sd,
               const double& sigma_V,
               pcg64& eng)
        : N(N_), V(V_), Vp(V_), spp(N_.size()),
          t(), N_t(), V_t(),
          A(V_.size()),
          ss_mat(),
          n(N_.size()),
          q(V_[0].n_elem),
          perturb_sd_(perturb_sd) {

        for (uint32_t i = 0; i < N_.size(); i++) spp[i] = i + 1;

        if (save_every > 0) {

            uint32_t n_saves = static_cast<uint32_t>(
                std::ceil(static_cast<double>(max_t - 1) /
                    static_cast<double>(save_every))) + 2U;
            t.reserve(n_saves);
            N_t.reserve(n_saves);
            V_t.reserve(n_saves);
            Vp_t.reserve(n_saves);
            spp_t.reserve(n_saves);

        }

        // adding stochasticity to starting phenotypes
        if (sigma_V > 0) {

            for (uint32_t i = 0; i < Vp.size(); i++) {
                for (uint32_t j = 0; j < Vp[i].n_elem; j++) {
                    Vp[i][j] *= std::exp(rnorm(eng) * sigma_V);
                }
            }

        }

        return;

    };
    OneRepInfo(const std::vector<double>& N_,
               const std::vector<arma::vec>& V_,
               const std::vector<arma::vec>& Vp_,
               const uint32_t& max_t,
               const uint32_t& save_every,
               const double& perturb_sd)
        : N(N_), V(V_), Vp(Vp_), spp(N_.size()),
          t(), N_t(), V_t(),
          A(V_.size()),
          ss_mat(),
          n(N_.size()),
          q(V_[0].n_elem),
          perturb_sd_(perturb_sd) {

        for (uint32_t i = 0; i < N_.size(); i++) spp[i] = i + 1;

        if (save_every > 0) {

            uint32_t n_saves = static_cast<uint32_t>(
                std::ceil(static_cast<double>(max_t - 1) /
                    static_cast<double>(save_every))) + 2U;
            t.reserve(n_saves);
            N_t.reserve(n_saves);
            V_t.reserve(n_saves);
            Vp_t.reserve(n_saves);
            spp_t.reserve(n_saves);

        }

        return;

    };


    /*
     One iteration that updates abundances and traits.
     It returns a boolean for whether all species are extinct.
     */
    bool iterate(const double& f,
                 const double& a0,
                 const arma::mat& C,
                 const double& r0,
                 const arma::mat& D,
                 const arma::vec& add_var,
                 const double& min_N,
                 const double& sigma_N,
                 const double& sigma_V,
                 pcg64& eng) {

        /*
         Update abundances
         */
        // Setting up vector of extinct clones (if any):
        std::vector<uint32_t> extinct;
        extinct.reserve(V.size());
        // Fill in density dependences:
        A_VN_<std::vector<double>>(A, Vp, N, a0, D);
        // Fill in abundances:
        for (uint32_t i = 0; i < A.size(); i++) {
            double r = r_V_<arma::vec>(Vp[i], f, C, r0);
            if (sigma_N <= 0) {
                N[i] *= std::exp(r - A[i]);
            } else N[i] *= std::exp(r - A[i] + rnorm(eng) * sigma_N);
            // See if it goes extinct:
            if (N[i] < min_N) extinct.push_back(i);
        }

        // If everything is gone, clear vectors and stop simulations:
        if (extinct.size() == N.size()) {
            N.clear();
            V.clear();
            Vp.clear();
            A.clear();
            spp.clear();
            return true;
        }

        /*
         Update traits
         */
        // Fill in selection-strength matrix:
        sel_str__(ss_mat, Vp, N, f, a0, C, r0, D);
        // Then include additive genetic variance when adding to trait values:
        if (sigma_V > 0) {
            for (uint32_t i = 0; i < V.size(); i++) {
                for (uint32_t j = 0; j < Vp[i].n_elem; j++) {
                    V[i][j] += (add_var(i) * ss_mat(j,i));
                    if (V[i][j] < 0) V[i][j] = 0; // <-- keeping traits >= 0
                    // including stochasticity:
                    Vp[i][j] = V[i][j] * std::exp(rnorm(eng) * sigma_V);
                }
            }
        } else {
            for (uint32_t i = 0; i < V.size(); i++) {
                for (uint32_t j = 0; j < Vp[i].n_elem; j++) {
                    V[i][j] += (add_var(i) * ss_mat(j,i));
                    if (V[i][j] < 0) V[i][j] = 0; // <-- keeping traits >= 0
                    Vp[i][j] = V[i][j];
                }
            }
        }

        /*
         Remove extinct clones (starting at the back):
         */
        for (uint32_t i = 0, j; i < extinct.size(); i++) {
            j = extinct.size() - i - 1;
            N.erase(N.begin() + extinct[j]);
            V.erase(V.begin() + extinct[j]);
            Vp.erase(Vp.begin() + extinct[j]);
            A.erase(A.begin() + extinct[j]);
            spp.erase(spp.begin() + extinct[j]);
        }

        return false;
    }

    // perturb trait values
    void perturb(const double& sigma_V, pcg64& eng) {
        if (sigma_V > 0) {
            for (uint32_t i = 0; i < V.size(); i++) {
                for (uint32_t j = 0; j < V[i].n_elem; j++) {
                    V[i][j] = trunc_rnorm_(V[i][j], perturb_sd_, eng);
                    // including stochasticity:
                    Vp[i][j] = V[i][j] * std::exp(rnorm(eng) * sigma_V);
                }
            }
        } else {
            for (uint32_t i = 0; i < V.size(); i++) {
                for (uint32_t j = 0; j < V[i].n_elem; j++) {
                    V[i][j] = trunc_rnorm_(V[i][j], perturb_sd_, eng);
                    Vp[i][j] = V[i][j];
                }
            }
        }
        return;
    }

    // save info for output
    void save_time(const uint32_t& t_) {
        t.push_back(t_);
        // If everything's extinct...
        if (N.size() == 0) {
            // Fill last set of N's with a zero:
            N_t.push_back(std::vector<double>(1, 0.0));
            // Fill last V and Vp with a `NaN` (closest to NA I know of):
            std::vector<arma::vec> V__(1, arma::vec(q));
            V__[0].fill(arma::datum::nan);
            V_t.push_back(V__);
            Vp_t.push_back(V__);
            // Fill last set of spp's with a zero:
            spp_t.push_back(std::vector<uint32_t>(1, 0U));
        } else {
            N_t.push_back(N);
            V_t.push_back(V);
            Vp_t.push_back(Vp);
            spp_t.push_back(spp);
        }
        return;
    }


private:

    std::vector<double> A;  // Density dependence
    arma::mat ss_mat;       // Selection strength
    uint32_t n;             // Starting # species
    uint32_t q;             // # traits
    double perturb_sd_;
    normal_distr rnorm = normal_distr(0, 1);


};




// Below is useful when doing plots of fitness landscapes:
// /*
//  Data for each combination of trial, time, and species.
//
//  Vectors for trait values are used to reconstruct vector of combinations of trait
//  values.
//  These vectors should all be the same length.
//  */
// class OneComboData {
// public:
//
//     uint32_t trial;
//     uint32_t time;
//     uint32_t sp;
//     std::vector<double> trait_mins;
//     std::vector<double> trait_maxs;
//     std::vector<double> trait_incr;
//     std::vector<double> fitness;
//
// };





#endif
