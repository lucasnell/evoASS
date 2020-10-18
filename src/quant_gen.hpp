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

    std::vector<double> N;          // abundances
    std::vector<arma::vec> V;       // traits - genotypes
    std::vector<arma::vec> Vp;      // traits - phenotypes
    std::vector<double> add_var;    // additive genetic variances
    std::vector<uint32_t> spp;      // species indexes (based on N0 and V0)
    uint32_t n;                     // Total # species added
    // Info for output if tracking through time:
    std::vector<double> t;
    std::vector<std::vector<double>> N_t;
    std::vector<std::vector<arma::vec>> V_t;
    std::vector<std::vector<arma::vec>> Vp_t;
    std::vector<std::vector<uint32_t>> spp_t;

    OneRepInfo () {};
    OneRepInfo(const std::deque<double>& N_,
               const std::deque<arma::vec>& V_,
               const std::deque<arma::vec>& Vp_,
               const std::deque<double>& add_var_)
        : N(N_.begin(), N_.end()),
          V(V_.begin(), V_.end()),
          Vp(Vp_.begin(), Vp_.end()),
          add_var(add_var_.begin(), add_var_.end()),
          spp(N_.size()),
          n(N_.size()),
          t(), N_t(), V_t(),
          A(N_.size()),
          ss_mat(),
          q(V_[0].n_elem) {

        if (V_.size() != N_.size()) stop("\nV_.size() != N_.size()");
        if (Vp_.size() != V_.size()) stop("\nVp_.size() != V_.size()");
        for (uint32_t i = 0; i < N_.size(); i++) {
            if (Vp_[i].n_elem != V_[i].n_elem) {
                stop(std::string("\nVp and V sizes don't match at index ") +
                    std::to_string(i));
            }
            spp[i] = i + 1;
        }

        return;

    };
    OneRepInfo(const double& N_,
               const arma::vec& V_,
               const arma::vec& Vp_,
               const double& add_var_)
        : N(1, N_),
          V(1, V_),
          Vp(1, Vp_),
          add_var(1, add_var_),
          spp(1, 1),
          n(1),
          t(), N_t(), V_t(),
          A(1),
          ss_mat(),
          q(V_.n_elem) {

        if (Vp_.n_elem != V_.n_elem) {
            stop("\nVp and V sizes don't match for first species");
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
                 const double& min_N,
                 const double& sigma_N,
                 const std::vector<double>& sigma_V,
                 const std::vector<double>& mu_V,
                 const bool& lnorm_V,
                 pcg64& eng) {

        uint32_t current_n = V.size(); // current # species (`n` is total added)

        /*
         This is for iterations where species will be later added, but
         all species that were already added have gone extinct:
         */
        if (current_n == 0) return true;

        /*
         Update abundances
         */
        // Setting up vector of extinct clones (if any):
        std::vector<uint32_t> extinct;
        extinct.reserve(current_n);
        // Fill in density dependences:
        A_VN_<std::vector<double>>(A, Vp, N, a0, D);
        // Fill in abundances:
        for (uint32_t i = 0; i < current_n; i++) {
            double r = r_V_<arma::vec>(Vp[i], f, C, r0);
            if (sigma_N <= 0) {
                N[i] *= std::exp(r - A[i]);
            } else N[i] *= std::exp(r - A[i] + rand_norm(eng) * sigma_N);
            // See if it goes extinct:
            if (N[i] < min_N) extinct.push_back(i);
        }

        // If everything is gone, clear vectors and stop simulations:
        if (extinct.size() == current_n) {
            rm_all();
            return true;
        }

        /*
         Update traits
         */
        // Fill in selection-strength matrix:
        sel_str__(ss_mat, Vp, N, f, a0, C, r0, D);

        /*
         Then include additive genetic variance when adding to trait values.
         Also add stochasticity to phenotypes if necessary.
         */
        change_V(sigma_V, mu_V, lnorm_V, eng);

        /*
         Remove extinct clones (starting at the back):
         */
        for (uint32_t i = 0, j; i < extinct.size(); i++) {
            j = extinct.size() - i - 1;
            rm_species(extinct[j]);
        }

        return false;
    }

    // add a new species to community
    void add_species(const double& new_N,
                     const arma::vec& new_V,
                     const arma::vec& new_Vp,
                     const double& new_add_var) {

        N.push_back(new_N);
        V.push_back(new_V);
        Vp.push_back(new_Vp);
        add_var.push_back(new_add_var);

        n++;
        spp.push_back(n);
        A.push_back(0);

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


    void reserve(const uint32_t& n_saves) {
        t.reserve(n_saves);
        N_t.reserve(n_saves);
        V_t.reserve(n_saves);
        Vp_t.reserve(n_saves);
        spp_t.reserve(n_saves);
        return;
    }


private:

    std::vector<double> A;  // Density dependence
    arma::mat ss_mat;       // Selection strength
    uint32_t q;             // # traits
    normal_distr rand_norm = normal_distr(0, 1);


    inline void change_V(const std::vector<double>& sigma_V,
                         const std::vector<double>& mu_V,
                         const bool& lnorm_V,
                         pcg64& eng) {
        for (uint32_t j = 0; j < q; j++) {
            if (sigma_V[j] > 0) {
                if (lnorm_V) {
                    change_V_lnorm(sigma_V, mu_V, j, eng);
                } else {
                    change_V_norm(sigma_V, j, eng);
                }
            } else {
                change_V_determ(j);
            }
        }
        return;
    }

    inline void change_V_lnorm(const std::vector<double>& sigma_V,
                               const std::vector<double>& mu_V,
                               const uint32_t& j,
                               pcg64& eng) {
        for (uint32_t i = 0; i < V.size(); i++) {
            V[i][j] += (add_var[i] * ss_mat(j,i));
            if (V[i][j] < 0) V[i][j] = 0; // <-- keeping traits >= 0
            Vp[i][j] = V[i][j];
            // including stochasticity:
            Vp[i][j] *= std::exp(rand_norm(eng) * sigma_V[j] + mu_V[j]);
        }
        return;
    }

    inline void change_V_norm(const std::vector<double>& sigma_V,
                              const uint32_t& j,
                              pcg64& eng) {
        for (uint32_t i = 0; i < V.size(); i++) {
            V[i][j] += (add_var[i] * ss_mat(j,i));
            if (V[i][j] < 0) V[i][j] = 0; // <-- keeping traits >= 0
            Vp[i][j] = trunc_rnorm_(V[i][j], sigma_V[j], eng);
        }
        return;
    }

    inline void change_V_determ(const uint32_t& j) {
        for (uint32_t i = 0; i < V.size(); i++) {
            V[i][j] += (add_var[i] * ss_mat(j,i));
            if (V[i][j] < 0) V[i][j] = 0; // <-- keeping traits >= 0
            Vp[i][j] = V[i][j];
        }
        return;
    }


    void rm_species(const uint32_t& idx) {

        N.erase(N.begin() + idx);
        V.erase(V.begin() + idx);
        Vp.erase(Vp.begin() + idx);
        add_var.erase(add_var.begin() + idx);
        A.erase(A.begin() + idx);
        spp.erase(spp.begin() + idx);

        return;

    }


    void rm_all() {

        N.clear();
        V.clear();
        Vp.clear();
        add_var.clear();
        A.clear();
        spp.clear();

        return;

    }


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
