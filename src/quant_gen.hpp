#ifndef __SAURON_QUANT_GEN_H
#define __SAURON_QUANT_GEN_H


#include <RcppArmadillo.h>
#include "sim.hpp"

using namespace Rcpp;


void sel_str__(arma::mat& ss_mat,
               const std::vector<arma::rowvec>& V,
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
    std::vector<arma::rowvec> V;    // traits
    std::vector<uint32_t> spp;      // species indexes (based on N0 and V0)
    double fitness;
    double selection;
    // Info for output if tracking through time:
    std::vector<double> t;
    std::vector<std::vector<double>> N_t;
    std::vector<std::vector<arma::rowvec>> V_t;
    std::vector<std::vector<uint32_t>> spp_t;

    OneRepInfo () {};
    OneRepInfo(const std::vector<double>& N_,
               const std::vector<arma::rowvec>& V_,
               const uint32_t& max_t,
               const uint32_t& save_every,
               const double& perturb_sd)
        : N(N_), V(V_), spp(N_.size()), fitness(-1), selection(-1),
          t(), N_t(), V_t(),
          A(V_.size()),
          ss_mat(V_.size(), V_[0].n_elem),
          n(N_.size()), q(V_[0].n_elem),
          norm_distr(0.0, perturb_sd) {

        for (uint32_t i = 0; i < N_.size(); i++) spp[i] = i;

        if (save_every > 0) {

            uint32_t n_saves = static_cast<uint32_t>(
                std::ceil(static_cast<double>(max_t - 1) /
                    static_cast<double>(save_every))) + 1U;
            t.reserve(n_saves);
            N_t.reserve(n_saves);
            V_t.reserve(n_saves);
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
                 const double& min_N) {

        /*
         Update abundances
         */
        // Setting up vector of extinct clones (if any):
        std::vector<uint32_t> extinct;
        extinct.reserve(V.size());
        // Fill in density dependences:
        A_VN_<std::vector<double>>(A, V, N, a0, D);
        // Fill in abundances:
        for (uint32_t i = 0; i < A.size(); i++) {
            double r = r_V_(V[i], f, C, r0);
            N[i] *= std::exp(r - A[i]);
            // See if it goes extinct:
            if (N[i] < min_N) extinct.push_back(i);
        }

        // If everything is gone, clear vectors and stop simulations:
        if (extinct.size() == N.size()) {
            N.clear();
            V.clear();
            A.clear();
            spp.clear();
            return true;
        }

        /*
         Update traits
         */
        // Fill in selection-strength matrix:
        sel_str__(ss_mat, V, N, f, a0, C, r0, D);
        // Then include additive genetic variance when adding to trait values:
        for (uint32_t i = 0; i < V.size(); i++) {
            V[i] += (add_var(i) * ss_mat.row(i));
            for (double& v : V[i]) if (v < 0) v *= -1; // <-- keeping traits >= 0
        }

        /*
         Remove extinct clones (starting at the back):
         */
        for (uint32_t i = 0, j; i < extinct.size(); i++) {
            j = extinct.size() - i - 1;
            N.erase(N.begin() + extinct[j]);
            V.erase(V.begin() + extinct[j]);
            A.erase(A.begin() + extinct[j]);
            spp.erase(spp.begin() + extinct[j]);
        }

        return false;
    }

    // perturb trait values
    void perturb(pcg64& eng) {
        for (uint32_t i = 0; i < V.size(); i++) {
            for (double& v : V[i]) {
                v += norm_distr(eng);
                if (v < 0) v *= -1; // <-- keeping traits >= 0
            }
        }
        return;
    }

    // save info for output
    void save_time(const uint32_t& t) {
        this->t.push_back(t);
        // If everything's extinct...
        if (N.size() == 0) {
            // Fill last set of N's with a zero:
            this->N_t.push_back(std::vector<double>(1, 0.0));
            // Fill last V with a `NaN` (closest to NA I know of):
            std::vector<arma::rowvec> V__(1, arma::rowvec(q));
            V__[0].fill(arma::datum::nan);
            this->V_t.push_back(V__);
            // Fill last set of spp's with a zero:
            this->spp_t.push_back(std::vector<uint32_t>(1, 0U));
        } else {
            this->N_t.push_back(N);
            this->V_t.push_back(V);
            this->spp_t.push_back(spp);
        }
        return;
    }


    void fitness_selection(const double& f,
                           const double& a0,
                           const arma::mat& C,
                           const double& r0,
                           const arma::mat& D) {

        // Temporary objects:
        arma::vec WN(this->V.size());
        arma::mat SV;
        // Filling in fitnesses and selection strengths:
        F_t__<arma::vec>(WN, this->V, this->N, f, a0, C, r0, D);
        sel_str__(SV, this->V, this->N, f, a0, C, r0, D);

        // Fill final values:
        this->fitness = arma::prod(WN);
        this->selection = arma::accu(SV);

        return;
    }


private:

    std::vector<double> A;  // Density dependence
    arma::mat ss_mat;       // Selection strength
    uint32_t n;             // Starting # species
    uint32_t q;             // # traits
    std::normal_distribution<double> norm_distr;


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
