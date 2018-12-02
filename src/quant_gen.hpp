#ifndef __EVOASS_QUANT_GEN_H
#define __EVOASS_QUANT_GEN_H


#include <RcppArmadillo.h>
#include "sim.hpp"

using namespace Rcpp;


void sel_str__(arma::mat& ss_mat,
               const std::vector<arma::rowvec>& V,
               const std::vector<double>& N,
               const double& f,
               const double& g,
               const arma::mat& C,
               const double& r0,
               const double& d);

/*
 Output info for one repetition:
 */
class OneRepInfo {
public:

    std::vector<double> N;
    std::vector<arma::rowvec> V;
    double fitness;
    double selection;
    // Info for output if tracking through time:
    std::vector<double> t;
    std::vector<std::vector<double>> N_t;
    std::vector<std::vector<arma::rowvec>> V_t;

    OneRepInfo () {};
    OneRepInfo(const std::vector<double>& N_,
               const std::vector<arma::rowvec>& V_)
        : N(N_), V(V_), fitness(-1), selection(-1),
          t(), N_t(), V_t(), A(V_.size()), ss_mat(V_.size(), V_[0].n_elem) {};
    // Don't use this way if save_every == 0!
    OneRepInfo(const std::vector<double>& N_,
               const std::vector<arma::rowvec>& V_,
               const uint32_t& max_t,
               const uint32_t& save_every)
        : N(N_), V(V_), fitness(-1), selection(-1),
          t(((max_t - 1) / save_every) + 2),
          N_t(((max_t - 1) / save_every) + 2, std::vector<double>(N.size())),
          V_t(((max_t - 1) / save_every) + 2,
              std::vector<arma::rowvec>(V_.size(), arma::rowvec(V_[0].n_elem))),
          A(V_.size()),
          ss_mat(V_.size(), V_[0].n_elem) {};


    /*
     One iteration that updates abundances and traits.
     It returns a boolean for whether all species are extinct.
     */
    bool iterate(const double& f,
                 const double& g,
                 const arma::mat& C,
                 const double& r0,
                 const double& d,
                 const arma::vec& add_var,
                 const double& min_N,
                 const bool& rm_extinct) {

        // Setting up vector of extinct clones (if any):
        std::vector<uint32_t> extinct;
        extinct.reserve(V.size());
        // Fill in density dependences:
        A_VN_<std::vector<double>>(A, V, N, g, d);
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
            return true;
        }

        /*
         Update traits
         */
        // Fill in selection-strength matrix:
        sel_str__(ss_mat, V, N, f, g, C, r0, d);
        // Then include additive genetic variance when adding to trait values:
        for (uint32_t i = 0; i < V.size(); i++) {
            V[i] += (add_var(i) * ss_mat.row(i));
        }

        /*
         Remove extinct clones (starting at the back):
         */
        if (rm_extinct) {
            for (uint32_t i = 0, j; i < extinct.size(); i++) {
                j = extinct.size() - i - 1;
                N.erase(N.begin() + extinct[j]);
                V.erase(V.begin() + extinct[j]);
                A.erase(A.begin() + extinct[j]);
            }
        }

        return false;
    }

    // perturb trait values
    void perturb(pcg64& eng, std::lognormal_distribution<double>& distr) {
        for (uint32_t i = 0; i < V.size(); i++) {
            for (uint32_t j = 0; j < V[i].n_elem; j++) {
                V[i](j) *= distr(eng);
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
            this->N_t.push_back(std::vector<double>(1, 0));
            // Fill last V with a `NaN` (closest to NA I know of):
            std::vector<arma::rowvec> V__(1, arma::rowvec(1));
            V__[0](0) = arma::datum::nan;
            this->V_t.push_back(V__);
        } else {
            this->N_t.push_back(N);
            this->V_t.push_back(V);
        }
        return;
    }


    void fitness_selection(const double& f,
                           const double& g,
                           const arma::mat& C,
                           const double& r0,
                           const double& d) {

        // Temporary objects:
        arma::vec WN(this->V.size());
        arma::mat SV;
        // Filling in fitnesses and selection strengths:
        F_t__<arma::vec>(WN, this->V, this->N, f, g, C, r0, d);
        sel_str__(SV, this->V, this->N, f, g, C, r0, d);

        // Fill final values:
        this->fitness = arma::prod(WN);
        this->selection = arma::accu(SV);

        return;
    }



private:

    std::vector<double> A;  // Density dependence
    arma::mat ss_mat;       // Selection strength


};




// Below is useful when doing plots of fitness landscapes:
// /*
//  Data for each combination of trial, time, and species.
//
//  Vectors for trait values are used to reconstruct vector of combinations of trait values.
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
