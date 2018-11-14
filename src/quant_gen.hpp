#ifndef __EVOASS_QUANT_GEN_H
#define __EVOASS_QUANT_GEN_H


#include <RcppArmadillo.h>
#include "sim.hpp"

using namespace Rcpp;


/*
 Output info for one repetition:
 */
class OneRepInfo {
public:

    std::vector<double> N;
    std::vector<arma::rowvec> V;
    double fitness;
    double selection;
    std::vector<double> t;
    std::vector<std::vector<double>> N_t;
    std::vector<std::vector<arma::rowvec>> V_t;

    OneRepInfo () {};
    OneRepInfo(const std::vector<double>& N_,
               const std::vector<arma::rowvec>& V_)
        : N(N_), V(V_), fitness(-1), selection(-1),
          t(), N_t(), V_t() {};
    // Don't use this way if save_every == 0!
    OneRepInfo(const std::vector<double>& N_,
               const std::vector<arma::rowvec>& V_,
               const uint32_t& max_t,
               const uint32_t& save_every)
        : N(N_), V(V_), fitness(-1), selection(-1),
          t(((max_t - 1) / save_every) + 2),
          N_t(((max_t - 1) / save_every) + 2, std::vector<double>(N.size())),
          V_t(((max_t - 1) / save_every) + 2,
              std::vector<arma::rowvec>(V_.size(), arma::rowvec(V_[0].n_elem))) {};

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
