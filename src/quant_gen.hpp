#ifndef __EASS_QUANT_GEN_H
#define __EASS_QUANT_GEN_H


#include <RcppArmadillo.h>
#include "sim.hpp"

using namespace Rcpp;

/*
 Data for each combination of trial, time, and species.

 Vectors for trait values are used to reconstruct vector of combinations of trait values.
 These vectors should all be the same length.
 */
class OneComboData {
public:

    uint32_t trial;
    uint32_t time;
    uint32_t sp;
    std::vector<double> trait_mins;
    std::vector<double> trait_maxs;
    std::vector<double> trait_incr;
    std::vector<double> fitness;

};



#endif
