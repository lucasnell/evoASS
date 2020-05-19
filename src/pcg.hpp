#ifndef __SAURON_PCG_H
#define __SAURON_PCG_H

#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <pcg/pcg_extras.hpp>  // pcg 128-bit integer type
#include <pcg/pcg_random.hpp>  // pcg prng



using namespace Rcpp;


typedef pcg_extras::pcg128_t uint128_t;



namespace pcg {
    const long double max = static_cast<long double>(pcg64::max());
}



/*
 ========================

 Seeding

 ========================
 */



// To sample for seeds before multi-core operations - one set of seeds per thread
inline std::vector<std::vector<uint64_t>> mc_seeds(const uint32_t& n_threads) {

    std::vector<std::vector<uint64_t>> sub_seeds(n_threads, std::vector<uint64_t>(8));

    for (uint32_t i = 0; i < n_threads; i++) {
        sub_seeds[i] = as<std::vector<uint64_t>>(Rcpp::runif(8,0,4294967296));
    }

    return sub_seeds;
}

// Fill two 128-bit seeds from 8 32-bit seeds (casted to 64-bit)
inline void fill_seeds(const std::vector<uint64_t>& sub_seeds,
                          uint128_t& seed1, uint128_t& seed2) {

    uint128_t seed64_1 = (sub_seeds[0]<<32) + sub_seeds[1];
    uint128_t seed64_2 = (sub_seeds[2]<<32) + sub_seeds[3];
    uint128_t seed64_3 = (sub_seeds[4]<<32) + sub_seeds[5];
    uint128_t seed64_4 = (sub_seeds[6]<<32) + sub_seeds[7];

    seed1 = (seed64_1<<64) + seed64_2;
    seed2 = (seed64_3<<64) + seed64_4;

    return;
}



// To sample for seeds before multi-core operations - one set of seeds per rep
inline std::vector<std::vector<uint128_t>> mc_seeds_rep(const uint32_t& n_reps) {

    std::vector<std::vector<uint128_t>> sub_seeds(n_reps, std::vector<uint128_t>(2));

    std::vector<uint64_t> tmp(8);

    for (uint32_t i = 0; i < n_reps; i++) {

        for (uint32_t j = 0; j < 8; j++) {
            tmp[j] = static_cast<uint64_t>(R::runif(0,4294967296));
        }
        fill_seeds(tmp, sub_seeds[i][0], sub_seeds[i][1]);
    }

    return sub_seeds;
}




/*
 For single-core operations, you can use R's RNG for 32-bit random number generation.
 */
inline pcg64 seeded_pcg() {

    // 32-bit seeds from unif_rand
    std::vector<uint64_t> sub_seeds(8);
    for (uint64_t& s : sub_seeds) s = R::runif(0, 4294967296);
    uint128_t seed1;
    uint128_t seed2;

    fill_seeds(sub_seeds, seed1, seed2);

    pcg64 out(seed1, seed2);
    return out;
}

/*
 For multi-core operations, you should call `mc_seeds` when outside multi-core mode,
 then input an inner `std::vector<uint32_t>` (from inside the object output from `mc_seeds`)
 to this function when in multi-core mode to seed the PRNG.
 */
// sub_seeds needs to be at least 8-long!
inline pcg64 seeded_pcg(const std::vector<uint64_t>& sub_seeds) {

    uint128_t seed1;
    uint128_t seed2;
    fill_seeds(sub_seeds, seed1, seed2);

    pcg64 out(seed1, seed2);

    return out;
}


/*
 ========================

 Number generation

 ========================
 */

// uniform in range [0,1]
inline long double runif_0011(pcg64& eng) {
    return static_cast<long double>(eng()) / pcg::max;
}
// uniform in range [0,1)
inline long double runif_001(pcg64& eng) {
    return static_cast<long double>(eng()) / (pcg::max + 1);
}
// uniform in range (0,1)
inline long double runif_01(pcg64& eng) {
    return (static_cast<long double>(eng()) + 1) / (pcg::max + 2);
}
// uniform in range (a,b)
inline long double runif_ab(pcg64& eng, const long double& a, const long double& b) {
    return a + ((static_cast<long double>(eng()) + 1) / (pcg::max + 2)) * (b - a);
}
// uniform in range [a,b]
inline uint64_t runif_aabb(pcg64& eng, const uint64_t& a, const uint64_t& b) {
    return a + (static_cast<long double>(eng()) / (pcg::max + 1)) * (b - a + 1);
}







#endif
