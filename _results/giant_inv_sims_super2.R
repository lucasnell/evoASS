
#'
#' This file does a "giant"  number of simulations across multiple parameters,
#' looking for how well species invade a 1-species community.
#'
#' `_super2` is for this being super-additive simulations with the second
#' trait state (V1 being high while V2 is low) being the resident
#' species starting trait values.
#'
#' Only d1 is changed for all these simulations
#'

suppressPackageStartupMessages({
    library(sauron)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(magrittr)
    library(parallel)
})
options(dplyr.summarise.inform = FALSE)


.N_THREADS <- 12

#'
#' These simulations are split into 7 sets.
#' Which set this one should conduct depends on the input to this script.
#'
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[[1]]) + 1L


seeds <- c(1721377756, 475139442, 683847556, 846947188,
           282136409, 215494285, 1338715051)


one_sim_combo <- function(.d1, .sigma_V, .sigma_N) {

    .eta = 0.6

    # .sigma_V = 0.05; .sigma_N = 0.0; .eta = -0.6; .d1 = -0.1

    # These simulations are totally deterministic, so no need to run >1 time
    if (.sigma_V == 0 && .sigma_N == 0) {
        .nreps <- 1
        .N_THREADS <- 1
    } else .nreps <- 96

    one_V12_combo <- function(.v1, .v2) {
        ..seed <- sample.int(2^31-1, 1)
        set.seed(..seed)
        X <- quant_gen(eta = .eta, d = c(.d1, 0.1), q = 2, n = 2,
                       V0 = cbind(t(stable_points(.eta)[2,]),
                                  c(.v1, .v2)),
                       sigma_V0 = 0,
                       N0 = c(1, 1),
                       spp_gap_t = 1e3L, final_t = 20e3L, save_every = 0L,
                       sigma_V = .sigma_V,
                       sigma_N = .sigma_N,
                       add_var = rep(0.05, 2),
                       n_reps = .nreps, n_threads = .N_THREADS,
                       show_progress = FALSE)
        X$seed = ..seed
        return(X)
    }

    Z <- crossing(V1 = seq(0, 4, 0.2),
                  V2 = seq(0, 4, 0.2)) %>%
        # bc `seq` makes weird numbers that are ~1e-15 from what they should be:
        mutate(across(.fns = round, digits = 1)) %>%
        mutate(sigma_V = .sigma_V, sigma_N = .sigma_N,
               eta = .eta, d1 = .d1)
    Z$qg = pmap(list(.v1 = Z$V1, .v2 = Z$V2), one_V12_combo)

    return(Z)

}




giant_sims <- crossing(.d1 = c(-0.1, -0.05, 0),
                       .sigma_V = seq(0, 0.3, 0.05),
                       .sigma_N = seq(0, 0.3, 0.05)) %>%
    # bc `seq` makes weird numbers that are ~1e-15 from what they should be:
    mutate(across(.fns = round, digits = 3))

n <- nrow(giant_sims)
k <- length(seeds)

i_rows <- ((i-1) * (n / k) + 1):(i * (n / k))

giant_sims <- giant_sims[i_rows,]

set.seed(seeds[i])
giant_sims <- giant_sims %>%
    pmap_dfr(one_sim_combo)


saveRDS(giant_sims, sprintf("giant_inv_sims_super2_%i.rds", i - 1L))



extract_qg <- function(j) {
    giant_sims[["qg"]][[j]]$nv %>%
        filter(trait == 1) %>%
        group_by(rep) %>%
        summarize(inv = any(spp == 2, na.rm = TRUE),
                  res = any(spp == 1, na.rm = TRUE),
                  coexist = inv & res,
                  replace = inv & !res,
                  reject = !inv & res,
                  extinct = any(is.na(spp)),
                  total = 1) %>%
        ungroup() %>%
        select(-rep, -inv, -res) %>%
        summarize(across(.fns = sum))
}

giant_sims[["qg"]] <- mclapply(1:nrow(giant_sims), extract_qg,
                               mc.cores = .N_THREADS)

giant_sims <- unnest(giant_sims, qg)

saveRDS(giant_sims, sprintf("giant_inv_sims_super2_outcomes_%i.rds", i - 1L))
