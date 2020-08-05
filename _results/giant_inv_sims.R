
#'
#' This file does a "giant"  number of simulations across multiple parameters,
#' looking for how well species invade a 1-species community.
#'
#'
#'

suppressPackageStartupMessages({
    library(sauron)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(magrittr)
})
options(dplyr.summarise.inform = FALSE)


.N_THREADS <- 12

#'
#' These simulations are split into 9 sets.
#' Which set this one should conduct depends on the input to this script.
#'
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[[1]]) + 1L


seeds <- c(668394673,  53842794,  1053078634, 777531456,
           2026768533, 152103205, 1786315847, 182943636)


one_sim_combo <- function(.sigma_V, .sigma_N, .eta, .d1) {

    # .sigma_V = 0.05; .sigma_N = 0.0; .eta = -0.6; .d1 = -0.1

    Z <- crossing(.v1 = seq(0, 3, 0.1),
                  .v2 = seq(0, 3, 0.1)) %>%
        pmap_dfr(function(.v1, .v2) {
            ..seed <- sample.int(2^31-1, 1)
            set.seed(..seed)
            X <- quant_gen(eta = .eta, d = c(.d1, 0.1), q = 2, n = 2,
                           V0 = cbind(t(stable_points(.eta)[1,]),
                                      c(.v1, .v2)),
                           sigma_V0 = 0,
                           N0 = c(1, 1),
                           spp_gap_t = 1e3L, final_t = 20e3L, save_every = 0L,
                           sigma_V = .sigma_V,
                           sigma_N = .sigma_N,
                           add_var = rep(0.05, 2),
                           n_reps = 96, n_threads = .N_THREADS,
                           show_progress = FALSE)
            X$nv %>%
                filter(trait == 1) %>%
                group_by(rep) %>%
                summarize(inv = any(spp == 2),
                          res = any(spp == 1)) %>%
                select(-rep) %>%
                summarise(across(.fns = sum)) %>%
                mutate(V1 = .v1, V2 = .v2, seed = ..seed)
        }) %>%
        mutate(sigma_V = .sigma_V, sigma_N = .sigma_N,
               eta = .eta, d1 = .d1)

}



giant_sims <- crossing(.sigma_V = seq(0, 0.2, 0.05),
                          .sigma_N = seq(0, 0.2, 0.05),
                          .eta = c(-0.6, 0.6),
                          .d1 = seq(-0.1, 0.05, 0.01)) %>%
    # bc `seq` makes weird numbers that are ~1e-15 from what they should be:
    mutate(across(where(is.numeric), .fns = round, digits = 2))

n <- nrow(giant_sims)
k <- length(seeds)

i_rows <- ((i-1) * (n / k) + 1):(i * (n / k))

giant_sims <- giant_sims[i_rows,]

set.seed(seeds[i])
giant_sims <- giant_sims %>%
    pmap_dfr(one_sim_combo)


saveRDS(giant_sims, sprintf("giant_inv_sims_%i.rds", i - 1L))

