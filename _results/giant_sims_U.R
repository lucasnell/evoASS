
#'
#' This file does a "giant"  number of simulations across multiple parameters.
#' The difference here is that species start with traits
#' from a uniform distribution.
#' That way, each invading species isn't as "high quality" as those in
#' `giant_sims.R` simulations.
#'
#' > This doesn't appear to be interesting at all.
#'
#' In all these simulations, both `d` values are changed.
#'
#'

library(sauron)
library(dplyr)
library(tidyr)
library(purrr)
library(magrittr)

.N_THREADS <- 12

#'
#' These simulations are split into 9 sets.
#' Which set this one should conduct depends on the input to this script.
#'
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[[1]]) + 1L


seeds <- c(56965945, 679702122, 1629399850, 457574606, 805131588,
           161027162, 1159880, 918587938, 989618742)




one_sim_combo <- function(.d, .eta, .add_var, .sigma_N, .sigma_V) {

    .n <- 100

    .V0 <- do.call(rbind, lapply(sapply(stable_points(.eta), max) * 2,
                                 runif, n = .n, min = 0))

    .seed <- sample.int(2^31 - 1, 1)
    set.seed(.seed)

    Z <- quant_gen(q = 2, eta = .eta, d = .d, n = .n,
                   add_var = rep(.add_var, .n),
                   V0 = .V0, sigma_V0 = 0,
                   spp_gap_t = 500L, final_t = 50e3L, n_reps = 12,
                   sigma_N = .sigma_N, sigma_V = .sigma_V,
                   save_every = 0L,
                   n_threads = .N_THREADS, show_progress = FALSE)

    if (.sigma_N == 0 && .sigma_V == 0) {
        Z$call[["eta"]] <- eval(.eta)
        Z$call[["d"]] <- eval(.d)
        Z$call[["n"]] <- eval(.n)
        Z$call[["add_var"]] <- eval(rep(.add_var, .n))
        jacs <- jacobians(Z)
    } else jacs <- NULL


    return(list(NV = Z$nv %>%
                    mutate(trait = paste0("V", trait)) %>%
                    spread(trait, geno) %>%
                    mutate(d = .d, eta = .eta, add_var = .add_var,
                           sigma_N = .sigma_N, sigma_V = .sigma_V),
                J = jacs,
                V0 = .V0,
                seed = .seed))
}

giant_sims <- crossing(.d = seq(-0.25, 2, length.out = 10),
                       .add_var = seq(0.01, 0.1, 0.01),
                       .eta = -1:1 * 0.6,
                       .sigma_N = c(0, 0.05, 0.1),
                       .sigma_V = c(0, 0.05, 0.1)) %>%
    # bc `seq` makes weird numbers that are ~1e-15 from what they should be:
    mutate(across(.fns = round, digits = 2))


i_rows <- ((i-1) * (nrow(giant_sims) / 9) + 1):(i * (nrow(giant_sims) / 9))

giant_sims <- giant_sims[i_rows,]

set.seed(seeds[i])
giant_sims <- giant_sims %>%
    pmap(one_sim_combo)


saveRDS(giant_sims, sprintf("giant_sims_U_%i.rds", i - 1L))



#' # ------------------------------------*
#' # __ vary d1 & d2; V0 ~ U ----
#' # ------------------------------------*
#'
#'
#'
#' #'
#' #' NOTE: this function always changes both `d` values
#' #'
#' one_coexist_U_combo <- function(.d, .eta, .add_var, .pb = NULL) {
#'
#'     if (!is.null(.pb)) .pb$tick(0)
#'
#'     # .d = 0.15; .eta = etas[[2]]; .add_var = 0.01
#'     # rm(.d, .eta, .add_var)
#'
#'     stopifnot(is.numeric(.d) && length(.d) == 1)
#'
#'     .d <- rep(.d, 2)
#'     .n <- 100
#'
#'     .V0 <- matrix(runif(.n * 2, 0, 3), nrow = .n)
#'
#'     .seed <- sample.int(2^31 - 1, 1)
#'     set.seed(.seed)
#'
#'     Z <- quant_gen(q = 2, eta = .eta, d = .d, n = .n,
#'                    V0 = .V0, sigma_V0 = 0,
#'                    add_var = rep(.add_var, .n),
#'                    spp_gap_t = 500L, final_t = 50e3L, n_reps = 12,
#'                    save_every = 0L,
#'                    n_threads = .N_THREADS, show_progress = FALSE)
#'
#'     Z$call[["eta"]] <- eval(.eta)
#'     Z$call[["d"]] <- eval(.d)
#'     Z$call[["n"]] <- eval(.n)
#'     Z$call[["add_var"]] <- eval(rep(.add_var, .n))
#'
#'     jacs <- jacobians(Z)
#'
#'     if (!is.null(.pb)) .pb$tick()
#'
#'     return(list(NV = Z$nv %>%
#'                     mutate(trait = paste0("V", trait)) %>%
#'                     spread(trait, geno) %>%
#'                     mutate(d = .d[1], eta = .eta, add_var = .add_var),
#'                 J = jacs,
#'                 V0 = .V0,
#'                 seed = .seed))
#' }
#'
#'
#' # Simulations varying both d values
#' if (.REDO_SIMS) {
#'     # Takes ~XX min w/ 3 threads
#'     t0 <- Sys.time()
#'     coexist_U_sims <- crossing(.d = seq(-0.25, 2, length.out = 10),
#'                              .eta = c(-1,1) * etas[[2]],
#'                              .add_var = seq(0.01, 0.1, 0.01)) %>%
#'         # bc `seq` makes weird numbers that are ~1e-15 from what they should be:
#'         mutate(across(.fns = round, digits = 2))
#'     pb <- progress_bar$new(format = " [:bar] :percent in :elapsed | eta: :eta",
#'                            total = nrow(coexist_U_sims), clear = FALSE,
#'                            show_after = 0)
#'     set.seed(2014151565)
#'     coexist_U_sims <- coexist_U_sims %>%
#'         pmap(one_coexist_U_combo, .pb = pb)
#'     saveRDS(coexist_U_sims, rds("coexist_U_sims"))
#'     t1 <- Sys.time()
#'     t1 - t0; rm(t0, t1)
#' } else {
#'     coexist_U_sims <- readRDS(rds("coexist_U_sims"))
#' }
