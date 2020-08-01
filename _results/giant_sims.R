
#'
#' This file does a "giant"  number of simulations across multiple parameters.
#'
#' In all these simulations, both `d` values are changed.
#'
#'

library(sauron)
library(dplyr)
library(tidyr)
library(purrr)
library(magrittr)

.N_THREADS <- 2

#'
#' These simulations are split into 9 sets.
#' Which set this one should conduct depends on the input to this script.
#'
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[[1]]) + 1L


seeds <- c(1909893419, 554504146, 884616499, 803234389, 1158971242,
           934673902, 1632225031, 1867003471, 651436069)


one_sim_combo <- function(.d, .eta, .add_var, .sigma_N, .sigma_V) {

    .n <- 100

    .seed <- sample.int(2^31 - 1, 1)
    set.seed(.seed)

    Z <- quant_gen(q = 2, eta = .eta, d = .d, n = .n,
                   add_var = rep(.add_var, .n),
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

giant_sims <- giant_sims %>%
    .[i_rows,]

set.seed(seeds[i])
giant_sims <- giant_sims %>%
    pmap(one_sim_combo)


saveRDS(giant_sims, sprintf("giant_sims_%i.rds", i - 1L))



#'
#' Installing `sauron` on cluster
#' =================
#'
#'
#' I mostly just followed the instructions at
#' <http://chtc.cs.wisc.edu/r-jobs.shtml>.
#'
#' To install dependencies I ran this in R:
#'
#' ```r
#' for (p in c("crayon", "dplyr", "ggplot2", "magrittr", "purrr", "Rcpp",
#'             "readr", "tibble", "tidyr", "RcppArmadillo", "RcppProgress")) {
#'     install.packages(p)
#' }
#' ```
#'
#' Because `tidyr` version `1.1.1` (the most recent as of this writing)
#' didn't install properly, I installed an older version:
#'
#' ```r
#' install.packages(paste0("https://cran.r-project.org/src/contrib/Archive/",
#'                         "tidyr/tidyr_1.1.0.tar.gz"),
#'     repos = NULL, type = "source")
#' ```
#'
#' Then I un-tared the `sauron` package `.tar.gz` file and commented the
#' lines under the `Using OpenMP` heading inside `src/Makevars`.
#' I replaced them with...
#'
#' ```
#' PKG_CXXFLAGS += $(SHLIB_OPENMP_CXXFLAGS)
#' PKG_CFLAGS += $(SHLIB_OPENMP_CFLAGS)
#' PKG_LIBS += $(SHLIB_OPENMP_CFLAGS)
#' ```
#'
#' Lastly I did this to install sauron (WITHOUT updating `tidyr`):
#'
#' ```r
#' devtools::install("sauron")
#' ```
