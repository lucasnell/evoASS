
#'
#'
#' This organizes the outcomes from `giant_inv_sims`.
#'
#'


suppressPackageStartupMessages({
    library(magrittr)
    library(dplyr)
    library(purrr)
    library(purrr)
    library(parallel)
})

options(dplyr.summarise.inform = FALSE)

.N_THREADS <- 12


# Where to find input rds files
in_rds <- function(.x) {
    sprintf("./giant_inv_sims/%s.rds",
            gsub("\\.rds$", "", gsub("^\\/", "", .x)))
}



giant_inv_sims <- map_dfr(0:17,
                          function(i) {
                              fn <- paste0("giant_inv_sims_", i)
                              X <- readRDS(in_rds(fn))
                              extract_qg <- function(j) {
                                  X[["qg"]][[j]]$nv %>%
                                      filter(trait == 1) %>%
                                      group_by(rep) %>%
                                      summarize(inv = any(spp == 2, na.rm = TRUE),
                                                res = any(spp == 1, na.rm = TRUE),
                                                ext = any(is.na(spp)),
                                                total = 1) %>%
                                      ungroup() %>%
                                      select(-rep) %>%
                                      summarize(across(.fns = sum))
                              }
                              # This part takes ~ 4 min per `i` w/ 3 threads
                              X[["qg"]] <- mclapply(1:nrow(X), extract_qg,
                                                    mc.cores = .N_THREADS)
                              X <- X %>%
                                  unnest(qg)
                              return(X)
                          })

saveRDS(giant_inv_sims, "giant_inv_sims_outcomes.rds")
