
suppressPackageStartupMessages({
    library(tidyverse)
    library(sauron)
    library(grid)
})

source(".Rprofile")



save_plot <- function(plot_obj, .width, .height, .prefix = NULL, .suffix = NULL) {
    fn <- gsub("_p$", "", paste(substitute(plot_obj)))
    if (!is.null(.prefix)) fn <- paste0(.prefix, "_", fn)
    if (!is.null(.suffix)) fn <- paste0(fn, "_", .suffix)
    fn <- sprintf("_simulations/zz-comm_plots/%s.pdf", fn)
    cat(fn, "\n")
    cairo_pdf(fn, width = .width, height = .height)
    plot(plot_obj)
    dev.off()
    invisible(NULL)
}

unq_spp_filter <- function(..., .prec = 0.001) {
    V <- list(...)
    V <- do.call(cbind, V)
    V <- split(V, row(V))
    return(as.logical(sauron:::unq_spp_cpp(V, precision = .prec)))
}


.q <- 3
etas <- runif(.q * ((.q - 1) / 2), 0.1, 0.4)

eta <- matrix(1/2, .q, .q)
eta[lower.tri(eta)] <- abs(etas) * sample(c(-1, 1), .q, TRUE)
eta[upper.tri(eta)] <- 0
eta <- eta + t(eta)

trait_ts <- quant_gen(q = .q, eta = eta, d = 0, max_t = 20e3L, n_reps = 24,
                      save_every = 0L, n = 100, N0 = rep(1, 100),
                      start_t = 0, perturb_sd = 2, n_cores = 3)

trait_ts$nv %>%
    mutate(trait = paste0("V", trait)) %>%
    spread(trait, value) %>%
    filter(unq_spp_filter(V1, V2, V3))

colSums(eta - diag(.q)); cat("\n"); eta


