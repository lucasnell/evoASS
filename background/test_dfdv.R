
suppressPackageStartupMessages({
    library(numDeriv)
    library(ggplot2)
    library(dplyr)
    library(purrr)
    library(tidyr)
})
# Used on my local machine for better plotting:
# source(".Rprofile")

# q is # traits, n is # species
test_dfdVi <- function(seed, q, n) {
    stopifnot(q >= 1 && n >= 2)
    set.seed(seed)
    V_i_ <- rbind(runif(q, 0, 10))
    V_nei_ <- lapply(1:(n-1), function(i) rbind(runif(q, 0, 10)))
    N_i_ <- exp(rnorm(1, log(10), 0.3))
    N_nei_ <- exp(rnorm(n-1, log(10), 0.5))
    f_ <- exp(rnorm(1, log(0.2), 0.5))
    g_ <- exp(rnorm(1, log(0.1), 0.5))
    C_ <- matrix(rnorm(1, 0, 0.1), q, q)
    diag(C_) <- 1
    r0_ <- exp(rnorm(1, log(0.25), 0.25))
    d_ <- rnorm(1, 0, 0.2)
    foo <- function(x) {
        evoASS:::F_t_deriv_cpp(V_i = rbind(x), V_nei_, N_i_, N_nei_, f_, g_, C_, r0_, d_)
    }
    num_diff <- grad(foo, x = as.numeric(V_i_))
    sym_diff <- evoASS:::dF_dVi_cpp(V_i_, V_nei_, N_i_, N_nei_, f_, g_, C_, r0_, d_)
    return(num_diff - sym_diff)
}

set.seed(78978906)
seeds <- sample.int(.Machine$integer.max, 1000, replace = TRUE)

diffs <- lapply(seeds, test_dfdVi, q = 3, n = 4)
diffs <- do.call(rbind, diffs)

# Highest |difference|:
diffs[abs(diffs) == max(abs(diffs))]

# Plot of differences:
as_data_frame(diffs) %>%
    set_names(paste(1:3)) %>%
    gather("column", "diff", factor_key = TRUE) %>%
    ggplot(aes(column, diff, color = column)) +
    theme_minimal() +
    theme(strip.text = element_blank(), strip.background = element_blank()) +
    facet_wrap(~ column, scales = "free") +
    geom_jitter(shape = 1, alpha = 0.5) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    NULL


