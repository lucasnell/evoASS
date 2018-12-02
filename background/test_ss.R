
# Testing that selection strength function (`sel_str_cpp`) works.


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
test_ss <- function(seed, q, n) {
    stopifnot(q >= 1 && n >= 2)
    set.seed(seed)
    V_ <- lapply(1:n, function(i) rbind(runif(q, 0, 10)))
    N_ <- exp(rnorm(n, log(10), 0.5))
    f_ <- exp(rnorm(1, log(0.2), 0.5))
    g_ <- exp(rnorm(1, log(0.1), 0.5))
    C_ <- matrix(rnorm(1, 0, 0.1), q, q)
    diag(C_) <- 1
    r0_ <- exp(rnorm(1, log(0.25), 0.25))
    d_ <- rnorm(1, 0, 0.2)
    F_ <- evoASS:::F_t_cpp(V_, N_, f_, g_, C_, r0_, d_)

    # Numerical derivatives:
    df_dVis <- lapply(1:n, function(i) {
        foo <- function(x) {
            V__ <- V_
            V__[[i]] <- rbind(x)
            evoASS:::F_t_cpp(V__, N_, f_, g_, C_, r0_, d_)[[i]]
        }
        return(grad(foo, x = as.numeric(V_[[i]])))
    })

    ss <- evoASS:::sel_str_cpp(V_, N_, f_, g_, C_, r0_, d_)

    diffs <- lapply(1:n, function(i) ss[i,] * F_[i] - df_dVis[[i]])
    diffs <- do.call(rbind, diffs)

    return(diffs)
}

set.seed(94156)
seeds <- sample.int(.Machine$integer.max, 1000, replace = TRUE)

diffs <- lapply(seeds, test_ss, q = 3, n = 4)
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




