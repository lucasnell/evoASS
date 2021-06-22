#
# Checks inputs to `adapt_dyn`, creates the D and C matrices, and edits
# the `n_threads` value if necessary.
#
#
check_adapt_dyn_args <- function(eta, d, q, n, V0, N0, f, a0, r0,
                                 mut_sd, mut_prob, max_clones,
                                 sigma_V0, sigma_N, sigma_V, n_reps, max_t,
                                 min_N, save_every, show_progress, n_threads) {


    stopifnot(sapply(list(eta, d, q, n, f, a0, r0, sigma_V0, sigma_N, sigma_V,
                          n_reps, max_t, min_N, save_every,
                          mut_sd, mut_prob, max_clones,
                          n_threads, N0), is.numeric))
    stopifnot(sapply(list(q, n, f, a0, r0, sigma_V0, sigma_N,
                          mut_sd, mut_prob, max_clones,
                          n_reps, max_t, save_every, n_threads),
                     length) == 1)
    stopifnot(length(sigma_V) %in% c(1, q))
    if (!is.null(V0)) {
        stopifnot(is.numeric(V0))
        stopifnot(inherits(V0, "matrix") || length(V0) %in% c(1, q))
        stopifnot(all(V0 >= 0))
        if (inherits(V0, "matrix")) stopifnot(nrow(V0) == q && ncol(V0) == n)
    }

    stopifnot(n >= 1 && q >= 1)
    stopifnot(N0 >= 0)
    stopifnot(c(n_reps, max_t, n_threads) >= 1)
    stopifnot(c(save_every, sigma_V0, min_N) >= 0)
    stopifnot(mut_sd > 0)
    stopifnot(mut_prob >= 0 && mut_prob <= 1)

    stopifnot(length(eta) %in% c(1, q^2))
    stopifnot(length(d) %in% c(1, q))

    if (n_threads > 1 && !using_openmp()) {
        message("\nOpenMP not enabled. Only 1 thread will be used.\n")
        n_threads <- 1
    }

    C <- matrix(eta[1], q, q)
    if (length(eta) == q^2) {
        stopifnot(inherits(eta, "matrix") &&
                      identical(dim(eta), as.integer(c(q,q))) &&
                      isSymmetric(eta))
        C <- eta
    }
    diag(C) <- 1

    D <- matrix(0, q, q)
    diag(D) <- d

    return(list(C = C, D = D, n_threads = n_threads))


}



#' Adaptive dynamics.
#'
#'
#' @param V0 Trait value(s) for each starting clone.
#'     For only one starting line, must be a numeric vector or a single
#'     matrix row or column.
#' @param N0 Abundance(s) for each starting clone. Must be a numeric vector or a single
#'     matrix row or column.
#' @param f A single number representing the cost of the trait on the growth rate.
#' @param a0 A single number representing the base density dependence.
#' @param eta Number(s) representing the non-additive effects of traits on the
#'     growth rate.
#'     Should be a single number or a symmetrical, numeric, `q` by `q` matrix.
#' @param r0 A single number representing the base growth rate.
#' @param d Number(s) that adjusts how the focal line is affected by
#'     other lines' trait values.
#'     Should be of length 1 or the same as the number of traits.
#'     If `d < 0`, then increases in `V_j` (trait that reduces competition
#'     experienced by clone `j`) increases competition experienced by clone `i`,
#'     thereby giving conflicting coevolution.
#'     Conversely, if `d > 0`, then increases in `V_j` decrease competition
#'     experienced by clone `i`, leading to nonconflicting coevolution.
#' @param max_t Maximum time simulated.
#' @param min_N Minimum N that's considered extant.
#' @param mut_sd Standard deviation of the normal distribution used to generate
#'     mutated (i.e., daughter) trait values.
#' @param mut_prob Probability of a mutation (creating a daughter lineage).
#' @param show_progress Boolean for whether to show a progress bar.
#' @param max_clones Maximum number of clones predicted. This is used only to reserve
#'     memory for some of the inner C++ objects, so when deciding on a value for this,
#'     you should choose a high value.
#' @param save_every Number of time steps between when saving information for output.
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr as_tibble
#' @importFrom tidyr gather
#' @importFrom dplyr mutate
#' @importFrom dplyr mutate_at
#'
#'
adapt_dyn <- function(
    eta,
    d,
    q,
    n_reps,
    n = 1,
    V0 = 1,
    N0 = rep(1, n),
    f = 0.1,
    a0 = 1e-4,
    r0 = 0.5,
    sigma_V0 = 1,
    sigma_N = 0,
    sigma_V = 0,
    max_t = 1e4L,
    min_N = 1,
    save_every = 10L,
    mut_sd = 0.1,
    mut_prob = 0.01,
    max_clones = 1e4,
    show_progress = TRUE,
    n_threads = 1) {


    call_ <- match.call()
    # So it doesn't show the whole function if using do.call:
    if (call_[1] != as.call(quote(adapt_dyn()))) {
        call_[1] <- as.call(quote(adapt_dyn()))
    }

    args <- check_adapt_dyn_args(eta, d, q, n, V0, N0, f, a0, r0,
                                 mut_sd, mut_prob, max_clones,
                                 sigma_V0, sigma_N, sigma_V, n_reps, max_t,
                                 min_N, save_every, show_progress, n_threads)

    C <- args$C
    D <- args$D
    n_threads <- args$n_threads

    if (length(sigma_V) == 1) sigma_V <- rep(sigma_V, q)

    if (max_clones < 100) max_clones <- 100

    if (is.null(V0)) {
        # Otherwise start at zero:
        if (sigma_V0 == 0 && all(sigma_V) == 0) {
            warning(paste("\nSimulations start with species having all axes",
                          "at zero, and you aren't providing stochasticity in",
                          "either the starting values or via phenotypes.",
                          "Because this is often an unstable equilibrium,",
                          "these simulations may be odd or boring.",
                          "Continuing anyway..."))
        }
        V0 <- matrix(0, q, n)
    } else if (!inherits(V0, "matrix")) {
        stopifnot(length(V0) == 1 || length(V0) == q)
        V0 <- matrix(V0, q, n)
    }


    sim_output <- adapt_dyn_cpp(n_reps = n_reps,
                                V0 = split(t(V0), 1:ncol(V0)),
                                N0 = N0,
                                f = f,
                                a0 = a0,
                                C = C,
                                r0 = r0,
                                D = D,
                                sigma_V0 = sigma_V0,
                                sigma_N = sigma_N,
                                sigma_V = sigma_V,
                                max_t = max_t,
                                min_N = min_N,
                                mut_sd = mut_sd,
                                mut_prob = mut_prob,
                                show_progress = show_progress,
                                max_clones = max_clones,
                                save_every = save_every,
                                n_threads = n_threads)

    colnames(sim_output) <- c("rep", "time", "clone", "N", sprintf("V%i", 1:q))

    if (show_progress) cat("Simulations finished...\n")

    NVt <- as_tibble(sim_output) %>%
        gather("trait", "V", dplyr::starts_with("V")) %>%
        mutate(trait = gsub("V", "", trait)) %>%
        mutate_at(dplyr::vars(rep, time, clone, trait), as.integer)

    ad_obj <- list(data = NVt, call = call_)

    class(ad_obj) <- "adapt_dyn"

    return(ad_obj)

}

#' Print a `adapt_dyn` object.
#'
#' @param x an object of class \code{adapt_dyn}.
#' @param digits the number of digits to be printed.
#' @param ... arguments passed to and from other methods.
#'
#' @export
#' @noRd
#'
print.adapt_dyn <- function(x, digits = max(3, getOption("digits") - 3),
                                    ...) {

    cat(crayon::cyan$bold("# Output from adaptive dynamics\n"))
    print(x$data, digits, ...)

    invisible(x)

}




