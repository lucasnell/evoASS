
#' Create list of arguments for running a set of simulations using `adapt_dyn`.
#'
#' @param eta_sign Sign of the `eta` parameter.
#' @param d_sign Sign of the `d` parameter.
#' @param q Number of traits.
#' @param ... Other arguments to add for `adapt_dyn` function, namely `n_reps`.
#'
#' @export
#'
adapt_dyn_args <- function(eta_sign, d_sign, q, ...) {
    stopifnot(is.numeric(eta_sign), is.numeric(d_sign), is.numeric(q))
    stopifnot(length(eta_sign) == 1, length(d_sign) == 1, length(q) == 1)
    # List for all parameters:
    args <- list(n_cores = 4, q = q)
    # the non-additive effects of traits on `r`:
    args$eta <- 0.01 * sign(eta_sign)
    # changes how the focal line's traits affect other lines' effects of competition:
    args$d <- 0.1 * sign(d_sign)
    # High d takes a very long time to run using quant_gen, so I reduced its value
    # when positive. I want to use the same value here.
    if (args$d > 0) args$d <- 1e-4
    others <- list(...)
    args[names(others)] <- others
    return(args)
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
#' @param g A single number representing the benefit of the trait on the density
#'     dependence.
#' @param eta A single number representing the non-additive effects of traits on the
#'     growth rate.
#' @param r0 A single number representing the base growth rate.
#' @param d A single number that adjusts how the focal line is affected by
#'     other lines' trait values.
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
#' @importFrom dplyr data_frame
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @importFrom purrr map
#' @importFrom purrr map_dbl
#'
#'
adapt_dyn <- function(
    V0,
    N0,
    f = 0.1,
    g = 0.1,
    eta = 0.2,
    r0 = 0.5,
    d = -0.01,
    max_t = 1e4,
    min_N = 1e-4,
    mut_sd = 0.1,
    mut_prob = 0.01,
    show_progress = FALSE,
    max_clones = 1e4,
    save_every = 100) {

    call_ <- match.call()
    # So it doesn't show the whole function if using do.call:
    if (call_[1] != as.call(quote(adapt_dyn()))) {
        call_[1] <- as.call(quote(adapt_dyn()))
    }

    if (inherits(V0, "numeric")) {
        V0 <- list(rbind(V0))
    } else if (inherits(V0, "matrix")) {
        V0 <- lapply(split(V0, 1:nrow(V0)), rbind)
    } else stop("V0 must be numeric or matrix")

    sim_output <- adapt_dyn_cpp(V0, N0, f, g, eta, r0, d, max_t, min_N, mut_sd,
                                mut_prob, show_progress, max_clones, save_every)

    time_pts <- sim_output[["T"]]
    if (length(sim_output$N) != length(time_pts)) {
        warning("should be the same length: ", length(sim_output$N), " ",
                length(time_pts))
        return(sim_output)
    }

    NVt <-
        data_frame(time = 1:length(sim_output$N) %>%
                       map(~ rep(sim_output[["T"]][.x],
                                 length(sim_output$N[[.x]]))) %>%
                       c(recursive = TRUE),
                   clone = c(sim_output$I, recursive = TRUE),
                   N = c(sim_output$N, recursive = TRUE)) %>%
        mutate(V = map(clone, ~ sim_output$V[[.x+1]])) %>%
        unnest() %>%
        group_by(time, clone) %>%
        mutate(trait = paste0("V", 1:n()) %>% factor()) %>%
        ungroup() %>%
        arrange(time, clone, trait) %>%
        dplyr::select(time, clone, trait, everything()) %>%
        identity()

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










#' @describeIn adapt_dyn Add a perturbation to an `adapt_dyn` object.
#'
#' @param obj An `adapt_dyn` object.
#' @param new_prop Proportion of new clones.
#' @param new_trait_means Mean(s) of new clones' trait(s).
#' @param new_N_mean Mean of new clones' abundances.
#' @param new_N_sd SD of new clones' abundances.
#' @param which_traits Traits to change. Defaults to `1:length(new_trait_means)`.
#' @param new_trait_sigmas SDs of new clones' trait(s).
#'     Defaults to SD of each trait at the end of the previous run of
#'     `adapt_dyn()`.
#' @param ... Arguments to pass to `adapt_dyn()`.
#'
#'
#' @export
#'
perturb.adapt_dyn <- function(obj, new_prop,
                              new_trait_means,
                              new_N_mean, new_N_sd,
                              which_traits = NULL,
                              new_trait_sigmas = NULL,
                              ...) {

    stopifnot(new_prop >= 0)

    if (is.null(which_traits)) which_traits <- 1:length(new_trait_means)
    if (length(which_traits) != length(new_trait_means)) {
        stop("\nwhich_traits, and delta_mus must all be the same length.")
    }

    if (ncol(obj$data) != 5) {
        stop("\nad_obj should not be edited before inputting here.")
    }
    old_clones <- obj %>%
        .[["data"]] %>%
        filter(time == max(time)) %>%
        spread(trait, V) %>%
        select(-time, -clone, -N) %>%
        as.matrix()

    if (is.null(new_trait_sigmas)) {
        new_trait_sigmas <- map_dbl(which_traits, ~ sd(old_clones[,.x]))
    } else if (length(which_traits) != length(new_trait_sigmas)) {
        stop("\nnew_trait_sigmas must be NULL or numeric vector the same ",
             "length as which_traits.")
    }

    n_old_clones <- nrow(old_clones)
    n_new_clones <- pmax(round(new_prop * n_old_clones), 1)
    # if (any(new_trait_means < 0)) stop("new_trait_means cannot be < 0")
    new_clones <- matrix(NA_real_, n_new_clones, ncol(old_clones))
    for (i in 1:ncol(new_clones)) {
        if (i %in% which_traits) {
            j <- which(which_traits == i)
            new_clones[,i] <- new_trait_means[j] +
                rnorm(n = n_new_clones) * new_trait_sigmas[j]
        } else {
            new_clones[,i] <- sample(old_clones[,i], n_new_clones, replace = TRUE)
        }
    }

    new_traits <- rbind(old_clones, new_clones)
    new_Ns <- c(obj %>%
                    .[["data"]] %>%
                    filter(time == max(time)) %>%
                    spread(trait, V) %>%
                    .[["N"]],
                trunc_rnorm_cpp(n_new_clones, new_N_mean, new_N_sd))

    new_call <- obj$call
    new_call$V0 <- quote(new_traits)
    new_call$N0 <- quote(new_Ns)
    new_args <- list(...)
    for (x in names(new_args)) {
        new_call[[x]] <- new_args[[x]]
    }
    new_ad <- eval(new_call)

    new_ad$data <- new_ad %>% .[["data"]] %>%
        mutate(time = time + max(obj$data$time))

    obj$data <- obj$data %>% filter(time < max(time))

    obj$data <- bind_rows(obj$data, new_ad$data)

    return(obj)

}
