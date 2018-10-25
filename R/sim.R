

#' r (growth rate) based on Vi (Vi = traits for clone i).
#'
#' @noRd
#'
r_V <- function(Vi, f, C, r0) {
    r <- r0 - f * as.numeric(Vi %*% C %*% t(Vi));
    return(r);
}


#' A (interspecific + intraspecific density dependence) for line i based on V and N
#'
#' @noRd
#'
A_i_VN <- function(V_i, V_nei, N_i, N_nei, g, d) {

    # Values of sum of squared trait values for each clone
    W <- numeric(nrow(V_nei) + 1)
    W[1] <- as.numeric(V_i %*% t(V_i))
    for (j in 1:nrow(V_nei)) {
        W[j+1] <- as.numeric(V_nei[j,,drop=FALSE] %*% t(V_nei[j,,drop=FALSE]))
    }

    # Effects of intra- and inter-specific competition
    A <- g * exp(-1 * W[1]) * N_i; # intra
    for (j in 1:nrow(V_nei)) {
        A <- A + g * exp(-1 * (W[1] + d * W[j+1])) * N_nei[j]
    }

    return(A)
}


#' A (interspecific + intraspecific density dependence) based on V and N
#'
#' @noRd
#'
A_VN <- function(V, N, g, d) {

    A <- numeric(nrow(V))

    # Values of sum of squared trait values for each clone
    W <- numeric(nrow(V));
    for (j in 1:nrow(V)) {
        W[j]  <- as.numeric(V[j,,drop=FALSE] %*% t(V[j,,drop=FALSE]));
    }

    for (i in 1:nrow(V)) {

        # Effects of intra- and inter-specific competition
        intra <- g * exp(-1 * W[i]) * N[i];
        inter = 0;
        for (j in 1:nrow(V)) {
            if (j == i) next;
            inter = inter + (g * exp(-1 * (W[i] + d * W[j])) * N[j]);
        }
        A[i] <- intra + inter
    }

    return(A)
}



#' Fitness at time t.
#'
#' @param V Matrix (dimensions: n x 2) of trait values at time t.
#' @param N Row vector of population abundances at time t.
#' @param f Effect of traits on growth rate.
#' @param g Effect of traits on density dependence.
#' @param C Matrix containing non-additive effects of traits on growth rate.
#' @param r0 Starting growth rate.
#' @param d Changes how the focal line is affected by other lines' trait values.
#'     If `d < 0`, then increases in `V_j` (trait that reduces competition
#'     experienced by clone `j`) increases competition experienced by clone `i`,
#'     thereby giving conflicting coevolution.
#'     Conversely, if `d > 0`, then increases in `V_j` decrease competition
#'     experienced by clone `i`, leading to nonconflicting coevolution.
#'
#' @export
#'
F_t <- function(V, N, f, g, C, r0, d) {

    A <- A_VN(V, N, g, d)
    F_ <- numeric(nrow(V))

    for (i in 1:nrow(V)) {
        r <- r_V(V[i,,drop=FALSE], f, C, r0)
        F_[i] = exp(r - A[i])
    }

    return(F_)
}


#' Same as above, but for computing partial derivatives.
#'
#' @param V_i Row vector of traits for focal clone.
#' @param V_nei List of row vectors of traits for non-focal clones.
#' @param N_i Abundance for focal clone.
#' @param N_nei Row vector of abundances for non-focal clones.
#' @inheritParams F_t
#'
#' @noRd
#'
F_t_deriv <- function(V_i, V_nei, N_i, N_nei, f, g, C, r0, d) {

    A <- A_i_VN(V_i, V_nei, N_i, N_nei, g, d)
    r <- r_V(V_i, f, C, r0)

    F_ = exp(r - A)

    return(F_);

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
#' @param min_N Minimum N that's considered extinct.
#' @param mut_sd Standard deviation for generating mutated trait values.
#' @param mut_prob Probability of a mutation.
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
adaptive_dynamics <- function(
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
    if (call_[1] != as.call(quote(adaptive_dynamics()))) {
        call_[1] <- as.call(quote(adaptive_dynamics()))
    }

    if (inherits(V0, "numeric")) {
        V0 <- list(rbind(V0))
    } else if (inherits(V0, "matrix")) {
        V0 <- lapply(split(V0, 1:nrow(V0)), rbind)
    } else stop("V0 must be numeric or matrix")

    sim_output <- adaptive_dynamics_(V0, N0, f, g, eta, r0, d, max_t, min_N, mut_sd,
                                     mut_prob, show_progress, max_clones, save_every)

    time_pts <- sim_output[["T"]]
    if (length(sim_output$N) != length(time_pts)) {
        warning("should be the same length: ", length(sim_output$N), " ",
                length(time_pts))
        return(sim_output)
    }

    NVt <-
        data_frame(time = 1:length(sim_output$N) %>%
                       map(~ rep(time_pts[.x], length(sim_output$N[[.x]]))) %>%
                       c(recursive = TRUE),
                   clone = c(sim_output$I, recursive = TRUE),
                   N = c(sim_output$N, recursive = TRUE)) %>%
        mutate(V1 = map_dbl(clone, ~ sim_output$V[[.x+1]][1]),
               V2 = map_dbl(clone, ~ sim_output$V[[.x+1]][2])) %>%
        gather("trait", "value", V1:V2, factor_key = TRUE) %>%
        mutate(trait = gsub("V", "", trait) %>%
                   as.integer() %>%
                   factor()) %>%
        arrange(time, clone, trait) %>%
        identity()

    ad_obj <- list(data = NVt, call = call_)

    class(ad_obj) <- "adaptive_dynamics"

    return(ad_obj)

}

#' Print a `adaptive_dynamics` object.
#'
#' @param x an object of class \code{adaptive_dynamics}.
#' @param digits the number of digits to be printed.
#' @param ... arguments passed to and from other methods.
#'
#' @export
#' @noRd
#'
print.adaptive_dynamics <- function(x, digits = max(3, getOption("digits") - 3),
                                    ...) {

    cat("< Output from adaptive dynamics >\n")
    print(x$data, digits, ...)

}

