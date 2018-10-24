

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
#' @param V List of row vectors, each containing trait values at time t (1x2 vector)
#'     for a particular clone.
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


