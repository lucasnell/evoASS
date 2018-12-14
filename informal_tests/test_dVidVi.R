
# Testing correctness of function for the partial derivative of trait values
# with respect to other trait values.
# This include trait values for the focal species and non-focal species.

suppressPackageStartupMessages({
    library(numDeriv)
    library(ggplot2)
    library(dplyr)
    library(purrr)
    library(tidyr)
})
source("informal_tests/get_sim_info.R")

one_analytical <- function(i, par_list) {

    an_deriv <- with(par_list, {
        n <- length(N)
        q <- ncol(V)
        eye <- diag(q)
        Z = N[i] + sum(sapply(1:n,
                              function(j) {
                                  if (j == i) return(0)
                                  exp(-d * V[j,,drop=F] %*% t(V[j,,drop=F])) * N[j]
                              }))
        Vi <- V[i,,drop=FALSE]

        eye + sigma2 * { 2 * g * Z * exp(- Vi %*% t(Vi))[[1]] *
                (eye + -2 * t(Vi) %*% Vi) - f * CCC }
    })

    return(an_deriv)
}
one_numeric <- function(i, par_list) {

    num_deriv <- with(par_list, {
        n <- length(N)
        q <- ncol(V)
        eye <- diag(q)
        Z = N[i] + sum(sapply(1:n,
                              function(j) {
                                  if (j == i) return(0)
                                  exp(-d * V[j,,drop=F] %*% t(V[j,,drop=F])) * N[j]
                              }))
        Vi <- V[i,,drop=FALSE]

        foo <- function(x) {
            Vi_ <- rbind(x)
            V_hat <- Vi_ + sigma2 * (-f * Vi %*% CCC + 2 * g * Z * Vi_ *
                                         exp(- Vi_ %*% t(Vi_))[[1]])
            return(as.numeric(V_hat))
        }
        numDeriv::jacobian(foo, x = Vi)
    })

    return(num_deriv)
}




sim_i <- 1
# bi <- 34:35
info <- get_sim_info(sim_i)

one_analytical(i = 1, info)
one_numeric(i = 1, info)

ans <- do.call(c, lapply(1:length(info$N), one_analytical, par_list = info))
nums <- do.call(c, lapply(1:length(info$N), one_numeric, par_list = info))

ans  # [bi]
nums  # [bi]

