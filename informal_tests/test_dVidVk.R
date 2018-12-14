
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



one_analytical <- function(i, k, par_list) {

    an_deriv <- with(par_list, {
        n <- length(N)
        q <- ncol(V)
        Vi <- V[i,,drop=FALSE]
        Ni <- N[i]
        Vk <- V[k,,drop=FALSE]
        Nk <- N[k]
        Z <- sum(sapply(1:n,
                        function(j) {
                            if (j == i || j == k) return(0)
                            Vj = V[j,,drop=FALSE]
                            exp(-d * Vj %*% t(Vj)) * N[j]
                        }))

        -4 * sigma2 * Nk * d * g * {
            (t(Vk) %*% exp(-d * Vk %*% t(Vk))) %*% (exp(- Vi %*% t(Vi)) %*% Vi)
        }
    })

    return(an_deriv)
}
one_numeric <- function(i, k, par_list, ...) {

    nd_pars <- list(...)

    num_deriv <- with(par_list, {
        n <- length(N)
        q <- ncol(V)
        Z <- sum(sapply(1:n,
                       function(j) {
                           if (j == i || j == k) return(0)
                           Vj = V[j,,drop=FALSE]
                           exp(-d * Vj %*% t(Vj)) * N[j]
                       }))
        Vi <- V[i,,drop=FALSE]
        Ni <- N[i]
        Vk <- V[k,,drop=FALSE]
        Nk <- N[k]

        foo <- function(x) {
            Vk_ <- rbind(x)
            V_hat <- Vi + sigma2 * {
                (Ni + Nk * exp(-d * Vk_ %*% t(Vk_)) + Z) %*% (2 * g * exp(- Vi %*% t(Vi)) %*% Vi) -
                    f * Vi %*% CCC
            }
            return(as.numeric(V_hat))
        }
        nd_pars <- c(list(func = foo, x = Vk), nd_pars)
        do.call(numDeriv::jacobian, nd_pars)
    })

    return(num_deriv)
}


sim_i <- 1

info <- get_sim_info(sim_i)

one_analytical(1, 2, info)
one_numeric(1, 2, info)

# method.args=list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2)


ans <- do.call(c, lapply(1:length(xN), one_analytical, N = xN, V = xV, f = xf, g = xg,
                         CCC = xCCC, d = xd, sigma2 = xsigma2))
nums <- do.call(c, lapply(1:length(xN), one_numeric, N = xN, V = xV, f = xf, g = xg,
                          CCC = xCCC, d = xd, sigma2 = xsigma2))


#
cat("\n")






























V = rbind(c(3.23531713336706, 1.02777096442878, 0.2787772170268),
          c(5.27526782592759, 1.44656382501125, 0.14410174684599),
          c(9.97409204719588, 1.06037032557651, 2.78145811287686),
          c(5.59907909948379, 2.24725058535114, 1.8701455486007))
N = c(9.4145709397, 7.34398103217961, 19.0046921810612, 10.0717215593192)
f = 0.2484041452;
g = 0.1476385485;
CC = rbind(c(1, -0.0450021495231871, -0.0450021495231871),
            c(-0.0450021495231871, 1, -0.0450021495231871),
            c(-0.0450021495231871, -0.0450021495231871, 1))
r0 = 0.4257492712;
d = -0.0917747103;
n = nrow(V);
sigma2 = 1;
CCC = CC + t(CC);
eye = diag(ncol(Vnei))

diffs <- lapply(1:n, function(i) {
    Z = N[i] + sum(sapply(1:n,
                          function(j) {
                              if (j == i) return(0)
                              exp(-d * V[j,,drop=F] %*% t(V[j,,drop=F])) * N[j]
                          }))
    Vi <- V[i,,drop=FALSE]

    an_deriv <- eye + sigma2 * { 2 * g * Z * exp(- Vi %*% t(Vi))[[1]] *
            (eye + -2 * t(Vi) %*% Vi) - f * CCC }

    foo <- function(x) {
        Vi_ <- rbind(x)
        V_hat <- Vi_ + sigma2 * (-f * Vi %*% CCC + 2 * g * Z * Vi_ *
                                     exp(- Vi_ %*% t(Vi_))[[1]])
        return(as.numeric(V_hat))
    }
    num_deriv <- jacobian(foo, x = Vi)

    # diffs_i <- as.numeric(an_deriv - num_deriv)
    diffs_i <- cbind(as.numeric(an_deriv), as.numeric(num_deriv))

    return(diffs_i)
})

diffs




# q is # traits, n is # species
# test_ss <- function(seed, q, n) {
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
    F_ <- sauron:::F_t_cpp(V_, N_, f_, g_, C_, r0_, d_)

    # Numerical derivatives:
    df_dVis <- lapply(1:n, function(i) {
        foo <- function(x) {
            V__ <- V_
            V__[[i]] <- rbind(x)
            sauron:::F_t_cpp(V__, N_, f_, g_, C_, r0_, d_)[[i]]
        }
        return(grad(foo, x = as.numeric(V_[[i]])))
    })

    ss <- sauron:::sel_str_cpp(V_, N_, f_, g_, C_, r0_, d_)

    diffs <- lapply(1:n, function(i) ss[i,] * F_[i] - df_dVis[[i]])
    diffs <- do.call(rbind, diffs)

    return(diffs)
# }





# ======================================================================================*
# ======================================================================================*
# ======================================================================================*
# ======================================================================================*
# ======================================================================================*
# ======================================================================================*
# ======================================================================================*
# ======================================================================================*
# ======================================================================================*
# ======================================================================================*









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
    F_ <- sauron:::F_t_cpp(V_, N_, f_, g_, C_, r0_, d_)

    # Numerical derivatives:
    df_dVis <- lapply(1:n, function(i) {
        foo <- function(x) {
            V__ <- V_
            V__[[i]] <- rbind(x)
            sauron:::F_t_cpp(V__, N_, f_, g_, C_, r0_, d_)[[i]]
        }
        return(grad(foo, x = as.numeric(V_[[i]])))
    })

    ss <- sauron:::sel_str_cpp(V_, N_, f_, g_, C_, r0_, d_)

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




