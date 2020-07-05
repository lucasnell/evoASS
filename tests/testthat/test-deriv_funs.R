
#'
#' Testing correctness of functions for the following partial derivatives:
#'
#' - Fitness in species i (at time t) with respect to its trait values (also at time t)
#' - Trait values in species i at time t+1 with respect to its trait values
#'   at time t
#' - Trait values in species i at time t+1 with respect to another species' trait values
#'   at time t
#'

# library(sauron)
# library(testthat)

context("derivative functions")


dir <- sprintf("%s/../../_check_derivs/", testthat::test_path())



# --------------------*
# Function to get simulated dataset
# --------------------*
get_sim_info <- function(sim_i) {

    sims <- readr::read_csv(paste0(dir, "simulated_data.csv"),
                            col_types = readr::cols(
                                .default = readr::col_double()))

    info <- list2env(as.list(sims[sim_i,c("f", "a0", "r0", "d", "eta")]))

    with(info, {
        N = as.numeric(sims[sim_i,colnames(sims)[grepl("^N", colnames(sims))]])
        V = matrix(as.numeric(sims[sim_i,colnames(sims)[grepl("^V", colnames(sims))]]),
                   ncol = length(N))
        V_ = lapply(split(V, col(V)), cbind)
        C = matrix(eta, nrow(V), nrow(V))
        diag(C) = 1
        add_var = 0.01
        D <- matrix(0, nrow(V), nrow(V))
        diag(D) <- d
        F_ <- sauron:::F_t_cpp(V_, N, f, a0, C, r0, D)
        Omegas <- sapply(1:n, function(i) {
            Njs <- sapply((1:n)[-i], function(j) {
                Vj <- V[,j,drop=FALSE]
                N[j] * exp(- t(Vj) %*% D %*% Vj)
            })
            return(N[i] + sum(Njs))
        })
    })

    return(as.list(info))

}


# --------------------*
# Functions to calculate sauron versions:
# --------------------*

calc_dF_dVi <- function(sim_info) {
    F_ <- with(sim_info, sauron:::F_t_cpp(V_, N, f, a0, C, r0, D))
    ss <- with(sim_info, {
        sauron:::sel_str_cpp(V_, N, f, a0, C, r0, D)
    })
    sauron_results <- ss %*% diag(as.numeric(F_))
    return(sauron_results)
}
calc_dVi_dVi <- function(sim_info) {
    n <- length(sim_info$N)
    lapply(1:n, function(i) {
        with(sim_info, {
            sauron:::dVi_dVi_cpp(i - 1, V, Omegas[i], C, f, a0, add_var)
        })
    })
}
calc_dVi_dVk <- function(sim_info) {
    n <- length(sim_info$N)
    mats <- lapply(1:n, function(i) {
        lapply((1:n)[1:n != i], function(k) {
            with(sim_info, {
                sauron:::dVi_dVk_cpp(i-1, k-1, N, V, D, a0, add_var)
            })
        })
    })
    unlist(mats, recursive = FALSE)
}

# NEW ONES -----

calc_dVi_dNi <- function(sim_info) {
    n <- length(sim_info$N)
    derivs <- lapply(1:n, function(i) {
        with(sim_info,
             {
                 Vi <- V[,i,drop=FALSE]
                 2 * add_var * a0 * Vi %*% exp(- t(Vi) %*% Vi)
             })
    })
    do.call(cbind, derivs)
}


calc_dVi_dNk <- function(sim_info) {

    n <- length(sim_info$N)
    mats <- lapply(1:n, function(i) {
        Ms <- lapply((1:n)[1:n != i], function(k) {
            with(sim_info,
                 {
                     Vi <- V[,i,drop=FALSE]
                     Vk <- V[,k,drop=FALSE]
                     2 * add_var * a0 * Vi %*% exp(- t(Vk) %*% D %*% Vk - t(Vi) %*% Vi)
                 })
        })
        do.call(cbind, Ms)
    })

    return(mats)

}


calc_dNi_dVi <- function(sim_info) {

    n <- length(sim_info$N)

    derivs <- lapply(1:n, function(i) {
        with(sim_info, {
            Vi <- V[,i,drop=FALSE]
            2 * F_[i] * N[i] * (a0 * Omegas[i] * exp(- t(Vi) %*% Vi) %*% t(Vi) - f * t(Vi) %*% C)
        })
    })

    return(do.call(rbind, derivs))
}


calc_dNi_dVk <- function(sim_info) {

    n <- length(sim_info$N)

    mats <- lapply(1:n, function(i) {
        Ms <- lapply((1:n)[1:n != i], function(k) {
            with(sim_info,
                 {
                     Vi <- V[,i,drop=FALSE]
                     Vk <- V[,k,drop=FALSE]
                     2 * F_[i] * N[i] * N[k] * a0 *
                         exp(- t(Vi) %*% Vi - t(Vk) %*% D %*% Vk) %*%
                         t(Vk) %*% D
                 })
        })
        do.call(rbind, Ms)
    })

    return(mats)

}


calc_dNi_dNi <- function(sim_info) {

    n <- length(sim_info$N)

    sapply(1:n, function(i) {
        with(sim_info, {
            Vi <- V[,i,drop=FALSE]
            F_[i] * (1 - N[i] * a0 * exp(- t(Vi) %*% Vi))
        })
    })

}


calc_dNi_dNk <- function(sim_info) {

    n <- length(sim_info$N)

    mats <- lapply(1:n, function(i) {
        Ms <- lapply((1:n)[1:n != i], function(k) {
            with(sim_info,
                 {
                     Vi <- V[,i,drop=FALSE]
                     Vk <- V[,k,drop=FALSE]
                     - F_[i] * N[i] * a0 *
                         exp(- t(Vi) %*% Vi - t(Vk) %*% D %*% Vk)
                 })
        })
        do.call(cbind, Ms)
    })

    do.call(rbind, mats)

}


check_results <- function(type) {

    # type = "dNi_dNk"

    poss_types <- c("dF_dVi", "dVi_dVi", "dVi_dVk",
                    "dVi_dNi", "dVi_dNk", "dNi_dVi", "dNi_dVk",
                    "dNi_dNi", "dNi_dNk")

    if (!type %in% poss_types) stop("\nERROR: type not recognized")

    py_results_df <- readr::read_csv(sprintf("%sresults/%s.csv", dir, type),
                                     col_names = FALSE,
                                     col_types = readr::cols(
                                         .default = readr::col_double()))

    nsims <- 100
    n <- 4
    q <- 3


    same_results <- logical(nsims)

    # sim_i = 1
    for (sim_i in 1:nsims) {


        py_results <- as.numeric(py_results_df[sim_i,])
        info <- get_sim_info(sim_i)

        if (type == "dF_dVi") {

            py_results <- matrix(py_results, q, n)
            sauron_results <- calc_dF_dVi(info)

        } else if (type == "dVi_dVi") {

            py_results <- lapply(split(py_results, rep(1:n, each = q^2)),
                                 matrix, nrow = q, ncol = q, byrow = TRUE)
            sauron_results <- calc_dVi_dVi(info)

        } else if (type == "dVi_dVk") {

            py_results <- lapply(split(py_results, rep(1:(n * (n-1)), each = q^2)),
                                 matrix, nrow = q, ncol = q, byrow = TRUE)
            sauron_results <- calc_dVi_dVk(info)

        } else if (type =="dVi_dNi") {

            py_results <- matrix(py_results, q, n)
            sauron_results <- calc_dVi_dNi(info)

        } else if (type =="dVi_dNk") {

            py_results <- lapply(split(py_results, rep(1:n, each = q * (n-1))),
                                 matrix, nrow = q)
            sauron_results <- calc_dVi_dNk(info)

        } else if (type =="dNi_dVi") {

            py_results <- matrix(py_results, nrow = n, byrow = TRUE)
            sauron_results <- calc_dNi_dVi(info)

        } else if (type =="dNi_dVk") {

            py_results <- lapply(split(py_results, rep(1:n, each = q * (n-1))),
                   matrix, ncol = q, byrow = TRUE)
            sauron_results <- calc_dNi_dVk(info)

        } else if (type =="dNi_dNi") {

            # py_results doesn't need to be manipulated
            sauron_results <- calc_dNi_dNi(info)

        } else if (type =="dNi_dNk") {

            py_results <- matrix(py_results, nrow = n, byrow = TRUE)
            sauron_results <- calc_dNi_dNk(info)

        } else {
            stop("\nERROR: type not recognized")
        }

        same_results[sim_i] <- isTRUE(all.equal(sauron_results, py_results,
                                                check.attributes = FALSE))

    }


    return(sum(same_results))
}




expect_equal(check_results("dF_dVi"), 100)
expect_equal(check_results("dVi_dVi"), 100)
expect_equal(check_results("dVi_dVk"), 100)
# NEW ONES ----
expect_equal(check_results("dVi_dNi"), 100)
expect_equal(check_results("dVi_dNk"), 100)
expect_equal(check_results("dNi_dVi"), 100)
expect_equal(check_results("dNi_dVk"), 100)
expect_equal(check_results("dNi_dNi"), 100)
expect_equal(check_results("dNi_dNk"), 100)
