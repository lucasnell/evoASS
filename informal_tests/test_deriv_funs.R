
#'
#' Testing correctness of functions for the following partial derivatives:
#'
#' - Fitness in species i (at time t) with respect to its trait values (also at time t)
#' - Trait values in species i at time t+1 with respect to its trait values
#'   at time t
#' - Trait values in species i at time t+1 with respect to another species' trait values
#'   at time t
#'



# --------------------*
# Function to get simulated dataset
# --------------------*
get_sim_info <- function(sim_i) {

    sims <- readr::read_csv("informal_tests/simulated_data.csv",
                            col_types = readr::cols(.default = readr::col_double()))

    info <- list2env(as.list(sims[sim_i,c("f", "g", "r0", "d", "eta")]))

    with(info, {
        N = as.numeric(sims[sim_i,colnames(sims)[grepl("^N", colnames(sims))]])
        V = matrix(as.numeric(sims[sim_i,colnames(sims)[grepl("^V", colnames(sims))]]),
                   length(N))
        C = matrix(eta, ncol(V), ncol(V))
        diag(C) = 1
        CCC = C + t(C)
        sigma2 = 0.01
    })

    return(as.list(info))
}


# --------------------*
# Functions to calculate sauron versions:
# --------------------*

calc_dF_dVi <- function(sim_info) {
    F_ <- with(sim_info, {
        sauron:::F_t_cpp(split(V, row(V)), N, f, g, C, r0, d)
    })
    ss <- with(sim_info, {
        sauron:::sel_str_cpp(split(V, row(V)), N, f, g, C, r0, d)
    })
    sauron_results <- diag(as.numeric(F_)) %*% ss
    return(sauron_results)
}
calc_dVi_dVi <- function(sim_info) {
    n <- length(sim_info$N)
    lapply(1:n, function(i) {
        Z <- with(sim_info, {
            N[i] + sum(sapply(1:length(N),
                              function(j) {
                                  if (j == i) return(0)
                                  exp(-d * V[j,,drop=F] %*% t(V[j,,drop=F])) * N[j]
                              }))
        })
        with(sim_info, {
            sauron:::dVi_dVi(i - 1, V, Z, CCC, f, g, sigma2)
        })
    })
}
calc_dVi_dVk <- function(sim_info) {
    n <- length(sim_info$N)
    mats <- lapply(1:n, function(i) {
        lapply((1:n)[1:n != i], function(k) {
            with(sim_info, {
                sauron:::dVi_dVk(i-1, k-1, N, V, d, g, sigma2)
            })
        })
    })
    unlist(mats, recursive = FALSE)
}



check_results <- function(type) {

    if (!type %in% c("dF_dVi", "dVi_dVi", "dVi_dVk")) stop("type not recognized")

    py_results_df <- readr::read_csv(sprintf("informal_tests/results/%s.csv", type),
                                     col_names = FALSE,
                                     col_types = readr::cols(.default = readr::col_double()))

    nsims <- nrow(py_results_df)

    info <- get_sim_info(1)

    n <- length(info$N)
    q <- ncol(info$V)

    same_results <- logical(nsims)

    for (sim_i in 1:nsims) {

        py_results <- as.numeric(py_results_df[sim_i,])
        info <- get_sim_info(sim_i)

        if (type == "dF_dVi") {
            py_results <- matrix(py_results, n, q, byrow = TRUE)
            sauron_results <- calc_dF_dVi(info)
        } else if (type == "dVi_dVi") {
            py_results <- lapply(split(py_results, rep(1:n, each = q^2)),
                                 matrix, nrow = q, ncol = q, byrow = TRUE)
            sauron_results <- calc_dVi_dVi(info)
        } else {
            py_results <- lapply(split(py_results, rep(1:(n * (n-1)), each = q^2)),
                                 matrix, nrow = q, ncol = q, byrow = TRUE)
            sauron_results <- calc_dVi_dVk(info)
        }

        same_results[sim_i] <- isTRUE(all.equal(sauron_results, py_results,
                                                check.attributes = FALSE))

    }

    return(sum(same_results) == nsims)
}



sapply(c("dF_dVi", "dVi_dVi", "dVi_dVk"), check_results)
