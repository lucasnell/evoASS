
# Testing that selection strength function (`sel_str_cpp`) works.


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

    df_dVis <- lapply(1:n, function(i) {
        V_i_ <- V_[[i]]
        V_nei_ <- V_[-i]
        N_i_ <- N_[i]
        N_nei_ <- N_[-i]
        as.numeric(evoASS:::dF_dVi_cpp(V_i_, V_nei_, N_i_, N_nei_,
                                       f_, g_, C_, r0_, d_))
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


