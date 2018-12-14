
# =================================================================
# ==========================

#   This file simulates random datasets to be used to test the accuracy of
#   analytical solutions for the Jacobian matrices

# ==========================
# =================================================================

library(readr)

q <- 3          # number of traits
n <- 4          # number of species
nsims <- 100    # number of simulated datasets

# Column names for output:
names_ <- c(paste0("V", 1:(n*q)), paste0("N", 1:n), "f", "g", "eta", "r0", "d")

sim_data <- function(i) { # `i` is just used so the fxn can be passed to lapply

    V <- rnorm(n*q, 0, 5)
    N <- rlnorm(n, log(10), 0.5)
    f <- rlnorm(1, log(0.2), 0.5)
    g <- rlnorm(1, log(0.1), 0.5)
    eta <- rnorm(1, 0, 0.1)
    r0 <- rlnorm(1, log(0.25), 0.25)
    d <- rnorm(1, 0, 0.2)

    out <- as.data.frame(rbind(c(V, N, f, g, eta, r0, d)))
    colnames(out) <- names_

    return(out)
}

set.seed(78667)
sims <- do.call(rbind, lapply(1:nsims, sim_data))

write_csv(sims, "informal_tests/simulated_data.csv")
