
# =================================================================
# ==========================

#   This file simulates random datasets to be used to test the accuracy of
#   analytical solutions for the Jacobian matrices

# ==========================
# =================================================================

library(readr)
library(sauron)

q <- 3          # number of traits
n <- 4          # number of species
nsims <- 100    # number of simulated datasets

# Column names for output:
names_ <- c(paste0("V", 1:(n*q)), paste0("N", 1:n), "f", "a0", "eta", "r0", "d")

sim_data <- function(i) { # `i` is just used so the fxn can be passed to lapply

    N <- rlnorm(n, log(10), 0.5)
    f <- runif(1, 0, 0.5)
    a0 <- runif(1, 0, 0.5)
    eta <- runif(1, -0.5, 0.5)
    r0 <- runif(1, f * (1 + eta), 2)
    d <- runif(1, -0.1, 0.1)
    .sd <- max(as.numeric(unlist(stable_points(eta, f, a0, r0))))
    V <- abs(rnorm(n*q, 0, .sd))

    out <- as.data.frame(rbind(c(V, N, f, a0, eta, r0, d)))
    colnames(out) <- names_

    return(out)
}

set.seed(673349653)
sims <- do.call(rbind, lapply(1:nsims, sim_data))

write_csv(sims, "_check_derivs/simulated_data.csv")
