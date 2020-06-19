
# i
# 0
#
# N
# [ 7.18015044  5.02769558 14.05186271  7.84885382]
#
# V
# [[4.94511003 5.62953205 4.26294319]
#     [2.86919893 1.8334118  1.88296186]
#     [6.74712574 3.67306824 2.71730932]
#     [6.14252194 3.10081916 2.42625001]]
#
# O
# 7040.411110648025
#
# C
# [[ 1.         -0.33114977 -0.33114977]
#     [-0.33114977  1.         -0.33114977]
#     [-0.33114977 -0.33114977  1.        ]]
#
# f
# 0.06888985494151711
#
# a0
# 0.11211295961402357
#
# s2
# 0.01
#
#
# auto
# [[9.98622203e-01 4.56257198e-04 4.56257198e-04]
#     [4.56257198e-04 9.98622203e-01 4.56257198e-04]
#     [4.56257198e-04 4.56257198e-04 9.98622203e-01]]
#
# sym
# [[9.98622203e-01 4.56257198e-04 4.56257198e-04]
#     [4.56257198e-04 9.98622203e-01 4.56257198e-04]
#     [4.56257198e-04 4.56257198e-04 9.98622203e-01]]




reticulate::use_python("/usr/local/bin/python3")
reticulate::source_python("_simulations/__derivs.py")


i = 0L
N = c(7.18015044, 5.02769558, 14.05186271, 7.84885382)

V = rbind(c(-1, 5.62953205, 4.26294319),
          c(2.86919893, 1.8334118 , 1.88296186),
          c(6.74712574, 3.67306824, 2.71730932),
          c(6.14252194, 3.10081916, 2.42625001))
O = 7040.411110648025
C = rbind(c(1, -0.33114977, -0.33114977),
          c(-0.33114977, 1, -0.33114977),
          c(-0.33114977, -0.33114977, 1))
f = 0.06888985494151711
a0 = 0.11211295961402357
s2 = 0.01

D <- matrix(0, 3, 3)
diag(D) <- 0.1

(dVi = dVi_dVi(i, V, O, C, f, a0, s2))

(dVk = dVi_dVk(i, 2, N, V, D, f, a0, C, s2))

library(sauron)
library(tidyverse)

#'
#' Scenarios from `__rough-results.Rmd` that weren't stable (e_min is min
#' leading eigenvalue from 24 reps):
#'
#' ```
#'   sign1 sign2 sign3     e_min     e_max
#' 1    -1    -1    -1 0.9963916 0.9963916
#' 2    -1    -1     0 0.9995038 0.9995038
#' 3    -1     0    -1 0.9984768 0.9984768
#' 4    -1     1    -1 0.9999428 0.9999428
#' 5     0    -1    -1 0.9987687 0.9987687
#' ```
#'

.q <- 3
set.seed(1652339501)
etas <- runif(.q * ((.q - 1) / 2), 0.1, 0.4)

unq_spp_filter <- function(..., .prec = 0.001) {
    V <- list(...)
    V <- do.call(cbind, V)
    V <- split(V, row(V))
    return(as.logical(sauron:::unq_spp_cpp(V, precision = .prec)))
}


get_jac <- function(r, obj) {

    Z <- obj$nv %>%
        filter(rep == r) %>%
        mutate(trait = paste0("V", trait)) %>%
        spread(trait, value)

    VV <- Z %>%
        select(starts_with("V")) %>%
        as.matrix()
    NN <- Z$N

    n <- length(NN)
    q <- ncol(VV)

    D <- matrix(0, q, q)
    diag(D) <- eval(obj$call[["d"]])

    jac <- sauron:::jacobian_cpp(lapply(split(VV, 1:n), rbind),
                                 N = NN,
                                 f = formals(quant_gen)[["f"]],
                                 a0 = formals(quant_gen)[["a0"]],
                                 D = D,
                                 C = eval(obj$call[["eta"]]),
                                 add_var = eval(formals(quant_gen)[["add_var"]]))

    return(jac)
}


C <- matrix(0, .q, .q)
C[lower.tri(C)] <- abs(etas) * c(1, -1, 1)
C <- C + t(C)
diag(C) <- 1


trait_to <- quant_gen(q = .q, eta = C, d = -0.1, max_t = 20e3L, n_reps = 24,
                      save_every = 0L, n = 100, N0 = rep(1, 100),
                      start_t = 0, perturb_sd = 2, n_threads = 3)


jacs <- lapply(1:24, get_jac, obj = trait_to)

map_dbl(jacs, ~ max(Re(eigen(.x)[["values"]])))

trait_to$nv %>%
    mutate(trait = paste0("V", trait)) %>%
    spread(trait, value) %>%
    filter(unq_spp_filter(V1, V2, V3))



C2 <- matrix(0, .q, .q)
C2[lower.tri(C2)] <- mean(abs(etas)) * c(-1, -1, -1)
C2 <- C2 + t(C2)
diag(C2) <- 1


trait_to2 <- quant_gen(q = .q, eta = C2, d = 1, max_t = 20e3L, n_reps = 24,
                      save_every = 0L, n = 100, N0 = rep(1, 100),
                      start_t = 0, perturb_sd = 2, n_threads = 3)


jacs2 <- lapply(1:24, get_jac, obj = trait_to2)

map_dbl(jacs2, ~ max(Re(eigen(.x)[["values"]])))


trait_to2$nv %>%
    mutate(trait = paste0("V", trait)) %>%
    spread(trait, value) %>%
    filter(unq_spp_filter(V1, V2, V3))






library(pbmcapply)

one_i <- function(i) {

    OO <- NN[i] + sum(sapply((1:length(NN))[-i], function(j) {
        Vj <- VV[j,,drop=FALSE]
        exp(-1 * Vj %*% D %*% t(Vj)) * NN[j]
    }))

    mats <- lapply(1:n, function(k) {

        if (k == i) {

            dV = dVi_dVi(i = i - 1, V = VV, O = OO, C = C, f = 0.1,
                         a0 = 0.5, s2 = 0.1)

        } else {

            dV = dVi_dVk(i = i - 1, k = k - 1, N = NN, V = VV, D = D, f = 0.1,
                         a0 = 0.5, C = C, s2 = 0.1)

        }

        return(dV)

    })

    M <- do.call(cbind, mats)

    return(M)

}

jac_mats <- pbmclapply(1:n, one_i, mc.cores = 4)


new_jac <- do.call(rbind, jac_mats)





(leading_eigen <- max(Re(eigen(jac)[["values"]])))

(new_leading_eigen <- max(Re(eigen(new_jac)[["values"]])))






pert <- perturb(trait_to, max_t = trait_to$call[["max_t"]], save_every = 0)


identical(pert$start$N, pert$end$N)
all.equal(pert$start$N, pert$end$N)

identical(pert$start$V, pert$end$V)
all.equal(pert$start$V, pert$end$V)
pert$start$V - pert$end$V
pert$start$V; pert$end$V

pV <- pert$end$V
pN <- pert$end$N
pjac <- sauron:::jacobian_cpp(lapply(split(pV, 1:n), rbind),
                      N = pN,
                      f = 0.1,
                      a0 = 0.5,
                      D = D,
                      C = C,
                      add_var = rep(0.1, n))

(p_leading_eigen <- max(Re(eigen(pjac)[["values"]])))
