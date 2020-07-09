
#' Create list of arguments for running a set of simulations using `quant_gen`.
#'
#' @param eta_sign Sign of the `eta` parameter.
#' @param d_sign Sign of the `d` parameter.
#' @param q Number of traits.
#'
#' @export
#'
quant_gen_args <- function(eta_sign, d_sign, q) {
    stopifnot(is.numeric(eta_sign), is.numeric(d_sign), is.numeric(q))
    stopifnot(length(eta_sign) == 1, length(d_sign) == 1, length(q) == 1)
    # List for all parameters:
    args <- list(n_threads = 4, q = q)
    # the non-additive effects of traits on `r`:
    args$eta <- 0.6 * sign(eta_sign)
    # changes how the focal line's traits affect other lines' effects of competition:
    args$d <- 1e-4 * sign(d_sign)
    # d > 0 has to run longer to reach equilibrium
    if (args$d > 0) {
        args$max_t <- 2e7L
        args$save_every <- 1e5L
    }
    return(args)
}



#' Quantitative genetics.
#'
#' @param n_reps Number of reps to perform.
#' @param perturb_sd Standard deviation of the perturbation.
#' @param sigma_N Standard deviation for stochasticity in population dynamics.
#' @param sigma_V Standard deviation for stochasticity in trait evolution.
#' @param add_var Vector of additive genetic variances for all starting species.
#' @param start_t Number of starting iterations (before the perturbation).
#' @param n_threads Number of cores to use. Defaults to 1.
#' @inheritParams adapt_dyn
#'
#' @return A `quant_gen` object with `nv` (for N and V output) and
#'     `call` (for original call) fields.
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom purrr set_names
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom dplyr starts_with
#' @importFrom dplyr arrange
#' @importFrom dplyr select
#' @importFrom tidyr gather
#' @importFrom tidyr spread
#' @importFrom tidyr extract
#'
#'
#'
quant_gen <- function(eta, d, q,
                      n = 100,
                      V0 = rep(list(matrix(0, q, 1)), n),
                      N0 = rep(1, n),
                      f = 0.1,
                      a0 = 0.5,
                      r0 = 0.5,
                      add_var = rep(0.01, n),
                      perturb_sd = 1,
                      sigma_N = 0,
                      sigma_V = 0,
                      n_reps = 100,
                      start_t = 0,
                      max_t = 1e6L,
                      min_N = 1e-4,
                      save_every = 1e4L,
                      show_progress = TRUE,
                      n_threads = 1) {

    stopifnot(inherits(V0, "list"))
    stopifnot(sapply(V0, inherits, what = c("numeric", "matrix", "array")))
    stopifnot(sapply(list(eta, d, q, n, f, a0, r0, perturb_sd, sigma_N, sigma_V,
                          n_reps, start_t, max_t, save_every,
                          n_threads, N0), is.numeric))
    stopifnot(sapply(list(q, n, f, a0, r0, perturb_sd, sigma_N, sigma_V,
                          n_reps, start_t, max_t, save_every, n_threads),
                     length) == 1)
    stopifnot(sapply(V0, function(x) all(x >= 0)))
    stopifnot(n >= 1 && q >= 1)
    stopifnot(N0 >= 0)
    stopifnot(c(n_reps, max_t, n_threads) >= 1)
    stopifnot(c(start_t, save_every, add_var, perturb_sd, min_N) >= 0)

    stopifnot(length(eta) %in% c(1, q^2))
    stopifnot(length(d) %in% c(1, q))

    if (n_threads > 1 && !using_openmp()) {
        message("\nOpenMP not enabled. Only 1 thread will be used.\n")
        n_threads <- 1
    }


    C <- matrix(eta[1], q, q)
    if (length(eta) == q^2) {
        stopifnot(inherits(eta, "matrix") && identical(dim(eta), as.integer(c(q,q))) &&
                      isSymmetric(eta))
        C <- eta
    }
    diag(C) <- 1

    D <- matrix(0, q, q)
    if (length(d) == 1) {
        diag(D) <- rep(d, q)
    } else if (length(d) == q) {
        diag(D) <- d
    }

    call_ <- match.call()
    # So it doesn't show the whole function if using do.call:
    if (call_[1] != as.call(quote(quant_gen()))) {
        call_[1] <- as.call(quote(quant_gen()))
    }

    qg <- quant_gen_cpp(n_reps = n_reps,
                        V0 = V0,
                        N0 = N0,
                        f = f,
                        a0 = a0,
                        C = C,
                        r0 = r0,
                        D = D,
                        add_var = add_var,
                        perturb_sd = perturb_sd,
                        sigma_N = sigma_N,
                        sigma_V = sigma_V,
                        start_t = start_t,
                        max_t = max_t,
                        min_N = min_N,
                        save_every = save_every,
                        show_progress = show_progress,
                        n_threads = n_threads)


    type_fmt <- "([[:alnum:]]+)_([[:digit:]]+)"

    if (save_every > 0) {
        colnames(qg) <- c("rep", "time", "spp", "N",
                          paste0("geno_", 1:q), paste0("pheno_", 1:q))
        qg <- qg %>%
            as_tibble() %>%
            gather(key, value, starts_with("geno_"), starts_with("pheno_")) %>%
            extract(key, c("type", "trait"), type_fmt) %>%
            spread(type, value) %>%
            mutate(across(c(rep, time, spp, trait),
                          # function(x) as.integer(x) %>%
                              factor(levels = 1:max(.)))) %>%
            select(rep, time, spp, trait, everything()) %>%
            arrange(rep, time, spp, trait) %>%
            mutate(across(c(geno, pheno), ~ ifelse(is.nan(.x), NA_real_, .x)))
    } else {
        colnames(qg) <- c("rep", "spp", "N", paste0("geno_", 1:q),
                          paste0("pheno_", 1:q))
        qg <- qg %>%
            as_tibble() %>%
            gather(key, value, starts_with("geno_"), starts_with("pheno_")) %>%
            extract(key, c("type", "trait"), type_fmt) %>%
            spread(type, value) %>%
            mutate(across(c(rep, spp, trait),
                          function(x) as.integer(x) %>%
                              factor(levels = 1:max(.)))) %>%
            select(rep, spp, trait, everything()) %>%
            arrange(rep, spp, trait) %>%
            mutate(across(c(geno, pheno), ~ ifelse(is.nan(.x), NA_real_, .x)))
    }
    if (sigma_V <= 0) {
        qg <- select(qg, -pheno)
    }


    qg_obj <- list(nv = qg, call = call_)

    class(qg_obj) <- "quant_gen"

    return(qg_obj)
}




#' Print a `quant_gen` object.
#'
#' @param x an object of class \code{quant_gen}.
#' @param digits the number of digits to be printed.
#' @param ... arguments passed to and from other methods.
#'
#' @export
#' @noRd
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr ungroup
#'
print.quant_gen <- function(x, digits = max(3, getOption("digits") - 3), ...) {

    blu_ <- crayon::make_style("dodgerblue")

    cat(crayon::inverse$bold(" -- Output from quant_gen -- \n"))
    if (is.null(x$call$save_every) || x$call$save_every > 0) {
        unq_nspp <- x$nv %>%
            filter(time == max(time), trait == levels(trait)[1]) %>%
            group_by(rep) %>%
            summarize(n_ = dplyr::n()) %>%
            ungroup() %>%
            .[["n_"]] %>%
            unique()
    } else {
        unq_nspp <- x$nv %>%
            filter(trait == levels(trait)[1]) %>%
            group_by(rep) %>%
            summarize(n_ = dplyr::n()) %>%
            ungroup() %>%
            .[["n_"]] %>%
            unique()
    }
    cat(blu_("* Coexistence:", any(unq_nspp > 1), "\n"))
    extinct_ <- nrow(x$nv) == 0 || sum(is.na(x$nv$geno)) > 0 ||
        length(unique(x$nv$rep)) < length(levels(x$nv$rep))
    if (is.null(x$call$save_every) || x$call$save_every > 0) {
        extinct_ <- extinct_ || length(unique(filter(x$nv,
                                                     time == max(time))$rep)) <
            length(levels(x$nv$rep))
    }
    cat(blu_("* Total extinction:", extinct_, "\n"))

    cat("\n\n")

    print(x$nv, n = 10)

    invisible(x)

}



#' Stable points from quantitative genetics simulations and analytical solutions.
#'
#' *Note:* Only works for 2 traits for now.
#'
#' @inheritParams adapt_dyn
#' @param return_geom Logical for whether to return a `ggplot2` `geom_*` object for
#'     plotting rather than a `tibble`. Defaults to `FALSE`.
#' @param line_n Number of points to use for lines. Defaults to `1000`.
#' @param ... Other arguments to pass to the call to `geom_point` or `geom_path`.
#'
#' @return A `tibble` if `return_geom` is `FALSE`, and a  `geom_*` object if it's `TRUE`.
#'
#' @export
#'
#' @importFrom tibble tibble
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 aes
#'
stable_points <- function(eta, f = 0.1, a0 = 0.5, r0 = 0.5,
                                 return_geom = FALSE, line_n = 1000, ...) {
    if (eta < 0) {
        xy <- sqrt(0.5 * ((r0 / (f * (1 + eta))) - 1))
        pts <- tibble(V1 = xy, V2 = xy)
        if (!return_geom) return(pts)
        geom <- geom_point(data = pts, aes(V1, V2), ...)
    } else if (eta > 0) {
        xy <- c(0, 1) * sqrt((r0 / f) - 1)
        pts <- tibble(V1 = xy, V2 = rev(xy))
        if (!return_geom) return(pts)
        geom <- geom_point(data = pts, aes(V1, V2), ...)
    } else {
        radius <- sqrt((r0 / f) - 1)
        pts <- tibble(V1 = seq(0, radius, length.out = line_n),
                      V2 = sqrt(radius ^2 - V1^2))
        if (!return_geom) return(pts)
        geom <- geom_path(data = pts, aes(V1, V2), ...)
    }
    return(geom)
}


#' @describeIn stable_points Unstable points
#'
#' @inheritParams stable_points
#'
#' @export
#'
#' @importFrom ggplot2 geom_blank
#'
unstable_points <- function(eta, f = 0.1, a0 = 0.5, r0 = 0.5, return_geom = FALSE, ...) {
    if (eta < 0) {
        xy <- c(0, 1) * sqrt((r0 / f) - 1)
        pts <- tibble(V1 = xy, V2 = rev(xy))
        if (!return_geom) return(pts)
        geom <- geom_point(data = pts, aes(V1, V2), ...)
    } else if (eta > 0) {
        xy <- sqrt(0.5 * ((r0 / (f * (1 + eta))) - 1))
        pts <- tibble(V1 = xy, V2 = xy)
        if (!return_geom) return(pts)
        geom <- geom_point(data = pts, aes(V1, V2), ...)
    } else {
        if (!return_geom) return(tibble(V1 = numeric(0), V2 = numeric(0)))
        geom <- geom_blank()
    }
    return(geom)
}




#' Population abundances from quantitative genetics analytical solutions.
#'
#' *Note:* Only works for 2 traits for now.
#'
#' @inheritParams adapt_dyn
#'
#' @return A numeric vector of abundances, of length `n`.
#'
#' @export
#'
pop_sizes <- function(n, eta, d, ...) {

    stopifnot(is.numeric(n) && length(eta) == 1)
    stopifnot(is.numeric(eta))
    stopifnot(is.numeric(d))

    if (length(d) > 2 || length(eta) > 1) {
        stop("`pop_sizes` is only programmed for 2 traits")
    }
    if (length(d) == 1) d <- rep(d, 2)

    LIST <- modifyList(formals(quant_gen), list(...))
    LIST <- modifyList(LIST,
                       list(eta = eta, n = n, d1 = d[1], d2 = d[2]))

    N <- with(LIST,
              {
                  num <- f * (1 + eta) * exp((r0 / (f * (1 + eta))) - 1)
                  denom <- a0 * (
                      1 + exp(-0.5 * (r0 / (f * (1 + eta)) - 1) * (d1 + d2)) *
                          (n - 1)
                  )
                  (num / denom)
              })

    rep(N, n)

}









#'
#' @importFrom magrittr %>%
#' @importFrom tidyr spread
#'
#' @noRd
#'
one_jacobian <- function(one_rep, qg_obj, evo_only) {

    if (!is.null(one_rep[["time"]])) {
        one_rep <- dplyr::filter(one_rep, time == max(time))
    }

    N <- dplyr::filter(one_rep, trait == 1)[["N"]]
    spp <- as.integer(paste(dplyr::filter(one_rep, trait == 1)[["spp"]]))
    V <- one_rep %>%
        dplyr::mutate(trait = paste0("V", trait)) %>%
        tidyr::spread(trait, geno) %>%
        dplyr::select(starts_with("V", ignore.case = FALSE)) %>%
        as.matrix() %>%
        t()


    if (is.null(qg_obj$call[["n"]])) {
        n <- eval(formals(quant_gen)[["n"]])
    } else n <- eval(qg_obj$call[["n"]], parent.frame(2L))
    if (is.null(qg_obj$call[["q"]])) {
        stop("arg `q` is NULL in quant_gen object call")# should never be NULL
    } else q <- eval(qg_obj$call[["q"]], parent.frame(2L))

    if (is.null(qg_obj$call[["f"]])) {
        f <- eval(formals(quant_gen)[["f"]])
    } else f <- eval(qg_obj$call[["f"]], parent.frame(2L))
    if (is.null(qg_obj$call[["a0"]])) {
        a0 <- eval(formals(quant_gen)[["a0"]])
    } else a0 <- eval(qg_obj$call[["a0"]], parent.frame(2L))
    if (is.null(qg_obj$call[["r0"]])) {
        r0 <- eval(formals(quant_gen)[["r0"]])
    } else r0 <- eval(qg_obj$call[["r0"]], parent.frame(2L))
    if (is.null(qg_obj$call[["d"]])) {
        stop("arg `d` is NULL in quant_gen object call")# should never be NULL
    } else d <- eval(qg_obj$call[["d"]], parent.frame(2L))
    if (is.null(qg_obj$call[["eta"]])) {
        stop("arg `eta` is NULL in quant_gen object call")# should never be NULL
    } else eta <- eval(qg_obj$call[["eta"]], parent.frame(2L))
    if (is.null(qg_obj$call[["add_var"]])) {
        add_var <- eval(formals(quant_gen)[["add_var"]])
    } else add_var <- eval(qg_obj$call[["add_var"]], parent.frame(2L))
    add_var <- add_var[spp]


    if (length(add_var) == 0) return(matrix(NA_real_, 0, 0))
    stopifnot(is.numeric(eta))

    C <- matrix(eta[1], q, q)
    if (length(eta) == q^2) {
        stopifnot(inherits(eta, "matrix") &&
                      identical(dim(eta), as.integer(c(q,q))) &&
                      isSymmetric(eta))
        C <- eta
    }
    diag(C) <- 1

    D <- matrix(0, q, q)
    if (length(d) == 1) {
        diag(D) <- rep(d, q)
    } else if (length(d) == q) {
        diag(D) <- d
    }

    jac <- jacobian_cpp(V, N, f, a0, r0, D, C, add_var, evo_only)

    # Now dealing with the step function that keeps traits >= 0
    # The ramp function is how we made traits >= 0, and the
    # Heaviside step function is the derivative of the ramp function
    heaviside <- function(x) ifelse(x > 0, 1, 0)
    if (length(add_var) == 1) {
        S <- matrix(add_var)
    } else S <- diag(add_var)
    deltaV <- sel_str_cpp(V = V, N = N, f = f, a0 = a0,
                          C = C, r0 = r0, D = D) %*% S
    newV <- as.numeric(V + deltaV)
    if (evo_only) {
        jac <- diag(heaviside(newV)) %*% jac
    } else {
        newN <- N * as.numeric(F_t_cpp(V, N, f, a0, C, r0, D))
        jac <- diag(heaviside(c(newV, newN))) %*% jac
    }



    return(jac)

}


#' Jacobian matrices, one per rep.
#'
#' @param qg_obj A `quant_gen` object from `quant_gen` function.
#'
#' @export
#'
#' @importFrom magrittr %>%
#'
jacobians <- function(qg_obj, evo_only = FALSE) {

    if (!inherits(qg_obj, "quant_gen")) {
        stop(paste("\nArgument `qg_obj` for function `jacobians` must",
                   "be of class \"quant_gen\"\n"))
    }

    if (length(evo_only) != 1 || !is.logical(evo_only)) {
        stop("\n`evo_only` argument must be a single logical")
    }


    jacs <- qg_obj$nv %>%
        split(.$rep) %>%
        lapply(one_jacobian, qg_obj = qg_obj, evo_only = evo_only)

    return(jacs)

}





#' @describeIn quant_gen Add a perturbation to a `quant_gen` object.
#'
#' @param obj A `quant_gen` object.
#' @inheritParams quant_gen
#' @param new_V List of trait values for the new species being added.
#'     Defaults to `NULL`, which causes no new species to be add.
#' @param new_N Numeric vector of new species to add. Defaults to `NULL`, which
#'     causes no new species to be add.
#' @param ... Other parameters to pass to `quant_gen`.
#'
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#'
perturb.quant_gen <- function(obj,
                              max_t,
                              save_every,
                              start_t = 0,
                              perturb_sd = 0,
                              new_V = NULL,
                              new_N = NULL,
                              ...) {

    if (!is.null(new_V)) {
        stopifnot(inherits(new_V, "list"))
        stopifnot(sapply(new_V, inherits, what = c("numeric", "matrix", "array")))
        stopifnot(sapply(new_V, function(x) all(x >= 0)))
        stopifnot(sapply(new_V, ncol) == 1)
        V <- new_V
    } else V <- list()

    if (!is.null(new_N)) {
        stopifnot(is.numeric(new_N))
        stopifnot(new_N >= 0)
        N <- new_N
    } else N <- numeric(0)

    stopifnot(sapply(list(start_t, max_t, save_every, perturb_sd), is.numeric))
    stopifnot(sapply(list(start_t, max_t, save_every, perturb_sd), length) == 1)
    stopifnot(max_t >= 1)
    stopifnot(c(start_t, save_every, perturb_sd) >= 0)

    obj$nv <- filter(obj$nv, rep == rep[1])
    if ("time" %in% colnames(obj$nv)) {
        obj$nv <- filter(obj$nv, time == max(time))
    }

    N0 <- obj$nv %>%
        filter(trait == 1) %>%
        .[["N"]]
    N0 <- c(N0, N)

    V0 <- obj$nv %>%
        mutate(trait = paste0("V", trait)) %>%
        spread(trait, geno) %>%
        select(starts_with("V", ignore.case = TRUE)) %>%
        as.matrix() %>%
        split(1:nrow(.)) %>%
        lapply(cbind)
    V0 <- c(V0, V)

    args <- obj$call %>%
        as.list() %>%
        .[-1]

    args[["max_t"]] <- max_t
    args[["save_every"]] <- save_every
    args[["start_t"]] <- start_t
    args[["perturb_sd"]] <- perturb_sd
    args[["V0"]] <- V0
    args[["N0"]] <- N0
    args[["n"]] <- length(N0)

    args[["n_reps"]] <- 1
    args[["n_threads"]] <- 1

    if (!is.null(args[["add_var"]])) {
        if (length(unique(args[["add_var"]])) > 1) {
            stop(paste("perturb.quant_gen not programmed for species with",
                       "differing additive genetic variations"))
        }
        args[["add_var"]] <- rep(unique(args[["add_var"]]), length(N0))
    }

    new_args <- list(...)
    for (n in names(new_args)) args[[n]] <- new_args[[n]]

    qg <- do.call(quant_gen, args)

    NV <- qg$nv
    if ("time" %in% colnames(NV)) {
        NV <- filter(NV, time == max(time))
    }
    Nt <- NV %>%
        filter(trait == 1) %>%
        .[["N"]]
    Vt <- NV %>%
        mutate(trait = paste0("V", trait)) %>%
        spread(trait, geno) %>%
        select(starts_with("V", ignore.case = TRUE)) %>%
        as.matrix() %>%
        split(1:nrow(.)) %>%
        lapply(cbind)


    out <- list(start = list(N = N0, V = do.call(cbind, V0)),
                end = list(N = Nt, V = do.call(cbind, Vt)),
                end_obj = qg)


    return(out)


}
