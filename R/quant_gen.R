
#
# Checks inputs to `quant_gen`, creates the D and C matrices, and edits
# the `n_threads` value if necessary.
#
#
check_quant_gen_args <- function(eta, d, q, n, V0, N0, f, a0, r0, add_var,
                                 sigma_V0, sigma_N, sigma_V, n_reps,
                                 spp_gap_t, final_t, min_N, save_every,
                                 show_progress, n_threads) {


    stopifnot(sapply(list(eta, d, q, n, f, a0, r0, sigma_V0, sigma_N, sigma_V,
                          n_reps, spp_gap_t, final_t, save_every,
                          n_threads, N0), is.numeric))
    stopifnot(sapply(list(q, n, f, a0, r0, sigma_V0, sigma_N, sigma_V,
                          n_reps, spp_gap_t, final_t, save_every, n_threads),
                     length) == 1)
    stopifnot(inherits(V0, "matrix") && is.numeric(V0))
    stopifnot(all(V0 >= 0))
    stopifnot(nrow(V0) == q && ncol(V0) == n)

    stopifnot(n >= 1 && q >= 1)
    stopifnot(N0 >= 0)
    stopifnot(c(n_reps, final_t, n_threads) >= 1)
    stopifnot(c(spp_gap_t, save_every, add_var, sigma_V0, min_N) >= 0)

    stopifnot(length(eta) %in% c(1, q^2))
    stopifnot(length(d) %in% c(1, q))

    if (n_threads > 1 && !using_openmp()) {
        message("\nOpenMP not enabled. Only 1 thread will be used.\n")
        n_threads <- 1
    }

    C <- matrix(eta[1], q, q)
    if (length(eta) == q^2) {
        stopifnot(inherits(eta, "matrix") &&
                      identical(dim(eta), as.integer(c(q,q))) &&
                      isSymmetric(eta))
        C <- eta
    }
    diag(C) <- 1

    D <- matrix(0, q, q)
    diag(D) <- d

    return(list(C = C, D = D, n_threads = n_threads))


}



#
# Turns the raw output from `quant_gen_cpp` into a `quant_gen` object
#
#
get_quant_gen_output <- function(qg, call_, save_every, q, sigma_V) {

    type_fmt <- "([[:alnum:]]+)_([[:digit:]]+)"

    if (save_every > 0) {
        colnames(qg) <- c("rep", "time", "spp", "N",
                          paste0("geno_", 1:q), paste0("pheno_", 1:q))
        qg <- qg %>%
            as_tibble() %>%
            gather(key, value, starts_with("geno_"), starts_with("pheno_")) %>%
            extract(key, c("type", "trait"), type_fmt) %>%
            spread(type, value) %>%
            mutate(across(c(rep, time, spp, trait), as.integer)) %>%
            mutate(across(c(rep, spp, trait),
                          ~ factor(.x, levels = 1:max(.x)))) %>%
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


    qg_obj <- structure(list(nv = qg, call = call_),
                        class = "quant_gen")

    return(qg_obj)
}



#' Quantitative genetics.
#'
#' @param n_reps Number of reps to perform.
#' @param final_t Length of final time period where all species are together.
#' @param sigma_V0 Standard deviation for normal distribution from which
#'     starting trait values can be derived.
#'     Set to 0 for species to start with the exact values of traits
#'     specified in `V0`.
#' @param sigma_N Standard deviation for stochasticity in population dynamics.
#' @param sigma_V Standard deviation for stochasticity in trait evolution.
#' @param add_var Vector of additive genetic variances for all starting species.
#' @param spp_gap_t Time period between each species introduction.
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
                      n = 10,
                      V0 = matrix(0, q, n),
                      N0 = rep(1, n),
                      f = 0.1,
                      a0 = 0.5,
                      r0 = 0.5,
                      add_var = rep(0.01, n),
                      sigma_V0 = 1,
                      sigma_N = 0,
                      sigma_V = 0,
                      n_reps = 100,
                      spp_gap_t = 5e3L,
                      final_t = 20e3L,
                      min_N = 1e-4,
                      save_every = 10L,
                      show_progress = TRUE,
                      n_threads = 1) {

    call_ <- match.call()
    # So it doesn't show the whole function if using do.call:
    if (call_[1] != as.call(quote(quant_gen()))) {
        call_[1] <- as.call(quote(quant_gen()))
    }


    args <- check_quant_gen_args(eta, d, q, n, V0, N0, f, a0, r0, add_var,
                                 sigma_V0, sigma_N, sigma_V, n_reps,
                                 spp_gap_t, final_t, min_N, save_every,
                                 show_progress, n_threads)

    invisible(list2env(args, environment()))

    qg <- quant_gen_cpp(n_reps = n_reps,
                        V0 = split(t(V0), 1:ncol(V0)),
                        Vp0 = list(),
                        N0 = N0,
                        f = f,
                        a0 = a0,
                        C = C,
                        r0 = r0,
                        D = D,
                        add_var = add_var,
                        sigma_V0 = sigma_V0,
                        sigma_N = sigma_N,
                        sigma_V = sigma_V,
                        spp_gap_t = spp_gap_t,
                        final_t = final_t,
                        min_N = min_N,
                        save_every = save_every,
                        show_progress = show_progress,
                        n_threads = n_threads)


    qg_obj <- get_quant_gen_output(qg, call_, save_every, q, sigma_V)

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

    cat("\n")

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
#' @param which_rep Integer for which repetition to perturb.
#' @param ... Other parameters to pass to `quant_gen`.
#'
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr bind_rows
#'
perturb.quant_gen <- function(obj,
                              final_t,
                              save_every,
                              spp_gap_t = 0,
                              sigma_V0 = 0,
                              new_V = NULL,
                              new_Vp = NULL,
                              new_N = NULL,
                              which_rep = 1,
                              ...) {

    call_ <- match.call()
    # So it doesn't show the whole function if using do.call:
    if (call_[1] != as.call(quote(perturb.quant_gen()))) {
        call_[1] <- as.call(quote(perturb.quant_gen()))
    }

    if (!is.null(new_V) && is.null(new_N)) {
        stop("if providing new_V, you must also provide new_N")
    }
    if (!is.null(new_N) && is.null(new_V)) {
        stop("if providing new_N, you must also provide new_V")
    }
    if (!is.null(new_Vp) && is.null(new_V)) {
        stop("if providing new_Vp, you must also provide new_V")
    }
    if (!is.null(new_Vp) && !("pheno" %in% colnames(obj$nv))) {
        stop(paste("you can't provide new_Vp if no trait stochasticity",
                   "in original call"))
    }


    stopifnot(sapply(list(spp_gap_t, final_t, save_every, sigma_V0), is.numeric))
    stopifnot(sapply(list(spp_gap_t, final_t, save_every, sigma_V0), length) == 1)
    stopifnot(final_t >= 1)
    stopifnot(c(spp_gap_t, save_every, sigma_V0) >= 0)

    if (sigma_V0 > 0 && !is.null(new_Vp)) {
        stop("\nwhen sigma_V0 > 0, it negates the provided new_Vp")
    }


    # Get parameters based on `...`, original call, or defaults (in that order)
    # If `incl_pars = FALSE`, then it only looks at last two
    get_par <- function(pname, incl_pars = TRUE) {
        if (incl_pars) {
            if (!is.null(list(...)[[pname]])) {
                list(...)[[pname]]
            } else if (!is.null(obj$call[[pname]])) {
                eval(obj$call[[pname]], parent.frame(3L))
            } else {
                eval(formals(quant_gen)[[pname]])
            }
        } else {
            if (!is.null(obj$call[[pname]])) {
                eval(obj$call[[pname]], parent.frame(3L))
            } else {
                eval(formals(quant_gen)[[pname]])
            }
        }
    }

    if (sign(get_par("sigma_V")) != sign(get_par("sigma_V", FALSE)) ||
        sign(get_par("sigma_N")) != sign(get_par("sigma_N", FALSE))) {
        stop(paste("\nperturb.quant_gen not designed to go from",
                   "deterministic to stochastic simulations, or vice versa"))
    }

    new_spp <- !is.null(new_V)

    if (new_spp) {

        stopifnot(inherits(new_V, "matrix") && is.numeric(new_V))
        stopifnot(all(new_V >= 0))
        stopifnot(nrow(new_V) == obj$call[["q"]])
        V <- lapply(split(t(new_V), 1:ncol(new_V)), cbind)

        stopifnot(is.numeric(new_N))
        stopifnot(new_N >= 0)
        N <- new_N

        if (!is.null(new_Vp)) {

            stopifnot(inherits(new_Vp, "matrix") && is.numeric(new_Vp))
            stopifnot(all(new_Vp >= 0))
            stopifnot(identical(dim(new_Vp), dim(new_V)))
            Vp <- lapply(split(t(new_Vp), 1:ncol(new_Vp)), cbind)

        } else {

            sigma_V <- get_par("sigma_V")

            if (sigma_V > 0) {
                Vp <- lapply(V, function(x) x * exp(rnorm(length(x),
                                                          sd = sigma_V)))
            } else Vp <- V

        }


    } else {

        V <- list()
        N <- numeric(0)
        Vp <- list()

    }


    NV <- filter(obj$nv, rep == which_rep) %>% select(-rep)
    if (nrow(NV) == 0) stop("\nThe rep you passed doesn't exist.")
    if ("time" %in% colnames(NV)) {
        NV <- filter(NV, time == max(time)) %>% select(-time)
    }

    N0 <- filter(NV, trait == 1)[["N"]]
    N0 <- c(N0, N)

    V0 <- NV %>%
        select(spp, trait, geno) %>%
        mutate(trait = paste0("V", trait)) %>%
        spread(trait, geno) %>%
        select(starts_with("V", ignore.case = TRUE)) %>%
        as.matrix() %>%
        split(1:nrow(.)) %>%
        lapply(cbind)
    V0 <- c(V0, V)


    if ("pheno" %in% colnames(NV)) {
        Vp0 <- NV %>%
            select(spp, trait, pheno) %>%
            mutate(trait = paste0("V", trait)) %>%
            spread(trait, pheno) %>%
            select(starts_with("V", ignore.case = TRUE)) %>%
            as.matrix() %>%
            split(1:nrow(.)) %>%
            lapply(cbind)
        Vp0 <- c(Vp0, Vp)
    } else Vp0 <- c(V0, Vp)


    # these are useful below
    n <- length(N0)
    q <- nrow(V0[[1]])

    call_args <- obj$call %>%
        as.list() %>%
        .[-1]

    new_args <- list(...)

    # Based on `...`, original call, or defaults (in that order)
    args <- lapply(names(formals(quant_gen)), get_par)

    # based on inputs to this function
    args[["final_t"]] <- final_t
    args[["save_every"]] <- save_every
    args[["spp_gap_t"]] <- spp_gap_t
    args[["sigma_V0"]] <- sigma_V0
    args[["V0"]] <- V0
    args[["N0"]] <- N0
    args[["n_reps"]] <- 1
    args[["n_threads"]] <- 1

    # C and D
    args[["n"]] <- n  # <-- prevents error in `check_quant_gen_args`
    C_and_D <- do.call(check_quant_gen_args, args)

    args[["C"]] <- C_and_D[["C"]]
    args[["D"]] <- C_and_D[["D"]]

    args[["Vp0"]] <- Vp0  # has to be added after `check_quant_gen_args`

    # Not needed for cpp version
    for (x in c("eta", "d", "q", "n")) args[[x]] <- NULL

    new_obj <- do.call(quant_gen_cpp, args)

    new_obj <- get_quant_gen_output(new_obj, NULL, save_every, q, args$sigma_V)

    out_NV <- bind_rows(obj$nv %>%
                            mutate(period = "before"),
                        new_obj$nv %>%
                            mutate(period = "after",
                                   time = time + max(obj$nv[["time"]])))

    out_obj <- structure(list(nv = out_NV,
                              calls = list(orig = obj$call, perturb = call_)),
                     class = "perturb_qg")


    return(out)


}
