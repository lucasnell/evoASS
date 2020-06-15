
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
#' @param add_var Vector of additive genetic variances for all starting species.
#' @param start_t Number of starting iterations (before the perturbation).
#' @param n_threads Number of cores to use. Defaults to 1.
#' @inheritParams adapt_dyn
#'
#' @return A `quant_gen` object with `nv` (for N and V output) and `fs` (for fitness
#'     and selection output) fields.
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @importFrom purrr set_names
#' @importFrom dplyr mutate_at
#' @importFrom tidyr gather
#' @importFrom dplyr starts_with
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom dplyr vars
#'
quant_gen <- function(eta, d, q,
                      n = 100,
                      V0 = rep(list(matrix(0, 1, q)), n),
                      N0 = rep(1, n),
                      f = 0.1, a0 = 0.5, r0 = 0.5,
                      add_var = rep(0.1, n), perturb_sd = 1,
                      n_reps = 100, start_t = 0, max_t = 1e6L,
                      min_N = 1e-4, save_every = 1e4L,
                      show_progress = TRUE, n_threads = 1) {

    stopifnot(inherits(V0, "list"))
    stopifnot(sapply(V0, inherits, what = c("numeric", "matrix", "array")))
    stopifnot(sapply(list(eta, d, q, n, f, a0, r0, n_reps, start_t, max_t, save_every,
                          n_threads, N0), is.numeric))
    stopifnot(sapply(list(q, n, f, a0, r0, n_reps, start_t, max_t, save_every,
                          n_threads), length) == 1)
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
                        start_t = start_t,
                        max_t = max_t,
                        min_N = min_N,
                        save_every = save_every,
                        show_progress = show_progress,
                        n_threads = n_threads)


    if (save_every > 0) {
        colnames(qg$NV) <- c("rep", "time", "spp", "N", paste0("trait_", 1:q))
        NV_ <- as_tibble(qg$NV) %>%
            mutate_at(vars(rep, time, spp), as.integer) %>%
            gather("trait", "value", starts_with("trait_"), factor_key = TRUE) %>%
            mutate(trait = as.integer(gsub("trait_", "", trait))) %>%
            mutate(rep = factor(rep, levels = 1:n_reps),
                   spp = factor(spp, levels = 1:n),
                   trait = factor(trait, levels = 1:q)) %>%
            arrange(rep, time, spp, trait) %>%
            mutate(value = ifelse(is.nan(value), NA, value))
    } else {
        colnames(qg$NV) <- c("rep", "spp", "N", paste0("trait_", 1:q))
        NV_ <- as_tibble(qg$NV) %>%
            mutate_at(vars(rep, spp), as.integer) %>%
            gather("trait", "value", starts_with("trait_")) %>%
            mutate(trait = as.integer(gsub("trait_", "", trait))) %>%
            mutate(rep = factor(rep, levels = 1:n_reps),
                   spp = factor(spp, levels = 1:n),
                   trait = factor(trait, levels = 1:q)) %>%
            arrange(rep, spp, trait) %>%
            mutate(value = ifelse(is.nan(value), NA, value))
    }

    colnames(qg$FS) <- c("fit", "sel")
    FS_ <- as_tibble(qg$FS) %>%
        mutate(rep = 1L:(dplyr::n())) %>%
        dplyr::select(rep, fit, sel)


    qg_obj <- list(nv = NV_, fs = FS_, call = call_)

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
            summarize(n_ = n()) %>%
            ungroup() %>%
            .[["n_"]] %>%
            unique()
    } else {
        unq_nspp <- x$nv %>%
            filter(trait == levels(trait)[1]) %>%
            group_by(rep) %>%
            summarize(n_ = n()) %>%
            ungroup() %>%
            .[["n_"]] %>%
            unique()
    }
    cat(blu_("* Coexistence:", any(unq_nspp > 1), "\n"))
    extinct_ <- nrow(x$nv) == 0 || sum(is.na(x$nv$value)) > 0 ||
        length(unique(x$nv$rep)) < length(levels(x$nv$rep))
    if (is.null(x$call$save_every) || x$call$save_every > 0) {
        extinct_ <- extinct_ || length(unique(filter(x$nv, time == max(time))$rep)) <
            length(levels(x$nv$rep))
    }
    cat(blu_("* Total extinction:", extinct_, "\n"))

    fs_ <- x$fs %>%
        gather("par","value", fit:sel) %>%
        group_by(par) %>%
        summarize(min = min(value), max = max(value))
    cat(blu_(sprintf("* Fitness range = %.3g, %.3g\n", fs_$min[1], fs_$max[1])))
    cat(blu_(sprintf("* Selection range = %.3g, %.3g\n", fs_$min[2], fs_$max[2])))

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



#'
#' @importFrom magrittr %>%
#'
#' @noRd
#'
one_jacobian <- function(one_rep, qg_obj) {


    if (!is.null(one_rep[["time"]])) one_rep <- dplyr::filter(one_rep, time == max(time))
    one_rep

    N <- dplyr::filter(one_rep, trait == 1)[["N"]]
    spp <- as.integer(paste(dplyr::filter(one_rep, trait == 1)[["spp"]]))
    V <- one_rep %>%
        mutate(trait = paste0("V", trait)) %>%
        spread(trait, value) %>%
        select(starts_with("V", ignore.case = FALSE)) %>%
        as.matrix() %>%
        split(1:nrow(.))


    if (is.null(qg_obj$call[["n"]])) {
        n <- eval(formals(quant_gen)[["n"]])
    } else n <- eval(qg_obj$call[["n"]])
    if (is.null(qg_obj$call[["q"]])) {
        q <- eval(formals(quant_gen)[["q"]])
    } else q <- eval(qg_obj$call[["q"]])

    if (is.null(qg_obj$call[["f"]])) {
        f <- eval(formals(quant_gen)[["f"]])
    } else f <- eval(qg_obj$call[["f"]])
    if (is.null(qg_obj$call[["a0"]])) {
        a0 <- eval(formals(quant_gen)[["a0"]])
    } else a0 <- eval(qg_obj$call[["a0"]])
    if (is.null(qg_obj$call[["d"]])) {
        d <- eval(formals(quant_gen)[["d"]])
    } else d <- eval(qg_obj$call[["d"]])
    if (is.null(qg_obj$call[["eta"]])) {
        eta <- eval(formals(quant_gen)[["eta"]])
    } else eta <- eval(qg_obj$call[["eta"]])
    if (is.null(qg_obj$call[["add_var"]])) {
        add_var <- eval(formals(quant_gen)[["add_var"]])
    } else add_var <- eval(qg_obj$call[["add_var"]])
    add_var <- add_var[spp]


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


    jac <- sauron:::jacobian_cpp(V, N, f, a0, D, C, add_var)

    return(jac)

}


#' Jacobian matrices, one per rep.
#'
#' @param qg_obj A `quant_gen` object from `quant_gen` function.
#'
#' @export
#'
#' @importFrom magrittr %>%
jacobians <- function(qg_obj) {

    if (!inherits(qg_obj, "quant_gen")) {
        stop(paste("\nArgument `qg_obj` for function `jacobians` must",
                   "be of class \"quant_gen\"\n"))
    }


    jacs <- qg_obj$nv %>%
        split(.$rep) %>%
        map(one_jacobian, qg_obj = qg_obj)

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
        spread(trait, value) %>%
        select(starts_with("V", ignore.case = TRUE)) %>%
        as.matrix() %>%
        split(1:nrow(.)) %>%
        lapply(rbind)
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
        spread(trait, value) %>%
        select(starts_with("V", ignore.case = TRUE)) %>%
        as.matrix() %>%
        split(1:nrow(.)) %>%
        lapply(rbind)


    out <- list(start = list(N = N0, V = do.call(rbind, V0)),
                end = list(N = Nt, V = do.call(rbind, Vt)),
                end_obj = qg)


    return(out)


}
