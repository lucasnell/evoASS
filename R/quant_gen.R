
#
# Checks inputs to `quant_gen`, creates the D and C matrices, and edits
# the `n_threads` value if necessary.
#
#
check_quant_gen_args <- function(eta, d, q, n, V0, N0, f, a0, r0, add_var,
                                 sigma_V0, sigma_N, sigma_V, phenos, n_reps,
                                 spp_gap_t, final_t, min_N,
                                 save_every, show_progress, n_threads) {


    stopifnot(sapply(list(eta, d, q, n, f, a0, r0, sigma_V0, sigma_N, sigma_V,
                          n_reps, spp_gap_t, final_t, min_N, save_every,
                          n_threads, N0), is.numeric))
    stopifnot(sapply(list(q, n, f, a0, r0, sigma_V0, sigma_N,
                          n_reps, spp_gap_t, final_t,
                          save_every, n_threads),
                     length) == 1)
    stopifnot(length(sigma_V) %in% c(1, q))
    stopifnot(length(phenos) == 1 && is.logical(phenos))
    if (!is.null(V0)) {
        stopifnot(is.numeric(V0))
        stopifnot(inherits(V0, "matrix") || length(V0) %in% c(1, q))
        stopifnot(all(V0 >= 0))
        if (inherits(V0, "matrix")) stopifnot(nrow(V0) == q && ncol(V0) == n)
    }

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
get_quant_gen_output <- function(qg, call_, save_every, q, n, sigma_V, phenos) {

    type_fmt <- "([[:alnum:]]+)_([[:digit:]]+)"

    if (save_every > 0) {
        colnames(qg) <- c("rep", "time", "spp", "N",
                          paste0("geno_", 1:q), paste0("pheno_", 1:q))
        qg <- qg %>%
            as_tibble() %>%
            gather(key, value, starts_with("geno_"), starts_with("pheno_")) %>%
            extract(key, c("type", "axis"), type_fmt) %>%
            spread(type, value) %>%
            mutate(across(c(rep, time, spp, axis), as.integer)) %>%
            mutate(rep = factor(rep, levels = 1:max(rep)),
                   spp = factor(spp, levels = 1:n),
                   axis = factor(axis, levels = 1:q)) %>%
            select(rep, time, spp, axis, everything()) %>%
            arrange(rep, time, spp, axis) %>%
            mutate(across(c(geno, pheno), ~ ifelse(is.nan(.x), NA_real_, .x)))
    } else {
        colnames(qg) <- c("rep", "spp", "N", paste0("geno_", 1:q),
                          paste0("pheno_", 1:q))
        qg <- qg %>%
            as_tibble() %>%
            gather(key, value, starts_with("geno_"), starts_with("pheno_")) %>%
            extract(key, c("type", "axis"), type_fmt) %>%
            spread(type, value) %>%
            mutate(across(c(rep, spp, axis), as.integer)) %>%
            mutate(rep = factor(rep, levels = 1:max(rep)),
                   spp = factor(spp, levels = 1:n),
                   axis = factor(axis, levels = 1:q)) %>%
            select(rep, spp, axis, everything()) %>%
            arrange(rep, spp, axis) %>%
            mutate(across(c(geno, pheno), ~ ifelse(is.nan(.x), NA_real_, .x)))
    }
    if (all(sigma_V <= 0) || !phenos) {
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
#'     starting axis values can be derived.
#'     Set to 0 for species to start with the exact values of axes
#'     specified in `V0`.
#' @param sigma_N Standard deviation for stochasticity in population dynamics.
#' @param sigma_V Standard deviation for stochasticity in axis evolution.
#' @param phenos Logical for whether axis evolution stochasticity produces
#'     phenotypes or whether it just adds variance to the genotype directly.
#'     Defaults to `TRUE`.
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
                      V0 = 1,
                      N0 = rep(1, n),
                      f = 0.1,
                      a0 = 1e-4,
                      r0 = 0.5,
                      add_var = rep(0.01, n),
                      sigma_V0 = 1,
                      sigma_N = 0,
                      sigma_V = 0,
                      phenos = TRUE,
                      n_reps = 10,
                      spp_gap_t = 500L,
                      final_t = 5e3L,
                      min_N = 1,
                      save_every = 10L,
                      show_progress = TRUE,
                      n_threads = 1) {

    call_ <- match.call()
    # So it doesn't show the whole function if using do.call:
    if (call_[1] != as.call(quote(quant_gen()))) {
        call_[1] <- as.call(quote(quant_gen()))
    }


    args <- check_quant_gen_args(eta, d, q, n, V0, N0, f, a0, r0, add_var,
                                 sigma_V0, sigma_N, sigma_V, phenos, n_reps,
                                 spp_gap_t, final_t, min_N,
                                 save_every, show_progress, n_threads)

    C <- args$C
    D <- args$D
    n_threads <- args$n_threads

    if (length(sigma_V) == 1) sigma_V <- rep(sigma_V, q)


    if (is.null(V0)) {
        # Otherwise start at zero:
        if (sigma_V0 == 0 && all(sigma_V) == 0) {
            warning(paste("\nSimulations start with species having all axes",
                          "at zero, and you aren't providing stochasticity in",
                          "either the starting values or via phenotypes.",
                          "Because this is often an unstable equilibrium,",
                          "these simulations may be odd or boring.",
                          "Continuing anyway..."))
        }
        V0 <- matrix(0, q, n)
    } else if (!inherits(V0, "matrix")) {
        stopifnot(length(V0) == 1 || length(V0) == q)
        V0 <- matrix(V0, q, n)
    }

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
                        phenos = phenos,
                        spp_gap_t = spp_gap_t,
                        final_t = final_t,
                        min_N = min_N,
                        save_every = save_every,
                        show_progress = show_progress,
                        n_threads = n_threads)


    qg_obj <- get_quant_gen_output(qg, call_, save_every, q, n, sigma_V, phenos)

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
            filter(time == max(time), axis == levels(axis)[1]) %>%
            group_by(rep) %>%
            summarize(n_ = dplyr::n()) %>%
            ungroup() %>%
            .[["n_"]] %>%
            unique()
    } else {
        unq_nspp <- x$nv %>%
            filter(axis == levels(axis)[1]) %>%
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

    N <- dplyr::filter(one_rep, axis == 1)[["N"]]
    spp <- as.integer(paste(dplyr::filter(one_rep, axis == 1)[["spp"]]))
    V <- one_rep %>%
        dplyr::mutate(axis = paste0("V", axis)) %>%
        tidyr::spread(axis, geno) %>%
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




