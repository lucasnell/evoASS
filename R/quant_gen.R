
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
    args <- list(n_cores = 4, q = q)
    # the non-additive effects of traits on `r`:
    args$eta <- 0.01 * sign(eta_sign)
    # changes how the focal line's traits affect other lines' effects of competition:
    args$d <- 0.1 * sign(d_sign)
    # High d takes a very long time to run, plus it has to run longer to reach equilibrium
    if (args$d > 0) {
        args$d <- 1e-4
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
#' @param n_cores Number of cores to use. Defaults to 1.
#' @inheritParams adapt_dyn
#'
#' @return A `quant_gen` object with `nv` (for N and V output) and `fs` (for fitness
#'     and selection output) fields.
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr as_data_frame
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
                      add_var = rep(0.5, n), perturb_sd = 1,
                      n_reps = 100, start_t = 0, max_t = 1e6L,
                      min_N = 1e-4, save_every = 1e4L,
                      show_progress = TRUE, n_cores = 1) {

    stopifnot(sapply(V0, inherits, what = c("numeric", "matrix", "array")))
    stopifnot(N0 >= 0)
    stopifnot(n >= 1 && q >= 1)
    stopifnot(sapply(list(n_reps, start_t, max_t, save_every, n_cores, N0), is.numeric))
    stopifnot(sapply(list(n_reps, start_t, max_t, save_every, n_cores), length) == 1)
    stopifnot(c(n_reps, max_t, n_cores) >= 1)
    stopifnot(c(start_t, save_every, add_var, perturb_sd, min_N) >= 0)

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
                        eta = eta,
                        r0 = r0,
                        d = d,
                        add_var = add_var,
                        perturb_sd = perturb_sd,
                        start_t = start_t,
                        max_t = max_t,
                        min_N = min_N,
                        save_every = save_every,
                        show_progress = show_progress,
                        n_cores = n_cores)


    if (save_every > 0) {
        NV_ <- as_data_frame(qg$NV) %>%
            set_names(c("rep", "time", "spp", "N",
                        paste0("trait_", 1:q))) %>%
            mutate_at(vars(rep, time, spp), ~ as.integer(.x + 1)) %>%
            gather("trait", "value", starts_with("trait_"), factor_key = TRUE) %>%
            mutate(trait = as.integer(gsub("trait_", "", trait))) %>%
            mutate(rep = factor(rep, levels = 1:n_reps),
                   spp = factor(spp, levels = 1:n),
                   trait = factor(trait, levels = 1:q)) %>%
            arrange(rep, time, spp, trait) %>%
            mutate(value = ifelse(is.nan(value), NA, value))
    } else {
        NV_ <- as_data_frame(qg$NV) %>%
            set_names(c("rep", "spp", "N",
                        paste0("trait_", 1:q))) %>%
            mutate_at(vars(rep, spp), ~ as.integer(.x + 1)) %>%
            gather("trait", "value", starts_with("trait_")) %>%
            mutate(trait = as.integer(gsub("trait_", "", trait))) %>%
            mutate(rep = factor(rep, levels = 1:n_reps),
                   spp = factor(spp, levels = 1:n),
                   trait = factor(trait, levels = 1:q)) %>%
            arrange(rep, spp, trait) %>%
            mutate(value = ifelse(is.nan(value), NA, value))
    }

    FS_ <- as_data_frame(qg$FS) %>%
        set_names(c("fit", "sel")) %>%
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



