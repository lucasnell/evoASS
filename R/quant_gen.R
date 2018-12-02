

#' Quantitative genetics.
#'
#' @param n_reps Number of reps to perform.
#' @param add_var Vector of additive genetic variances for all starting species.
#' @param delta Standard deviation of perturbation after starting iterations.
#' @param start_t Number of starting iterations.
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
#'
quant_gen <- function(n_reps, V0, N0, f, g, eta, r0, d, add_var, delta, start_t,
                      max_t, min_N, save_every,
                      show_progress = TRUE, n_cores = 1) {

    n <- length(N0)
    q <- length(V0[[1]])

    stopifnot(n >= 2 && q >= 2)
    stopifnot(sapply(list(n_reps, start_t, max_t, save_every, n_cores, N0), is.numeric))
    stopifnot(sapply(list(n_reps, start_t, max_t, save_every, n_cores), length) == 1)
    stopifnot(c(n_reps, max_t, n_cores) >= 1)
    stopifnot(c(start_t, save_every, add_var, delta, min_N, N0) >= 0)

    call_ <- match.call()
    # So it doesn't show the whole function if using do.call:
    if (call_[1] != as.call(quote(quant_gen()))) {
        call_[1] <- as.call(quote(quant_gen()))
    }

    qg <- quant_gen_cpp(n_reps = n_reps,
                        V0 = V0,
                        N0 = N0,
                        f = f,
                        g = g,
                        eta = eta,
                        r0 = r0,
                        d = d,
                        add_var = add_var,
                        delta = delta,
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
            mutate_at(vars(rep, spp, trait), factor) %>%
            arrange(rep, time, spp, trait) %>%
            identity()
    } else {
        NV_ <- as_data_frame(qg$NV) %>%
            set_names(c("rep", "spp", "N",
                        paste0("trait_", 1:q))) %>%
            mutate_at(vars(rep, spp), ~ as.integer(.x + 1)) %>%
            gather("trait", "value", starts_with("trait_")) %>%
            mutate(trait = as.integer(gsub("trait_", "", trait))) %>%
            mutate_at(vars(rep, spp, trait), factor) %>%
            arrange(rep, spp, trait) %>%
            identity()
    }

    FS_ <- as_data_frame(qg$FS) %>%
        set_names(c("fit", "sel")) %>%
        mutate(rep = 1L:n()) %>%
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
print.quant_gen <- function(x, digits = max(3, getOption("digits") - 3), ...) {

    cat(crayon::cyan$bold("# Output from quantitative genetics\n"))
    cat(crayon::green("< N and V >\n"))
    print(x$nv, digits, ...)
    cat(crayon::green("< Fitness and selection >\n"))
    print(x$fs, digits, ...)

    invisible(x)

}



