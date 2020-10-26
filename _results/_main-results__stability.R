


# load packages
suppressPackageStartupMessages({
    library(sauron)
    library(tidyverse)
    library(grid)
    library(gridExtra)
    library(scales)
    library(pbmcapply)
    library(parallel)
    library(egg)
    library(cowplot)
    library(progress)
})

options(dplyr.summarise.inform = FALSE)

# This sets plotting device on LAN computer:
if (file.exists(".Rprofile")) source(".Rprofile")


# Clean captions
cc <- function(.x) {
    .x <- gsub("\n", "", .x)
    .x <- gsub("\\s+", " ", .x)
    return(.x)
}

# whether to re-do simulations (use rds files otherwise)
.REDO_SIMS <- TRUE
# whether to re-save plots
.RESAVE_PLOTS <- FALSE
# number of threads to use for simulations
.N_THREADS <- 3



# Where to save rds files
rds <- function(.x) {
    sprintf("~/GitHub/Wisconsin/sauron/_results/results-rds/%s.rds",
            gsub("\\.rds$", "", gsub("^\\/", "", .x)))
}

print_big_nums <- function(x, ...) {
    stopifnot(inherits(x, "tbl_df") || is.numeric(x))
    if (inherits(x, "tbl_df")) {
        .opt <- getOption("pillar.sigfig")
        options(pillar.sigfig = 10)
        print(x, ...)
        options(pillar.sigfig = .opt)
    } else print(x, digits = 10, ...)
    invisible(NULL)
}



save_plot <- function(plot_obj, .width, .height,
                      .name = NULL, .prefix = NULL, .suffix = NULL) {
    if (is.null(.name)) {
        fn <- gsub("_p$", "", paste(substitute(plot_obj)))
        if (!is.null(.prefix)) fn <- paste0(.prefix, fn)
        if (!is.null(.suffix)) fn <- paste0(fn, .suffix)
    } else fn <- gsub(".pdf$", "", .name)
    fn <- sprintf("_results/results-plots/%s.pdf", fn)
    message(fn)
    plot_fun <- ifelse(inherits(plot_obj, "egg"), print, plot)
    cairo_pdf(fn, width = .width, height = .height)
    plot_fun(plot_obj)
    dev.off()
    invisible(NULL)
}

unq_spp_filter <- function(..., .prec = 0.001) {
    V <- list(...)
    V <- do.call(cbind, V)
    V <- split(V, row(V))
    return(as.logical(sauron:::unq_spp_cpp(V, precision = .prec)))
}

group_spp <- function(..., .prec = 0.001) {
    V <- list(...)
    V <- do.call(cbind, V)
    V <- split(V, row(V))
    return(sauron:::group_spp_cpp(V, precision = .prec))
}

set.seed(1898348146)
etas <- map(1:6, ~ with(list(.q = .x), runif(.q * ((.q - 1) / 2), 0.1, 0.4)))
# With just one eta, it can be a simple number:
etas[[2]] <- 0.6





# =============================================================================*
# =============================================================================*

# 1. 2-trait outcomes ----

# =============================================================================*
# =============================================================================*


eta_sims <- readRDS(rds("eta_sims"))

eta_sim_df <- map_dfr(eta_sims, ~.x$ts)

#'
#' Based on these eigenvalues...
#'   * When the tradeoff is additive, the state is neutrally stable
#'   * Everything else is stable
#'
#'
eta_sim_eigs <- map_dfr(1:length(eta_sims),
                        function(i) {
                            eigs <- map_dbl(eta_sims[[i]][["jacs"]],
                                            function(.x){
                                                if (any(is.na(.x))) return(NA)
                                                z <- eigen(.x,
                                                           only.values = TRUE)
                                                z <- z$values
                                                maxIm <- max(abs(Im(z)))
                                                if (maxIm > 1e-2) return(NA_real_)
                                                max(Re(z))
                                            })
                            eta_sims[[i]]$ts %>%
                                distinct(eta1) %>%
                                mutate_all(sign) %>%
                                rename_all(~ gsub("^eta", "sign", .x)) %>%
                                mutate(e_min = min(eigs),
                                       e_max = max(eigs))
                        })


print_big_nums(eta_sim_eigs)






# =============================================================================*
# =============================================================================*

# 2. Conditional coexistence ----

# =============================================================================*
# =============================================================================*


cond_coexist <- readRDS(rds("cond_coexist"))



#'
#' They're all stable except for  sub-additive and conflicting trait 1,
#' which is neutrally stable.
#'
crossing(.V0 = c("restricted", "unrestricted"),
         .eta_sign = c(-1,1),
         .d1 = c(-0.1, 0.1)) %>%
    filter(!(.V0 == "restricted" & .eta_sign < 0)) %>%
    mutate(eigen = cond_coexist %>%
               map(~ .x$jacs[[1]]) %>%
               map(~ eigen(.x, only.values = TRUE)[["values"]]) %>%
               map_dbl(~ if (is.complex(.x)) { NA_real_ } else { max(.x) })) %>%
    print_big_nums()










# =============================================================================*
# =============================================================================*

# S1. Coexistence ----

# =============================================================================*
# =============================================================================*


#'
#' Code below shows that...
#' * Most are stable, except for the neutrally stable equilibria when
#'   `d1 + d2` is close to zero.
#' * Some reps had complex eigenvalues, none with complex
#'   leading eigenvalues when very small imaginary components weren't allowed
#'   (anything less than 1e-10).
#'   A rep with a complex 2nd eigenvalue was plotted and didn't appear to
#'   have anything resembling fluctuations.
#'
#'
#'
#'
# Takes ~1 min
giant_sim_jacs <- map_dfr(0:11,
                          function(i) {
                              fn <- rds(sprintf("giant_sims/giant_sims_%i", i))
                              .d <- readRDS(fn)
                              map_dfr(.d,
                                      function(.x) {
                                          .x[["NV"]][1,] %>%
                                              select(d, eta, add_var, sigma_N,
                                                     sigma_V, vary_d2) %>%
                                              mutate(J = list(.x[["J"]]),
                                                     seed = .x[["seed"]])
                                      })
                      })

# Makes no sense to talk about Jacobians when there's stochasticity
giant_sim_jacs <- giant_sim_jacs %>%
    filter(sigma_V == 0, sigma_N == 0) %>%
    select(-sigma_V, -sigma_N)

# Compute eigenvalues and remove Jacobians from tibble
giant_sim_jacs <- giant_sim_jacs %>%
    mutate(eigens = map(J, ~ map(.x, function(z) {
        eigen(z, only.values = TRUE)[["values"]]
    }))) %>%
    select(-J) %>%
    unnest(eigens)

# Add info about eigenvalues
giant_sim_jacs <- giant_sim_jacs %>%
    mutate(max_eigen = map_dbl(eigens, ~ max(Re(.x))),
           which_complex = map_int(eigens, function(x) {
               maxIm <- max(abs(Im(x)))
               if (maxIm > 1e-10) {
                   return(min(which(abs(Im(x)) == maxIm)))
               } else return(NA_integer_)
           }))


#'
#' Only a few reps have eigenvalues > 1, and these are very close to 1.
#' Thus I'm saying these are neutrally stable.
#' (Everything else is stable.)
#'
giant_sim_jacs %>%
    group_by(d, eta, add_var, vary_d2) %>%
    summarize(max_eigen = max(max_eigen)) %>%
    ungroup() %>%
    filter(max_eigen > 1) %>%
    mutate(above1 = max_eigen - 1)


#'
#' When `d = 1.5`, `eta = 0.6`, `add_var = 0.05`, and `vary_d2 = TRUE`,
#' 3 reps return a complex 2nd eigenvalue.
#' I'm going to plot these below.
#'
#'
giant_sim_jacs %>%
    filter(!is.na(which_complex)) %>%
    filter(which_complex == min(which_complex)) %>%
    mutate(max_im = map_dbl(eigens, ~ max(abs(Im(.x)))))




.d = 1.5
.eta = 0.6
.add_var = 0.05
.sigma_N = 0
.sigma_V = 0
.vary_d2 = TRUE
.seed <- 858577784


.n <- 100

if (.vary_d2) {
    .ds <- rep(.d, 2)
} else .ds <- c(.d, 0.1)

set.seed(.seed)

Z <- quant_gen(q = 2, eta = .eta, d = .ds, n = .n,
               add_var = rep(.add_var, .n),
               spp_gap_t = 500L, final_t = 50e3L, n_reps = 12,
               sigma_N = .sigma_N, sigma_V = .sigma_V,
               save_every = 1000L,
               n_threads = .N_THREADS, show_progress = FALSE)

if (.sigma_N == 0 && .sigma_V == 0) {
    Z$call[["eta"]] <- eval(.eta)
    Z$call[["d"]] <- eval(.ds)
    Z$call[["n"]] <- eval(.n)
    Z$call[["add_var"]] <- eval(rep(.add_var, .n))
    jacs <- jacobians(Z)
} else jacs <- NULL


map_int(jacs, function(x) {
    z <- abs(Im(eigen(x, only.values = TRUE)[["values"]]))
    min(which(z == max(z)))
})



Z$nv %>%
    filter(rep == 12, trait == 1) %>%
    # filter(time > 90e3) %>%
    ggplot(aes(time, N)) +
    geom_line(aes(color = spp), na.rm = TRUE) +
    facet_wrap(~ rep, nrow = 3) +
    scale_color_viridis_d(begin = 0.1, end = 0.9, guide = FALSE)




