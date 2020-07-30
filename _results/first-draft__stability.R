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
                                                max(eigen(.x)$values)
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
               map_dbl(~ if (is.complex(.x)) { NA_real_ } else { max(.x) }))










# =============================================================================*
# =============================================================================*

# S1. Coexistence ----

# =============================================================================*
# =============================================================================*




vary_d1_coexist_sims <- readRDS(rds("vary_d1_coexist_sims"))

vary_d1_coexist_df <- map_dfr(vary_d1_coexist_sims, ~ .x[["NV"]])



coexist_sims <- readRDS(rds("coexist_sims"))

coexist_df <- map_dfr(coexist_sims, ~ .x[["NV"]])







#'
#' Code below shows that...
#' * Most are stable, except for the neutrally stable equilibria when
#'   d = 0.
#' * Some reps had complex eigenvalues, none with complex
#'   leading eigenvalues when very small imaginary components weren't allowed
#'   (anything less than 1e-10).
#'   A rep with a complex 15th eigenvalue (the closest to leading I could find)
#'   was plotted and didn't appear to have anything resembling fluctuations.
#'
#'
coexist_sims_eigens <- map_dfr(1:length(coexist_sims),
                               function(i) {
                                   eigs <- map_dbl(coexist_sims[[i]][["J"]],
                                                   function(.x) {
                                                       if (any(is.na(.x))) return(NA)
                                                       return(max(Re(eigen(.x, only.values = TRUE)$values)))
                                                   })
                                   coexist_sims[[i]]$NV %>%
                                       distinct(rep, d1, eta, add_var) %>%
                                       arrange(rep) %>%
                                       mutate(eigval = eigs)
                               })
coexist_sims_eigens %>%
    group_by(d1, eta, add_var) %>%
    summarize(n_unstable = sum(eigval > 1), max_eigen = max(eigval)) %>%
    ungroup() %>%
    filter(n_unstable > 0) %>%
    mutate(one_diff = max_eigen - 1) %>%
    filter(one_diff > 1e-10) %>%
    print_big_nums()




# --------------------*
# complex ----

coexist_sims_eigen_vals_cmplx <- map_dfr(1:length(coexist_sims),
                                         function(i) {
                                             eig_comps <- map_dbl(coexist_sims[[i]][["J"]],
                                                                  function(.x) {
                                                                      if (any(is.na(.x))) return(NA)
                                                                      ev <- eigen(.x, only.values = TRUE)$values
                                                                      if (is.complex(ev)) {
                                                                          return(max(abs(Im(ev))))
                                                                      } else return(NA_real_)
                                                                  })
                                             names(eig_comps) <- NULL
                                             cmplx_reps <- which(!is.na(eig_comps))
                                             if (length(cmplx_reps) == 0) {
                                                 out <- coexist_sims[[i]]$NV %>%
                                                     .[1,] %>%
                                                     select(rep, d1, eta, add_var) %>%
                                                     mutate(cmplx = 0.0) %>%
                                                     .[0,]
                                             } else {
                                                 out <- coexist_sims[[i]]$NV %>%
                                                     filter(rep %in% cmplx_reps) %>%
                                                     distinct(rep, d1, eta, add_var) %>%
                                                     mutate(cmplx = eig_comps[!is.na(eig_comps)])
                                             }
                                             return(out)
                                         })
coexist_sims_eigen_vals_cmplx %>%
    arrange(desc(cmplx))




coexist_sims_eigens_cmplx <- map_dfr(1:length(coexist_sims),
                                     function(i) {
                                         .t <- 1e-10
                                         eig_comps <- map_int(coexist_sims[[i]][["J"]],
                                                              function(.x) {
                                                                  if (any(is.na(.x))) return(NA)
                                                                  ev <- eigen(.x, only.values = TRUE)$values
                                                                  if (is.complex(ev)) {
                                                                      if (all(abs(Im(ev)) < .t)) {
                                                                          return(NA_integer_)
                                                                      }
                                                                      return(min(which(abs(Im(ev)) > .t)))
                                                                  } else return(NA_integer_)
                                                              })
                                         names(eig_comps) <- NULL
                                         cmplx_reps <- which(!is.na(eig_comps))
                                         if (length(cmplx_reps) == 0) {
                                             out <- coexist_sims[[i]]$NV %>%
                                                 .[1,] %>%
                                                 select(rep, d1, eta, add_var) %>%
                                                 mutate(cmplx = 1L) %>%
                                                 .[0,]
                                         } else {
                                             out <- coexist_sims[[i]]$NV %>%
                                                 filter(rep %in% cmplx_reps) %>%
                                                 distinct(rep, d1, eta, add_var) %>%
                                                 mutate(cmplx = eig_comps[!is.na(eig_comps)])
                                         }
                                         return(out)
                                     })
coexist_sims_eigens_cmplx %>%
    filter(cmplx == min(cmplx)) %>%
    group_by(d1, eta, add_var) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    print_big_nums()



.d1 = 0.9; .eta = 0.6; .add_var = 0.01


i <- map_lgl(coexist_sims,
             ~ .x[["NV"]] %>%
                 .[1,] %>%
                 with(d1 == .d1 & eta == .eta & add_var == .add_var)) %>%
    which()




.d <- c(.d1, 0.1)
.n <- 100

.seed <- coexist_sims[[i]][["seed"]]
set.seed(.seed)

Z <- quant_gen(q = 2, eta = .eta, d = .d, n = .n,
               add_var = rep(.add_var, .n),
               spp_gap_t = 500L, final_t = 50e3L, n_reps = 12,
               save_every = 10L,
               n_threads = .N_THREADS, show_progress = FALSE)

Z$call[["eta"]] <- eval(.eta)
Z$call[["d"]] <- eval(.d)
Z$call[["n"]] <- eval(.n)
Z$call[["add_var"]] <- eval(rep(.add_var, .n))

jacs <- jacobians(Z)
first_complex <- map_int(jacs,
                         function(x) {
                             eigs <- eigen(x, only.values = TRUE)[["values"]]
                             if (is.complex(eigs) && any(abs(Im(eigs)) > 1e-10)) {
                                 return(min(which(abs(Im(eigs)) > 1e-10)))
                             } else return(NA_integer_)
                         })
# map_dbl(jacs,
#         function(x) {
#             eigs <- eigen(x, only.values = TRUE)[["values"]]
#             max(Re(eigs))
#         })


Z %>%
    .[["nv"]] %>%
    # filter(rep %in% which(first_complex == 1), trait == 1) %>%
    # filter(rep %in% which(first_complex == 3), trait == 1) %>%
    filter(rep == 5) %>%
    filter(trait == 1) %>%
    filter(time > 75e3) %>%
    filter(N > 2) %>% .[["N"]] %>% range() %>% `-`(2.922151 - 3.008999e-07)

ggplot(aes(time, N)) +
    geom_line(aes(color = spp), na.rm = TRUE) +
    facet_wrap(~ rep, nrow = 3) +
    scale_color_viridis_d(begin = 0.1, end = 0.9, guide = FALSE)

# Z %>%
#     .[["nv"]] %>%
#     filter(rep %in% which(first_complex == 3), trait == 1) %>%
#     ggplot(aes(time, N)) +
#     geom_line(aes(color = spp), na.rm = TRUE) +
#     facet_wrap(~ rep, nrow = 3) +
#     scale_color_viridis_d(begin = 0.1, end = 0.9, guide = FALSE)


