#
#
# qg <- quant_gen(q = 2, eta = 0.6, d = 0,
#                 n_reps = 1, n = 1,
#                 V0 = cbind(c(0, 2)),
#                 sigma_V0 = 0,
#                 spp_gap_t = 0L,
#                 final_t = 200e3L,
#                 save_every = 0,
#                 add_var = rep(0.1, 1))
#
# qg$nv %>%
#     mutate(trait = paste0("V", trait)) %>%
#     spread(trait, geno) %>%
#     # .[["N"]]
#     identity()
#
#
# max(eigen(jacobians(qg)[[1]])$values)
#
#
# # .C <- diag(2)
# # .C[lower.tri(.C)] <- .C[upper.tri(.C)] <- 0.6
# #
# #
# max(eigen(sauron:::dVi_dVi_cpp(0, rbind(1.030776, 1.030776),
#                            2.679327,
#                            .C, formals(quant_gen)$f, formals(quant_gen)$a0,
#                            0.05))$values)
#
# deltaV <- sauron:::sel_str_cpp(V = rbind(0, 2), N = 10.91963, f = formals(quant_gen)$f,
#                      a0 = formals(quant_gen)$a0,
#                      C = .C, r0 = formals(quant_gen)$r0, D = diag(2) * 0) %*%
#     matrix(0.05)
#
# V <- rbind(0, 2)
# newV <- as.numeric(V + deltaV)
# heaviside <- function(x) ifelse(x > 0, 1, 0)
#
# jac <- sauron:::dVi_dVi_cpp(0, V,
#                             10.9196300066,
#                             .C, formals(quant_gen)$f, formals(quant_gen)$a0,
#                             0.05)
#
# max(eigen(diag(heaviside(newV)) %*% jac)$values)
#
#
# sauron:::F_t_cpp(rbind(0, 2),
#                  10.9196300066,
#                  formals(quant_gen)$f,
#                  formals(quant_gen)$a0,
#                  .C,
#                  formals(quant_gen)$r0,
#                  diag(2) * 0)
#
# max(eigen(sauron:::dVi_dVi_cpp(0, rbind(0, 2),
#                                10.9196300066,
#                                .C, formals(quant_gen)$f, formals(quant_gen)$a0,
#                                0.05))$values)
#
#
#
#
#
#
#
#
# # %>%
# #     ggplot(aes(V1, V2, color = spp)) +
# #     geom_path() +
# #     coord_equal(ylim = c(0, 2.5), xlim = c(0, 2.5))


# start ----

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
    library(viridisLite)
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
.REDO_SIMS <- FALSE
# whether to re-save plots
.RESAVE_PLOTS <- TRUE
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




grab_sims <- function(.d, .eta, .add_var, .sigma_N, .sigma_V, .vary_d2) {

    .d <- round(.d, 3)
    .eta <- round(.eta, 2)
    .add_var <- round(.add_var, 2)
    .sigma_N <- round(.sigma_N, 2)
    .sigma_V <- round(.sigma_V, 2)

    giant_sims <- mclapply(0:17,
                           function(i) {
                               fn <- rds(sprintf("giant_sims/giant_sims_%i", i))
                               dd <- readRDS(fn)
                               map_dfr(dd, function(Z) {
                                   Z[["NV"]] %>%
                                       filter(d %in% .d,
                                              eta %in% .eta,
                                              add_var %in% .add_var,
                                              sigma_N %in% .sigma_N,
                                              sigma_V %in% .sigma_V,
                                              vary_d2 %in% .vary_d2)
                               })
                           },
                           mc.cores = .N_THREADS)

    giant_sims <- bind_rows(giant_sims)

    return(giant_sims)
}

set.seed(1898348146)
etas <- map(1:6, ~ with(list(.q = .x), runif(.q * ((.q - 1) / 2), 0.1, 0.4)))
# With just one eta, it can be a simple number:
etas[[2]] <- 0.6


# =============================================================================*
# =============================================================================*

# Fig 2. 2-trait outcomes ----

# =============================================================================*
# =============================================================================*



one_eta_combo <- function(signs, .d = 1, ...) {


    # signs = -1; .d = 1
    # rm(signs, .d)

    stopifnot(is.numeric(signs))
    stopifnot(sum(is.na(signs)) == 0)
    stopifnot(all(signs %in% c(-1, 0, 1)))

    .q <- 1/2 * (1 + sqrt(1 + 8 * length(signs)))

    stopifnot(.q %in% 1:length(etas))

    .etas <- etas[[.q]]

    stopifnot(length(signs) == length(.etas))

    stopifnot(is.numeric(.d))
    stopifnot(sum(is.na(.d)) == 0)
    stopifnot(length(.d) %in% c(1L, .q))

    C <- matrix(0, .q, .q)
    C[lower.tri(C)] <- abs(.etas) * signs
    C <- C + t(C)
    diag(C) <- 1

    args <- list(q = .q, eta = C, d = .d, n_reps = 24,
                 spp_gap_t = 500L, final_t = 20e3L, save_every = 0L,
                 sigma_V0 = 1, n_threads = .N_THREADS,
                 show_progress = FALSE)

    other_args <- list(...)
    if (length(other_args) > 0) {
        stopifnot(!is.null(names(other_args)))
        stopifnot(all(names(other_args) != ""))
        for (n in names(other_args)) args[[n]] <- other_args[[n]]
    }

    trait_to <- do.call(quant_gen, args)


    NV <- trait_to$nv %>%
        mutate(trait = paste0("V", trait)) %>%
        select(rep, spp, trait, geno) %>%
        spread(trait, geno)
    cn <- colnames(NV)[grepl("^V", colnames(NV))]

    NV <- NV %>%
        mutate(unq = do.call(unq_spp_filter, !!rlang::syms(cn))) %>%
        filter(unq) %>%
        select(!!!cn) %>%
        identity()

    for (i in 1:length(.etas)) {
        NV[,paste0("eta",i)] <- abs(.etas[i]) * signs[i]
    }

    output <- list(jacs = jacobians(trait_to),
                   ts = NV)


    return(output)
}



if (.REDO_SIMS) {
    # Takes ~5 sec with q=2 and 3 threads
    set.seed(145746935)
    eta_sims <- list(-1, 0, 1) %>%
        map(one_eta_combo)
    saveRDS(eta_sims, rds("eta_sims"))
} else {
    eta_sims <- readRDS(rds("eta_sims"))
}

eta_sim_df <- map_dfr(eta_sims, ~.x$ts)







#'
#' Based on eigenvalues (see `_main-results__stability.R`)...
#'   * When the tradeoff is additive, the state is neutrally stable
#'   * Everything else is stable
#'
#'




outcomes_q2_p <- eta_sim_df %>%
    mutate(eta1 = factor(eta1, levels = sort(unique(eta1)),
                         labels = c("sub-additive", "additive",
                                    "super-additive"))) %>%
    group_by(eta1) %>%
    filter(unq_spp_filter(V1, V2, .prec = 0.1)) %>%
    ungroup() %>%
    ggplot(aes(V1, V2)) +
    geom_path(data = stable_points(0) %>% mutate(eta1 = factor("additive")),
              linetype = 1, color = "gray40") +
    geom_point(shape = 21, color = "dodgerblue", fill = "dodgerblue",
               alpha = 0.75, size = 3) +
    scale_x_continuous("Axis 1", breaks = 0:2) +
    scale_y_continuous("Axis 2", breaks = 0:2) +
    coord_equal(xlim = c(0, 2.5), ylim = c(0, 2.5)) +
    facet_wrap(~ eta1, nrow = 1) +
    theme(strip.text = element_text(size = 10),
          panel.border = element_rect(size = 0.5, fill = NA)) +
    NULL



if (.RESAVE_PLOTS) save_plot(outcomes_q2_p, 5, 2, .prefix = "2-")





# =============================================================================*
# =============================================================================*

# Fig 3. Coexistence ----

# =============================================================================*
# =============================================================================*




# ------------------------------------*
# __ vary d1 & d2 ----
# ------------------------------------*

# Simulations varying both d values

coexist_spp_df <- grab_sims(.d = seq(-0.25, 2, length.out = 10),
                            .eta = -1:1 * etas[[2]],
                            .add_var = c(0.01, 0.05, 0.1),
                            .sigma_N = 0,
                            .sigma_V = 0,
                            .vary_d2 = TRUE) %>%
    group_by(d, eta, add_var, rep) %>%
    summarize(n_spp = spp[N > 0] %>% unique() %>% length()) %>%
    ungroup() %>%
    mutate(eta = factor(eta, levels = -1:1 * etas[[2]],
                        labels = paste0(c("sub-", "", "super-"), "additive")))


coexist_spp_p1 <- coexist_spp_df %>%
    mutate(n_spp = n_spp / 100,
           add_var = factor(add_var, levels = sort(unique(add_var)))) %>%
    ggplot(aes(d, n_spp, color = add_var)) +
    geom_vline(xintercept = 0, color = "gray80", linetype = 1) +
    geom_hline(yintercept = 0, color = "gray80", linetype = 1) +
    geom_jitter(width = 0.02, height = 0, shape = 1, size = 0.5) +
    # geom_text(data = coexist_spp_df %>%
    #               filter(d == 0, eta == "sub-additive") %>%
    #               group_by(eta, add_var) %>%
    #               summarize(n_spp = median(n_spp / 100)) %>%
    #               ungroup() %>%
    #               mutate(d = -0.1,
    #                      lab = sprintf("italic(sigma[i])^2 == %.2f", add_var),
    #                      add_var = factor(add_var)),
    #           aes(label = lab), parse = TRUE, hjust = 1, vjust = 0.5,
    #           size = 10 * 0.352778) +
    facet_wrap(~ eta, ncol = 1) +
    scale_color_viridis_d(expression(italic(sigma[i])^2),
                          begin = 0.1, end = 0.85, guide = FALSE) +
    scale_y_continuous("Proportion of species that survive",
                       breaks = c(0, 0.4, 0.8), limits = c(0, 1)) +
    scale_x_continuous(expression(italic(d[1]) ~ "and" ~ italic(d[2]))) +
    NULL


#'
#' Code in `_main-results__stability.R` shows that...
#' * Most are stable, except for the neutrally stable equilibria when
#'   d = 0.
#' * Some reps had complex eigenvalues, none with complex
#'   leading eigenvalues when very small imaginary components weren't allowed
#'   (anything less than 1e-10).
#'   A rep with a complex 15th eigenvalue (the closest to leading I could find)
#'   was plotted and didn't appear to have anything resembling fluctuations.
#'
#'



# ------------------------------------*
# __ vary d1 only ----
# ------------------------------------*



vary_d1_coexist_spp_df <- grab_sims(.d = seq(-0.15, 0.05, 0.05),
                                    .eta = -1:1 * etas[[2]],
                                    .add_var = c(0.01, 0.05, 0.1),
                                    .sigma_N = 0,
                                    .sigma_V = 0,
                                    .vary_d2 = FALSE) %>%
    group_by(d, eta, add_var, rep) %>%
    summarize(n_spp = spp[N > 0] %>% unique() %>% length()) %>%
    ungroup() %>%
    mutate(eta = factor(eta, levels = -1:1 * etas[[2]],
                        labels = paste0(c("sub-", "", "super-"), "additive")))




coexist_spp_p2 <- vary_d1_coexist_spp_df %>%
    mutate(n_spp = n_spp / 100,
           add_var = factor(add_var)) %>%
    ggplot(aes(d, n_spp, color = add_var)) +
    geom_vline(xintercept = 0, color = "gray80", linetype = 1) +
    geom_hline(yintercept = 0, color = "gray80", linetype = 1) +
    geom_jitter(width = 0.002, height = 0, shape = 1, size = 1) +
    geom_text(data = vary_d1_coexist_spp_df %>%
                  filter(d == -0.1, eta == "sub-additive") %>%
                  group_by(eta, add_var) %>%
                  summarize(n_spp = median(n_spp / 100)) %>%
                  ungroup() %>%
                  mutate(d = -0.11,
                         lab = sprintf("italic(sigma[i])^2 == %.2f", add_var),
                         add_var = factor(add_var)),
              aes(label = lab), parse = TRUE, hjust = 1, vjust = 0.5,
              size = 10 * 0.352778) +
    facet_wrap(~ eta, ncol = 1) +
    scale_color_viridis_d(expression(italic(sigma[i])^2),
                          begin = 0.1, end = 0.85, guide = FALSE) +
    scale_y_continuous("Proportion of species that survive",
                       breaks = c(0, 0.4, 0.8), limits = c(0, 1)) +
    scale_x_continuous(expression(italic(d[1]) ~ ("with" ~ italic(d[2]) == 0.1)),
                       limits = c(-0.17, 0.06))







coexist_spp_p <- plot_grid(coexist_spp_p1,
                           coexist_spp_p2 +
                               theme(axis.title.y = element_blank()),
                           nrow = 1, labels = LETTERS[1:2], align = "vh",
                           label_fontface = "plain")

# coexist_spp_p


if (.RESAVE_PLOTS) {
    save_plot(coexist_spp_p, 6, 6, .prefix = "3-")
}












# =============================================================================*
# =============================================================================*

# Figs 4,S1-S2 Conditional coexistence ----

# =============================================================================*
# =============================================================================*




cond_coexist_test <- function(.V0, .eta_sign, .d1, .sigma_V = 0, .sigma_N = 0) {

    # .V0 = "restricted"; .eta_sign = 0; .d1 = -0.1; .sigma_V = 0; .sigma_N = 0
    # rm(.V0, .eta_sign, .d1, .sigma_V, .sigma_N)

    .dist <- 0.6
    which_switch <- 3 # which spp to switch for unrestricted

    .n <- 5
    .q = 2
    stopifnot(is.numeric(.d1) && length(.d1) == 1)
    .ds <- c(.d1, abs(.d1))
    stopifnot(is.numeric(.eta_sign) && length(.eta_sign) == 1 &&
                  .eta_sign %in% -1:1)
    .V0 <- match.arg(.V0, c("restricted", "unrestricted"))
    .lab <- .V0

    .eta <- .eta_sign * etas[[2]]

    if (.lab == "restricted" && .eta < 0) {
        stop(paste("\n'restricted' with sub-additivity is not programmed bc",
                   "there's only one stable trait state.",
                   "Thus, there is no way to restrict starting axes to",
                   "be outside of that state's basin of attraction."))
    }

    if (.eta < 0) {
        .V0 <- rbind(seq(3, 2, length.out = 5),
                     seq(2, 3, length.out = 5))
    } else if (.eta == 0) {
        .V0 <- rbind(seq(1.5, 0, length.out = 5),
                     seq(1.6, 3.1, length.out = 5))
        if (.lab == "unrestricted") {
            v2 <- .V0[2, which_switch]
            .V0[2, which_switch] <- .V0[1, which_switch]
            .V0[1, which_switch] <- v2
        }
    } else {
        .V0 <- rbind(seq(1.2, 0, length.out = 5),
                     seq(1.3, 2.5, length.out = 5))
        if (.lab == "unrestricted") {
            v2 <- .V0[2, which_switch]
            .V0[2, which_switch] <- .V0[1, which_switch]
            .V0[1, which_switch] <- v2
        }
    }
    .V0 <- round(.V0, 3)

    if (.sigma_V == 0 && .sigma_N == 0) {
        .nreps <- 1
        .N_THREADS <- 1
    } else .nreps <- 12


    qg <- quant_gen(q = .q, eta = .eta, d = .ds,
                   n_reps = .nreps, n = ncol(.V0),
                   V0 = .V0,
                   sigma_V0 = 0,
                   sigma_V = .sigma_V,
                   sigma_N = .sigma_N,
                   spp_gap_t = 500L,
                   final_t = 20e3L,
                   add_var = rep(0.05, .n),
                   n_threads = .N_THREADS,
                   show_progress = FALSE)


    if (.sigma_V == 0) {

        out <- list(nv = qg %>%
                        .[["nv"]] %>%
                        mutate(trait = paste0("V", trait)) %>%
                        spread(trait, geno))
        if (.sigma_N == 0) {
            out$nv <- select(out$nv, -rep)

            qg$call[["q"]] <- eval(.q)
            qg$call[["n"]] <- eval(.n)
            qg$call[["d"]] <- eval(.ds)
            qg$call[["eta"]] <- eval(.eta)
            qg$call[["add_var"]] <- eval(rep(0.05, .n))
            out$jacs <- jacobians(qg)
        }


    } else {

        out <- list(nv = qg %>%
                        .[["nv"]] %>%
                        mutate(trait = paste0("V", trait)) %>%
                        nest(value = c(geno, pheno)) %>%
                        spread(trait, value) %>%
                        unnest(c(V1, V2), names_sep = "_") %>%
                        rename(V1 = V1_geno,
                               V2 = V2_geno,
                               Vp1 = V1_pheno,
                               Vp2 = V2_pheno))

    }

    out$nv <- out$nv %>%
        mutate(V0 = .lab, eta = .eta, d1 = .d1,
               sigma_V = .sigma_V, sigma_N = .sigma_N)

    return(out)
}


if (.REDO_SIMS) {

    # Takes ~6 sec
    cond_coexist <- crossing(.V0 = c("restricted", "unrestricted"),
                             .eta_sign = -1:1,
                             .d1 = c(-0.1, 0.1)) %>%
        filter(!(.V0 == "restricted" & .eta_sign < 0)) %>%
        pmap(cond_coexist_test)
    saveRDS(cond_coexist, rds("cond_coexist"))

    # Takes ~3.5 min
    set.seed(351367879)
    cond_coexist_sV <- crossing(.V0 = c("restricted", "unrestricted"),
                             .eta_sign = -1:1,
                             .d1 = c(-0.1, 0.1),
                             .sigma_V = c(0.05, 0.1, 0.2)) %>%
        filter(!(.V0 == "restricted" & .eta_sign < 0)) %>%
        pmap(cond_coexist_test)
    saveRDS(cond_coexist_sV, rds("cond_coexist_sV"))

    # Takes ~1.7 min
    set.seed(602553504)
    cond_coexist_sN <- crossing(.V0 = c("restricted", "unrestricted"),
                             .eta_sign = -1:1,
                             .d1 = c(-0.1, 0.1),
                             .sigma_N = c(0.05, 0.1, 0.2)) %>%
        filter(!(.V0 == "restricted" & .eta_sign < 0)) %>%
        pmap(cond_coexist_test)
    saveRDS(cond_coexist_sN, rds("cond_coexist_sN"))

} else {
    cond_coexist <- readRDS(rds("cond_coexist"))
    cond_coexist_sV <- readRDS(rds("cond_coexist_sV"))
    cond_coexist_sN <- readRDS(rds("cond_coexist_sN"))
}


cond_coexist_df_prep <- function(.dd) {
    .dd %>%
        map_dfr(~ .x$nv) %>%
        # filter(time < 4e3L) %>%
        mutate(V0 = factor(V0, levels = c("restricted", "unrestricted")),
               eta = factor(eta, levels = -1:1 * etas[[2]],
                            labels = c("sub-additive", "additive",
                                       "super-additive")),
               d1 = factor(d1, levels = c(-0.1, 0.1),
                           labels = c("conflicting",
                                      "ameliorative"))) %>%
        # `trait_space` is a combination of starting trait values and eta:
        mutate(trait_space =
                   case_when(V0 == "unrestricted" & eta == "sub-additive" ~
                                 "i",
                             V0 == "restricted" & eta == "super-additive" ~
                                 "ii",
                             V0 == "restricted" & eta == "additive" ~
                                 "iii",
                             V0 == "unrestricted" & eta == "super-additive" ~
                                 "iv",
                             V0 == "unrestricted" & eta == "additive" ~
                                 "v",
                             TRUE ~ NA_character_) %>%
                   factor(levels = c("i", "ii", "iii", "iv", "v")))
}



cond_coexist_df <- cond_coexist %>%
    cond_coexist_df_prep() %>%
    select(-starts_with("sigma_"))

# Reps with stochasticity:
cond_coexist_stoch_df <- map_dfr(list(cond_coexist_sV,
                                       cond_coexist_sN),
                                  cond_coexist_df_prep)


if (any(is.na(cond_coexist_df$trait_space))) {
    stop("\nERROR: unknown combination of V0 and eta in `cond_coexist_df`")
}
if (any(is.na(cond_coexist_stoch_df$trait_space))) {
    stop("\nERROR: unknown combination of V0 and eta in `cond_coexist_stoch_df`")
}



# Starting conditions and trajectories:


cond_coexist_sc_p_fun <- function(.d1) {
    .dd <- cond_coexist_df %>%
        filter(d1 == .d1) %>%
        split(interaction(.$trait_space, .$spp, drop = TRUE)) %>%
        map_dfr(~ mutate(.x,
                         first = time == min(time),
                         last = time == max(time))) %>%
        select(trait_space, time, spp, V1, V2, first, last) %>%
        arrange(trait_space, time)
    .dd %>%
        ggplot(aes(V1, V2)) +
        geom_abline(data = tibble(trait_space = factor(c("ii", "iv")),
                                  slp = 1, int = 0),
                    aes(slope = slp, intercept = int), linetype = 3, color = "gray70") +
        geom_line(data = bind_rows(stable_points(0), stable_points(0)) %>%
                      mutate(trait_space = factor(rep(c("iii", "v"),
                                                      each = stable_points %>%
                                                          formals() %>%
                                                          .[["line_n"]]))),
                  linetype = 2, color = "black") +
        geom_point(data = .dd %>% filter(first),
                   aes(color = spp), size = 1.5, na.rm = TRUE) +
        geom_point(data = map_dfr(etas[[2]] * c(-1, 1, 1), ~ stable_points(.x)) %>%
                       mutate(trait_space = map2(c("i", "ii", "iv"), c(1,2,2), rep) %>%
                                  do.call(what = c) %>%
                                  factor(),
                              shp = case_when(V1 > 0 & V2 > 0 ~ 3L,
                                              V1 == 0 ~ 2L,
                                              TRUE ~ 1L) %>%
                                  factor(levels = 1:3)),
                   aes(shape = shp), size = 3, color = "black") +
        geom_point(data = .dd %>%
                       filter(time < max(time), last),
                   aes(color = spp), shape = 4, size = 1.5) +
        geom_path(aes(color = spp)) +
        scale_color_brewer(palette = "Dark2", guide = FALSE) +
        scale_shape_manual(values = c(5,1,2), guide = FALSE) +
        scale_size_continuous(range = c(0.1, 1)) +
        facet_grid(~ trait_space) +
        coord_equal(xlim = c(-0.1, 3.15), ylim = c(-0.1, 3.15)) +
        ylab("Axis 2") +
        xlab("Axis 1") +
        theme(plot.margin = margin(0,0,0,b=6),
              plot.title = element_text(hjust = 0.5,
                                        margin = margin(0,0,0,b=6)))
}

cond_coexist_sc_p <- cond_coexist_sc_p_fun("conflicting")



# Time series of abundances
cc_N_p <- cond_coexist_df %>%
    filter(d1 == "conflicting") %>%
    filter(time < 7e3) %>%
    ggplot(aes(time / 1000L, N, color = spp)) +
    geom_line() +
    geom_point(data = cond_coexist_df %>%
                   filter(d1 == "conflicting") %>%
                   group_by(trait_space, spp) %>%
                   filter(time == max(time)) %>%
                   ungroup() %>%
                   filter(time < max(time)),
               shape = 4, size = 1.5) +
    facet_wrap(~ trait_space, nrow = 1) +
    scale_color_brewer(palette = "Dark2",
                       guide = FALSE) +
    scale_y_continuous("Abundance", trans = "log",
                       breaks = 10^(c(-3, 0, 3)),
                       labels = parse(
                           text = sprintf("10^{%i}",
                                          c(-3, 0, 3)))) +
    xlab("Time (× 1,000)") +
    theme(plot.margin = margin(0,r=12,t=10,b=10),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    margin = margin(0,0,0,b=6)))





stable_state_df <- map_dfr(c(1:2, 4),
                           function(i) {
                               ts <- cond_coexist_df$trait_space %>%
                                   unique() %>%
                                   sort() %>%
                                   .[i]
                               stable_points((etas[[2]] * c(-1,1,0,1,0))[i]) %>%
                                   mutate(trait_space = ts)
                           }) %>%
    mutate(time = max(cond_coexist_df$time[cond_coexist_df$time < 7e3]),
           shp = case_when(V1 > 0 & V2 > 0 ~ 3L,
                           V1 == 0 ~ 2L,
                           TRUE ~ 1L) %>%
               factor(levels = 1:3)) %>%
    gather(trait, value, V1:V2) %>%
    mutate(trait = gsub("^V", "trait ", trait))



dist_from_equil <- function(V1, V2, eta) {
    .eta <- case_when(eta[1] == "additive" ~ 0,
                      eta[1] == "sub-additive" ~ -etas[[2]],
                      TRUE ~ etas[[2]])
    if (.eta == 0) {
        eq_radius <- with(formals(quant_gen), sqrt((r0 / f) - 1))
        R <- sqrt(V1^2 + V2^2)
        dist <- abs(R - eq_radius)
    } else {
        equil <- stable_points(.eta)
        dist <- matrix(NA_real_, length(V1), nrow(equil))
        for (j in 1:nrow(equil)) {
            dist[,j] <- sqrt((V1 - equil$V1[j])^2 + (V2 - equil$V2[j])^2)
        }
        dist <- apply(dist, 1, min)
    }
    return(dist)
}



cc_V_p <- cond_coexist_df %>%
    filter(d1 == "conflicting") %>%
    filter(time < 7e3) %>%
    split(.$eta) %>%
    map_dfr(function(.dd) {
        mutate(.dd, dist = dist_from_equil(V1, V2, eta))
    }) %>%
    mutate(dist = ifelse(eta == "additive", dist, mean(dist))) %>%
    gather(trait, value, V1:V2) %>%
    mutate(trait = gsub("^V", "trait ", trait)) %>%
    ggplot(aes(time / 1000L, value)) +
    geom_hline(yintercept = 0, size = 0.5,
               linetype = 1, color = "gray70") +
    geom_vline(xintercept = 0, size = 0.5,
               linetype = 1, color = "gray70") +
    geom_point(data = stable_state_df,
               aes(shape = shp), size = 3) +
    geom_line(aes(color = spp, size = dist)) +
    geom_line(aes(color = spp)) +
    geom_point(data = cond_coexist_df %>%
                   filter(d1 == "conflicting") %>%
                   gather(trait, value, V1:V2) %>%
                   mutate(trait = gsub("^V", "trait ", trait)) %>%
                   group_by(trait_space, spp) %>%
                   filter(time == max(time)) %>%
                   ungroup() %>%
                   filter(time < max(time)),
               aes(color = spp), shape = 4, size = 1.5) +
    facet_grid(trait ~ trait_space) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    scale_shape_manual(values = c(5,1,2), guide = FALSE) +
    scale_size_continuous(range = c(0.4, 2), guide = FALSE) +
    scale_y_continuous("Axis value", limits = c(-0.2, NA)) +
    scale_x_continuous("Time (× 1,000)",
                       limits = c(0, 7.2)) +
    theme(plot.margin = margin(0,0,0,t=10),
          plot.title = element_text(hjust = 0.5,
                                    margin = margin(0,0,0,b=6)))








cond_coexist_p <- plot_grid(cond_coexist_sc_p +
                                xlab(sprintf("Axis 1\n(%s)", "conflicting")) +
                                ylab("(ameliorative)\nAxis 2") +
                                ggtitle(paste("Starting conditions and",
                                              "trajectories")) +
                                theme(plot.margin = margin(0,0,0,r=12)),
                            cc_N_p +
                                ggtitle("Abundance time series"),
                            cc_V_p +
                                ggtitle("Axis time series"),
                            align = "v", axis = "l",
                            ncol = 1, rel_heights = c(3, 2, 4),
                            labels = LETTERS[1:3],
                            label_fontface = "plain", label_x = 0.06)



if (.RESAVE_PLOTS) {
    save_plot(cond_coexist_p, 6, 7, .name = "4-cond_coexist")
}




#'
#' They're all stable except for sub-additive and conflicting trait 1,
#' which is neutrally stable (see `_main-results__stability.R`).
#'




# ... with stoch. ----


cc_N_stoch_plot_fun <- function(.V_stoch, .ts = FALSE) {

    .d1 = "conflicting"

    # .V_stoch = TRUE; .ts = FALSE
    # rm(.d1, .V_stoch, .ts, .dd, .dd2)

    stopifnot(is.logical(.V_stoch) && length(.V_stoch) == 1)

    if (.V_stoch) {
        .sigma_V <- 0.1
        .sigma_N <- 0
    } else {
        .sigma_V <- 0
        .sigma_N <- 0.1
    }

    stopifnot((.sigma_V > 0 || .sigma_N) > 0 && !(.sigma_V > 0 && .sigma_N > 0))

    .dd <- cond_coexist_stoch_df %>%
        filter(d1 == .d1,
               sigma_N == .sigma_N,
               sigma_V == .sigma_V)

    if (.ts) {
        .dd2 <- cond_coexist_df %>%
            filter(d1 == .d1) %>%
            mutate(rep = 0L) %>%
            select(trait_space, rep, time, spp, N)
        .dd %>%
            mutate(rep = rep %>% paste() %>% as.integer()) %>%
            bind_rows(.dd2) %>%
            mutate(rep = factor(rep, levels = 0:12)) %>%
            mutate(id = interaction(spp, rep)) %>%
            ggplot(aes(time / 1000L, N, color = spp)) +
            geom_line(aes(group = id)) +
            geom_point(data = .dd %>%
                           group_by(trait_space, rep, spp) %>%
                           filter(time == max(time)) %>%
                           ungroup() %>%
                           filter(time < max(time)),
                       shape = 4, size = 1.5) +
            facet_grid(rep ~ trait_space) +
            scale_color_brewer(palette = "Dark2",
                               guide = FALSE) +
            scale_y_continuous("Abundance", trans = "log",
                               breaks = 10^(c(-3, 0, 3)),
                               labels = parse(
                                   text = sprintf("10^{%i}",
                                                  c(-3, 0, 3)))) +
            xlab("Time (× 1,000)") +
            theme(strip.text.y = element_blank())
    } else {

        .dd <- .dd %>%
            filter(time == max(time)) %>%
            select(-time)

        .dd2 <- cond_coexist_df %>%
            filter(d1 == .d1, time == max(time)) %>%
            group_by(trait_space) %>%
            mutate(tN = sqrt(N) / sum(sqrt(N))) %>%
            ungroup() %>%
            mutate(rep = 0L) %>%
            select(trait_space, rep, spp, N, tN)

        .dd3 <- map_dfr(list(.dd, .dd2),
                        ~ .x %>%
                            group_by(trait_space, rep) %>%
                            summarize(n_spp = sum(N > 0)) %>%
                            ungroup() %>%
                            mutate(rep = factor(rep %>% paste() %>% as.integer(),
                                                levels = 0:12)))

        .dd %>%
            group_by(trait_space, rep) %>%
            mutate(tN = sqrt(N) / sum(sqrt(N))) %>%
            # mutate(tN = N / sum(N)) %>%
            ungroup() %>%
            select(trait_space, rep, spp, N, tN) %>%
            mutate(rep = rep %>% paste() %>% as.integer()) %>%
            bind_rows(.dd2) %>%
            mutate(rep = factor(rep, levels = 0:12)) %>%
            ggplot(aes(rep)) +
            geom_bar(aes(weight = tN, fill = spp), color = "white", size = 0.1) +
            geom_text(data = .dd3 %>% filter(n_spp > 1),
                      aes(label = n_spp, y = 1.05),
                      size = 8 / 2.83465) +
            geom_vline(xintercept = 1.5, size = 0.5) +
            facet_grid( ~ trait_space) +
            scale_fill_brewer(palette = "Dark2", guide = FALSE) +
            scale_y_continuous("Scaled relative abundance",
                               breaks = c(0, 0.5, 1)) +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  plot.title = element_text(hjust = 0.5,
                                            margin = margin(0,0,0,b=6),
                                            size = 11)) +
            NULL
    }
}


cond_coexist_stoch_ps <- map(c(TRUE, FALSE), cc_N_stoch_plot_fun)
names(cond_coexist_stoch_ps) <- c("V_stoch", "N_stoch")

cond_coexist_stoch_p <- ggarrange(cond_coexist_stoch_ps[["N_stoch"]] +
                                      ggtitle("With abundance stochasticity") +
                                      theme(plot.margin = margin(0,0,0,b=12),
                                            axis.title.y = element_blank()),
                                  cond_coexist_stoch_ps[["V_stoch"]] +
                                      ggtitle("With trait stochasticity") +
                                      theme(plot.margin = margin(0,0,0,0),
                                            axis.title.y = element_blank()),
                                  ncol = 1,
                                  padding = unit(1, "lines"),
                                  left = "Scaled relative abundance",
                                  draw = FALSE, labels = LETTERS[1:2],
                                  label.args = list(gp = grid::gpar(
                                      fontface = "plain", fontsize = 12),
                                      hjust = -2))

# cond_coexist_stoch_p

#'
#' This helps explain why you get totally different outcomes for situation v
#' when sigma_V > 0, compared to when sigma_V = 0.
#'
#' This is NOT explained by the distribution for V stochasticity
#' being lognormal and having a mean > 1.
#' I tried with simulations where the distribution was corrected for this
#' (and did indeed have a mean of 1).
#' The results were the same.
#'
cc_sigmaV_sit_v_p <- cond_coexist_stoch_df %>%
    filter(d1 == "conflicting",
           sigma_N == 0,
           sigma_V == 0.1) %>%
    filter(trait_space == "v") %>%
    mutate(id = interaction(spp, rep)) %>%
    arrange(time) %>%
    ggplot(aes(V1, V2, color = spp)) +
    geom_path(aes(group = id)) +
    geom_point(data = cond_coexist_stoch_df %>%
                   filter(d1 == "conflicting",
                          sigma_N == 0,
                          sigma_V == 0.1) %>%
                   filter(trait_space == "v") %>%
                   group_by(spp) %>%
                   summarize(V1 = V1[time == min(time)][1],
                             V2 = V2[time == min(time)][1])) +
    facet_wrap(~ rep, nrow = 3) +
    coord_equal(xlim = c(0, 3.11), ylim = c(0, 3.11)) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    ylab("(ameliorative)\nAxis 2") +
    xlab("Axis 1\n(conflicting)")




if (.RESAVE_PLOTS) {
    save_plot(cond_coexist_stoch_p, 6, 4, .prefix = "S1-")
    save_plot(cc_sigmaV_sit_v_p, 6, 5, .prefix = "S2-")
}








# =============================================================================*
# =============================================================================*

# Fig S3 Stoch. - # species ----

# =============================================================================*
# =============================================================================*


# Simulations varying both d values


stoch_coexist_spp_df <- grab_sims(.d = seq(-0.25, 2, length.out = 10),
                                  .eta = -1:1 * etas[[2]],
                                  .add_var = c(0.01, 0.05, 0.1),
                                  .sigma_N = c(0, 0.05, 0.1, 0.2, 0.3),
                                  .sigma_V = c(0, 0.05, 0.1, 0.2, 0.3),
                                  .vary_d2 = TRUE) %>%
    group_by(d, eta, add_var, sigma_N, sigma_V, rep) %>%
    summarize(n_spp = spp[N > 0] %>% unique() %>% length()) %>%
    ungroup() %>%
    mutate(eta = factor(eta, levels = -1:1 * etas[[2]],
                        labels = c("sub-additive", "additive", "super-additive")))


stoch_coexist_p_fun <- function(.x) {
    stoch_coexist_spp_df %>%
        filter(eta == .x) %>%
        filter(sigma_N %in% c(0, 0.1, 0.2),
               sigma_V %in% c(0, 0.1, 0.2),
               add_var == 0.05) %>%
        mutate(sigma_N = factor(sigma_N, levels = sort(unique(sigma_N)),
                                labels = sprintf("sigma[N] == %.2f",
                                                 sort(unique(sigma_N)))),
               sigma_V = factor(sigma_V, levels = sort(unique(sigma_V)),
                                labels = sprintf("sigma[V] == %.2f",
                                                 sort(unique(sigma_V)))),
               add_var = factor(add_var, levels = sort(unique(add_var))),
               extinct = factor(n_spp == 0)) %>%
        mutate(n_spp = n_spp / 100) %>%
        # ggplot(aes(d, n_spp, color = add_var)) +
        ggplot(aes(d, n_spp)) +
        geom_hline(yintercept = 0, color = "gray80") +
        geom_vline(xintercept = 0, color = "gray80") +
        geom_jitter(aes(shape = extinct, size = extinct),
                    color = "dodgerblue",
                    width = 0.02, height = 0) +
        ggtitle(.x) +
        facet_grid(sigma_N ~ sigma_V, labeller = label_parsed) +
        # scale_color_viridis_d(expression(italic(sigma[i])^2),
        #                       begin = 0.1, end = 0.85) +
        scale_shape_manual(values = c(1, 4), guide = FALSE) +
        scale_size_manual(values = c(0.5, 2), guide = FALSE) +
        scale_y_continuous("Proportion of species that survive",
                           breaks = c(0, 0.4, 0.8), limits = c(0, 1)) +
        scale_x_continuous(expression(italic(d[1]) ~ "and" ~ italic(d[2])),
                           breaks = c(0, 1, 2)) +
        guides(color = guide_legend(override.aes = list(size = 2, shape = 19))) +
        theme(strip.text = element_text(size = 10),
              strip.text.y = element_text(angle = 0),
              plot.title = element_text(hjust = 0.5,
                                        margin = margin(0,0,0,b=12))) +
        NULL
}

stoch_coexist_ps <- map(c("sub-additive", "additive", "super-additive"),
                        stoch_coexist_p_fun)
names(stoch_coexist_ps) <- c("sub", "add", "super")

# stoch_coexist_ps[["sub"]]
# stoch_coexist_ps[["add"]]
# stoch_coexist_ps[["super"]]



# Now looking at it near the boundaries, only varying d1:


stoch_vary_d1_coexist_spp_df <- grab_sims(.d = seq(-0.15, 0.05, 0.025),
                                          .eta = -1:1 * etas[[2]],
                                          .add_var = c(0.01, 0.05, 0.1),
                                          .sigma_N = c(0, 0.05, 0.1, 0.2, 0.3),
                                          .sigma_V = c(0, 0.05, 0.1, 0.2, 0.3),
                                          .vary_d2 = FALSE) %>%
    group_by(d, eta, add_var, sigma_N, sigma_V, vary_d2, rep) %>%
    summarize(n_spp = spp[N > 0] %>% unique() %>% length()) %>%
    ungroup() %>%
    mutate(eta = factor(eta, levels = -1:1 * etas[[2]],
                        labels = c("sub-additive", "additive", "super-additive")))


stoch_coexist_d1_p_fun <- function(.x) {

    stoch_vary_d1_coexist_spp_df %>%
        filter(eta == .x) %>%
        filter(sigma_N %in% c(0, 0.1, 0.2),
               sigma_V %in% c(0, 0.1, 0.2),
               add_var == 0.05) %>%
        mutate(sigma_N = factor(sigma_N, levels = sort(unique(sigma_N)),
                                labels = sprintf("sigma[N] == %.2f",
                                                 sort(unique(sigma_N)))),
               sigma_V = factor(sigma_V, levels = sort(unique(sigma_V)),
                                labels = sprintf("sigma[V] == %.2f",
                                                 sort(unique(sigma_V)))),
               add_var = factor(add_var, levels = sort(unique(add_var))),
               extinct = factor(n_spp == 0)) %>%
        mutate(n_spp = n_spp / 100) %>%
        # ggplot(aes(d, n_spp, color = add_var)) +
        ggplot(aes(d, n_spp)) +
        ggtitle(.x) +
        geom_hline(yintercept = 0, color = "gray80") +
        geom_vline(xintercept = 0, color = "gray80") +
        geom_jitter(aes(shape = extinct, size = extinct), width = 0.002,
                    height = 0, color = "firebrick") +
        facet_grid(sigma_N ~ sigma_V, labeller = label_parsed) +
        # scale_color_viridis_d(expression(italic(sigma[i])^2),
        #                       begin = 0.1, end = 0.85) +
        scale_shape_manual(values = c(1, 4), guide = FALSE) +
        scale_size_manual(values = c(1, 3), guide = FALSE) +
        scale_y_continuous("Proportion of species that survive",
                           breaks = c(0, 0.4, 0.8), limits = c(-0.01, 1)) +
        scale_x_continuous(expression(italic(d[1]) ~
                                          ("with" ~ italic(d[2]) == 0.1)),
                           breaks = c(-0.1, 0)) +
        guides(color = guide_legend(override.aes = list(shape = 19, size = 2))) +
        theme(strip.text = element_text(size = 10),
              strip.text.y = element_text(angle = 0),
              legend.position = "top",
              plot.title = element_text(hjust = 0.5,
                                        margin = margin(0,0,0,b=12)))

}

stoch_coexist_d1_ps <- map(c("sub-additive", "additive", "super-additive"),
                           stoch_coexist_d1_p_fun)
names(stoch_coexist_d1_ps) <- c("sub", "add", "super")

# stoch_coexist_d1_ps[["sub"]]
# stoch_coexist_d1_ps[["add"]]
# stoch_coexist_d1_ps[["super"]]

stoch_coexist_ps <- map(stoch_coexist_ps,
                        ~.x + theme(axis.title.y = element_blank(),
                                    strip.text.y = element_blank()))
stoch_coexist_d1_ps <- map(stoch_coexist_d1_ps,
                           ~.x + theme(axis.title.y = element_blank(),
                                       axis.text.y = element_blank()))
for (f in c("sub", "add")) {
    stoch_coexist_ps[[f]] <- stoch_coexist_ps[[f]] +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank())
    stoch_coexist_d1_ps[[f]] <- stoch_coexist_d1_ps[[f]] +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank())
}
for (f in c("add", "super")) {
    stoch_coexist_ps[[f]] <- stoch_coexist_ps[[f]] +
        theme(strip.text.x = element_blank())
    stoch_coexist_d1_ps[[f]] <- stoch_coexist_d1_ps[[f]] +
        theme(strip.text.x = element_blank())
}

stoch_coexist_p <- ggarrange(stoch_coexist_ps[["sub"]],
          stoch_coexist_d1_ps[["sub"]],
          stoch_coexist_ps[["add"]],
          stoch_coexist_d1_ps[["add"]],
          stoch_coexist_ps[["super"]],
          stoch_coexist_d1_ps[["super"]],
          ncol = 2, left = "Proportion of species that survive",
          draw = FALSE)


if (.RESAVE_PLOTS) {
    save_plot(stoch_coexist_p, 6.5, 7, .prefix = "S3-")
}







# =============================================================================*
# =============================================================================*

# Figs 5,S4 Stoch. - heatmaps ----

# =============================================================================*
# =============================================================================*

#'
#' For each set of parameter values, how many (out of 96 total) resulted in
#' the invading (`inv`) and the resident (`res`) species surviving?
#'

determ_inv_sims <- map_dfr(0:20,
                           function(i) {
                               fn <- rds(paste0("giant_inv_sims/giant_inv",
                                                "_sims_outcomes_", i))
                               readRDS(fn) %>%
                                   filter(sigma_V == 0, sigma_N == 0,
                                          eta != 0) %>%
                                   select(-sigma_V, -sigma_N)
                           }) %>%
    bind_rows(map_dfr(0:6,
                      function(i) {
                          fn <- rds(paste0("giant_inv_sims_add-mid/",
                                           "giant_inv_sims_add-mid_",
                                           "outcomes_", i))
                          readRDS(fn) %>%
                              filter(sigma_V == 0, sigma_N == 0) %>%
                              select(-sigma_V, -sigma_N)
                      })) %>%
    mutate(d1 = factor(d1, levels = sort(unique(d1)),
                       labels = sprintf("d[1] == %.4f", sort(unique(d1)))),
           eta = factor(eta, levels = c(-0.6, 0, 0.6),
                        labels = paste0(c("sub-", "", "super-"),
                                        "additive"))) %>%
    select(-total)

get_determ <- function(.par, .V1, .V2, .eta, .d1) {
    filter(determ_inv_sims, V1 == .V1[1], V2 == .V2[1],
           eta == .eta[1], d1 == .d1[1])[[.par]]
}

stoch_inv_sims <- map_dfr(0:20,
                          function(i) {
                              fn <- rds(paste0("giant_inv_sims/giant_inv",
                                               "_sims_outcomes_", i))
                              readRDS(fn) %>%
                                  filter((sigma_V == 0.1 & sigma_N == 0) |
                                             (sigma_V == 0 & sigma_N == 0.1)) %>%
                                  filter(eta != 0)
                          }) %>%
    bind_rows(map_dfr(0:6,
                      function(i) {
                          fn <- rds(paste0("giant_inv_sims_add-mid/",
                                           "giant_inv_sims_add-mid_",
                                           "outcomes_", i))
                          readRDS(fn) %>%
                              filter((sigma_V == 0.1 & sigma_N == 0) |
                                         (sigma_V == 0 & sigma_N == 0.1))
                      })) %>%
    mutate(d1 = factor(d1, levels = sort(unique(d1)),
                       labels = sprintf("d[1] == %.4f", sort(unique(d1)))),
           sigma_V = factor(sigma_V, levels = sort(unique(sigma_V)),
                            labels = sprintf("sigma[V] == %.4f",
                                             sort(unique(sigma_V)))),
           sigma_N = factor(sigma_N, levels = sort(unique(sigma_N)),
                            labels = sprintf("sigma[N] == %.4f",
                                             sort(unique(sigma_N)))),
           eta = factor(eta, levels = c(-0.6, 0, 0.6),
                            labels = paste0(c("sub-", "", "super-"),
                                            "additive"))) %>%
    mutate(across(coexist:extinct, ~ .x / total)) %>%
    select(-total) %>%
    group_by(V1, V2, eta, d1) %>%
    mutate(coexist_diff = coexist - get_determ("coexist", V1, V2, eta, d1),
           replace_diff = replace - get_determ("replace", V1, V2, eta, d1),
           reject_diff = reject - get_determ("reject", V1, V2, eta, d1),
           extinct_diff = extinct - get_determ("extinct", V1, V2, eta, d1)) %>%
    ungroup()






inv_sims_one_p_fun <- function(.par = "coexist",
                               .which_sigma = "V",
                               .fill_low = "#ca0020",
                               .fill_high = "#0571b0",
                               .fill_limits = c(-1, 1)) {

    # .par = "coexist"; .which_sigma = "V"
    # .fill_low = "#ca0020"; .fill_high = "#0571b0"; .fill_limits = c(-1, 1)
    # rm(.par, .which_sigma, .fill_low, .fill_high, .fill_limits)

    .par <- match.arg(.par, c("coexist", "replace", "reject", "extinct"))

    .which_sigma <- match.arg(.which_sigma, c("V", "N"))

    fill_scale <- scale_fill_gradient2(bquote("Effect of" ~
                                                  sigma[.(.which_sigma)]),
                                       low = .fill_low,
                                       mid = "white",
                                       high = .fill_high,
                                       midpoint = 0,
                                       limits = .fill_limits)

    make_eta_fct <- function(..x) {
        factor(..x, levels = c(-0.6, 0, 0.6),
               labels = paste0(c("sub-", "", "super-"),
                               "additive"))
    }
    midi <- ceiling(formals(stable_points)[["line_n"]] / 2)

    .sigma <- sprintf("sigma_%s", .which_sigma)

    # p <-
    stoch_inv_sims %>%
        filter(!!as.name(.sigma) == sprintf("sigma[%s] == %.4f",
                                            .which_sigma, 0.1)) %>%
        ggplot(aes(V1, V2)) +
        geom_raster(aes_string(fill = paste0(.par, "_diff"))) +
        geom_tile(data = determ_inv_sims %>% filter(!!as.name(.par) > 0),
                  fill = NA, color = "black", size = 0.1) +
        facet_grid(d1 ~ eta, label = label_parsed) +
        geom_point(data = map_dfr(c(-0.6, 0.6), ~ stable_points(.x) %>%
                                      mutate(eta = make_eta_fct(.x))),
                   aes(shape = factor(V1)),
                   size = 2, color = "black") +
        geom_path(data = stable_points(0) %>%
                      mutate(eta = make_eta_fct(0)),
                  color = "gray50", size = 1) +
        geom_point(data = tibble(V1 = sqrt(2), V2 = sqrt(2),
                                 eta = make_eta_fct(0)),
                   shape = 16, size = 2, color = "black") +
        scale_x_continuous("Axis 1\n(varying)", breaks = c(0, 2, 4)) +
        scale_y_continuous("(ameliorative)\nAxis 2", breaks = c(0, 2, 4)) +
        scale_shape_manual(values = c(16, 17, 16), guide = FALSE) +
        coord_equal() +
        NULL +
        theme(strip.text = element_text(size = 8),
              strip.text.y = element_text(angle = 0, hjust = 0,
                                          margin = margin(0,0,0,l=3)),
              strip.text.x = element_text(margin = margin(0,0,0,b=3)),
              panel.border = element_rect(size = 0.5, fill = NA),
              plot.title = element_text(size = 12, hjust = 0.5,
                                        margin = margin(0,0,0,b=6)),
              legend.title = element_text(size = 10),
              plot.margin = margin(0,0,0,0)) +
        fill_scale

}




inv_sims_ps <- list()
inv_sims_ps[["coexist"]] <- inv_sims_one_p_fun("coexist")
inv_sims_ps[["replace"]] <- inv_sims_one_p_fun("replace",
                                               .fill_low = "#e66101",
                                               .fill_high = "#5e3c99")
inv_sims_ps[["extinct"]] <- inv_sims_one_p_fun("extinct",
                                               .fill_low = "#d01c8b",
                                               .fill_high = "#4dac26",
                                               .fill_limits = NULL)

# inv_sims_ps[["coexist"]]
# inv_sims_ps[["replace"]]
# inv_sims_ps[["extinct"]]  # <-- super boring

# Also boring:
# inv_sims_one_p_fun("coexist", .which_sigma = "N")




if (.RESAVE_PLOTS) {
    save_plot(inv_sims_ps[["coexist"]], 5, 4.5, .name = "5-inv_sims_hm_coexist")
    save_plot(inv_sims_ps[["replace"]], 5, 4.5, .name = "S4-inv_sims_hm_replace")
}












# =============================================================================*
# =============================================================================*

# Additive --> Sub-additive ?? ----

# =============================================================================*
# =============================================================================*



# What direction does selection move them?

ss_t <- sauron:::sel_str_cpp
trnorm <- sauron:::trunc_rnorm_cpp

D <- matrix(0, 2, 2)
diag(D) <- c(-0.1, 0.1)
C <- diag(2)

V <- rbind(c(1.367653, 0.9889897),
           c(1.458830, 1.7362264))

N <- c(9.6430886, 0.9798823, 0.3046547)


dir_df <- crossing(.v1 = seq(0, 4, 0.2),
                   .v2 = seq(0, 4, 0.2)) %>%
    mutate(across(.fns = round, digits = 3)) %>%
    pmap_dfr(function(.v1, .v2) {
        ss <- ss_t(cbind(V, c(.v1, .v2)),
                   N,
                   formals(quant_gen)$f,
                   formals(quant_gen)$a0,
                   C,
                   formals(quant_gen)$r0,
                   D)
        tibble(V1 = .v1, V2 = .v2, V1_delta = ss[1,3], V2_delta = ss[2,3])
    })

dir_df %>%
    ggplot(aes(V1, V2)) +
    stable_points(0, return_geom = TRUE, color = "gray70") +
    geom_segment(aes(xend = V1 + 0.15 * V1_delta, yend = V2 + 0.15 * V2_delta),
                 arrow = arrow(length = unit(0.05, "inches"))) +
    scale_fill_viridis_c() +
    scale_x_continuous("Axis 1\n(conflicting)", breaks = c(0, 2, 4)) +
    scale_y_continuous("(ameliorative)\nAxis 2", breaks = c(0, 2, 4)) +
    coord_equal() +
    NULL


# set.seed(043506945)
maxt <- 2e4
.V_sigmas <- c(0.1, 0.2)

V_test <- matrix(0, maxt+1, 4)
V_test[1,1:2] <- c(0.4, sqrt(4-0.4^2))
# V_test[1,1:2] <- c(sqrt(4-0.4^2), 0.4)

for (k in 1:2) {
    if (.V_sigmas[k] > 0) {
        V_test[1,k+2] <- V_test[1,k] * rlnorm(1, - .V_sigmas[k]^2 / 2,
                                              .V_sigmas[k])
    } else V_test[1,k+2] <- V_test[1,k]
}


for (t in 1:maxt) {

    deltaV <- 0.05 * ss_t(cbind(V_test[t,3:4]),
                          pop_sizes(1, 0, c(-0.1, 0.1)),
                          formals(quant_gen)$f,
                          formals(quant_gen)$a0,
                          C,
                          formals(quant_gen)$r0,
                          D)

    V_test[t+1,1:2] <- pmax(0, V_test[t,1:2] + deltaV[,ncol(deltaV)])

    for (k in 1:2) {
        if (.V_sigmas[k] > 0) {
            V_test[t+1,k+2] <- V_test[t+1,k] * rlnorm(1, - .V_sigmas[k]^2 / 2,
                                                      .V_sigmas[k])
        } else V_test[t+1,k+2] <- V_test[t+1,k]
    }

}

colnames(V_test) <- c("V1", "V2", "Vp1", "Vp2")


# Standard deviation of difference between genotype and phenotype:
sqrt(sum((V_test[,"Vp1"] - V_test[,"V1"])^2) / nrow(V_test))
#' Method 1:
#'   - 0.6783879 for .sigma_V = 0.4
#'   - 0.2892102 for .sigma_V = 0.2
#'   - 0.1341219 for .sigma_V = 0.1
#' Method 2:
#'   - 0.411026 for .sigma_V = 0.4
#'   - 0.1908039 for .sigma_V = 0.2
#'   - 0.09974028 for .sigma_V = 0.1
#' Method 3:
#'   - 0.6127865 for .sigma_V = 0.4
#'   - 0.2876067 for .sigma_V = 0.2
#'   - 0.1336868 for .sigma_V = 0.1




V_test %>%
    # .[1:20e3L,] %>%
    # .[seq(1, nrow(.), 100),] %>%
    as_tibble() %>%
    ggplot(aes(V1, V2)) +
    geom_path(aes(Vp1, Vp2), color = "dodgerblue", alpha = 0.2) +
    geom_path() +
    geom_point(data = V_test[1,,drop=FALSE] %>% as_tibble(),
               shape = 1, size = 4, color = "firebrick") +
    geom_point(data = V_test[nrow(V_test),,drop=FALSE] %>% as_tibble(),
               shape = 4, size = 4, color = "firebrick") +
    stable_points(0, return_geom = TRUE, color = "gray40") +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray70") +
    coord_equal() #xlim = c(0, 3), ylim = c(0, 3))



V_test %>%
    # .[1:100,] %>%
    as_tibble() %>%
    mutate(angle = atan(Vp1 / Vp2) * 180 / pi,
           radius = sqrt(Vp1^2 + Vp2^2),
           time = 1:n()) %>%
    gather(type, value, angle:radius) %>%
    ggplot(aes(time, value)) +
    geom_line(color = "dodgerblue", alpha = 0.5) +
    geom_line(data = V_test %>%
                  # .[1:100,] %>%
                  as_tibble() %>%
                  mutate(angle = atan(V1 / V2) * 180 / pi,
                         radius = sqrt(V1^2 + V2^2),
                         time = 1:n()) %>%
                  gather(type, value, angle:radius),
              color = "black") +
    geom_hline(data = tibble(type = c("angle", "radius"),
                             yint = c(45, 2)),
               aes(yintercept = yint), color = "firebrick") +
    geom_hline(data = tibble(yint = c(0, 0)),
               aes(yintercept = yint), color = "gray70") +
    facet_wrap(~ type, scales = "free", ncol = 1) +
    theme(strip.text = element_text(size = 16),
          panel.spacing = unit(2, "lines"))


#'
#' The simulations above show that when stochasticity for evolution is
#' a normal distribution instead of lognormal, it doesn't inextricably move
#' toward the center point on the ring (i.e., it doesn't become sub-additive).
#'



V_test %>%
    as_tibble() %>%
    mutate(out = sqrt(Vp1^2 + Vp2^2) > 2,
           diff1 = V1 - lag(V1),
           diff2 = V2 - lag(V2)) %>%
    filter(!is.na(diff1)) %>%
    gather("axis", "value", diff1, diff2) %>%
    ggplot(aes(value, ..density..)) +
    geom_vline(xintercept = 0, linetype = 2, color = "gray70") +
    geom_freqpoly(aes(color = out), bins = 50) +
    facet_grid(~ axis, scales = "free") +
    scale_color_manual(values = c("firebrick", "dodgerblue"))




.V0 = "unrestricted"
.eta_sign = 0
.d1 = -0.1
.sigma_V = 0.1
.sigma_N = 0



which_switch <- 3 # which spp to switch for unrestricted

.n <- 5
.q = 2
stopifnot(is.numeric(.d1) && length(.d1) == 1)
.ds <- c(.d1, abs(.d1))
stopifnot(is.numeric(.eta_sign) && length(.eta_sign) == 1 &&
              .eta_sign %in% -1:1)
.V0 <- match.arg(.V0, c("restricted", "unrestricted"))
.lab <- .V0

.eta <- .eta_sign * etas[[2]]

if (.lab == "restricted" && .eta < 0) {
    stop(paste("\n'restricted' with sub-additivity is not programmed bc",
               "there's only one stable trait state.",
               "Thus, there is no way to restrict starting axes to",
               "be outside of that state's basin of attraction."))
}

if (.eta < 0) {
    .V0 <- rbind(seq(3, 2, length.out = 5),
                 seq(2, 3, length.out = 5))
} else if (.eta == 0) {
    .V0 <- rbind(seq(1.5, 0, length.out = 5),
                 seq(1.6, 3.1, length.out = 5))
    if (.lab == "unrestricted") {
        v2 <- .V0[2, which_switch]
        .V0[2, which_switch] <- .V0[1, which_switch]
        .V0[1, which_switch] <- v2
    }
} else {
    .V0 <- rbind(seq(1.2, 0, length.out = 5),
                 seq(1.3, 2.5, length.out = 5))
    if (.lab == "unrestricted") {
        v2 <- .V0[2, which_switch]
        .V0[2, which_switch] <- .V0[1, which_switch]
        .V0[1, which_switch] <- v2
    }
}
.V0 <- round(.V0, 3)

if (.sigma_V == 0 && .sigma_N == 0) {
    .nreps <- 1
    .N_THREADS <- 1
} else .nreps <- 12



set.seed(1415272918)
qg <- quant_gen(q = .q, eta = .eta, d = .ds,
                n_reps = .nreps, n = ncol(.V0),
                V0 = .V0,
                sigma_V0 = 0,
                sigma_V = .sigma_V,
                sigma_N = .sigma_N,
                spp_gap_t = 500L,
                final_t = 20e3L,
                save_every = 1L,
                add_var = rep(0.05, .n),
                n_threads = .N_THREADS,
                show_progress = FALSE)
qg0 <- quant_gen(q = .q, eta = .eta, d = .ds,
                n_reps = 1, n = ncol(.V0),
                V0 = .V0,
                sigma_V0 = 0,
                sigma_V = 0,
                sigma_N = 0,
                spp_gap_t = 500L,
                final_t = 20e3L,
                save_every = 1L,
                add_var = rep(0.05, .n),
                n_threads = 1,
                show_progress = FALSE)


nv <- qg %>%
    .[["nv"]] %>%
    mutate(trait = paste0("V", trait)) %>%
    nest(value = c(geno, pheno)) %>%
    spread(trait, value) %>%
    unnest(c(V1, V2), names_sep = "_") %>%
    rename(V1 = V1_geno,
           V2 = V2_geno,
           Vp1 = V1_pheno,
           Vp2 = V2_pheno)

nv %>%
    mutate(id = interaction(spp, rep)) %>%
    arrange(time) %>%
    # ggplot(aes(V1, V2, color = spp)) +
    ggplot(aes(Vp1, Vp2, color = spp)) +
    geom_path(aes(group = id)) +
    geom_point(data = nv %>%
                   group_by(spp) %>%
                   summarize(V1 = V1[time == min(time)][1],
                             V2 = V2[time == min(time)][1])) +
    facet_wrap(~ rep, nrow = 3) +
    coord_equal(xlim = c(0, 3.11), ylim = c(0, 3.11)) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    ylab("(ameliorative)\nAxis 2") +
    xlab("Axis 1\n(conflicting)")

nv %>%
    mutate(id = interaction(spp, rep)) %>%
    ggplot(aes(time / 1000, N, color = spp)) +
    geom_path(aes(group = id)) +
    facet_wrap(~ rep, nrow = 3) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    ylab("Abundance") +
    xlab("Time (× 1000)")
