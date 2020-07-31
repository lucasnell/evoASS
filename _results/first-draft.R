
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

set.seed(1898348146)
etas <- map(1:6, ~ with(list(.q = .x), runif(.q * ((.q - 1) / 2), 0.1, 0.4)))
# With just one eta, it can be a simple number:
etas[[2]] <- 0.6


# =============================================================================*
# =============================================================================*

# 1. 2-trait outcomes ----

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
#' Based on eigenvalues (see `first-draft__stability.R`)...
#'   * When the tradeoff is additive, the state is neutrally stable
#'   * Everything else is stable
#'
#'




outcomes_q2_p <- eta_sim_df %>%
    mutate(eta1 = factor(eta1, levels = sort(unique(eta1)),
                         labels = c("sub-additive", "additive",
                                    "super-additive"))) %>%
    group_by(eta1) %>%
    filter(unq_spp_filter(V1, V2, .prec = 0.05)) %>%
    ungroup() %>%
    ggplot(aes(V1, V2)) +
    geom_point(shape = 21, color = "dodgerblue", fill = "dodgerblue",
               alpha = 0.5, size = 3) +
    scale_x_continuous("Trait 1", breaks = 0:2) +
    scale_y_continuous("Trait 2", breaks = 0:2) +
    coord_equal(xlim = c(0, 2.5), ylim = c(0, 2.5)) +
    facet_wrap(~ eta1, nrow = 1) +
    theme(strip.text = element_text(size = 10),
          panel.border = element_rect(size = 0.5, fill = NA)) +
    NULL



if (.RESAVE_PLOTS) save_plot(outcomes_q2_p, 5, 2, .prefix = "1-")





# =============================================================================*
# =============================================================================*

# 2. Coexistence ----

# =============================================================================*
# =============================================================================*



one_coexist_combo <- function(.d1, .eta, .add_var, .pb = NULL, .vary_d2 = FALSE) {

    if (!is.null(.pb)) .pb$tick(0)

    # .d1 = 0.15; .eta = etas[[2]]; .add_var = 0.01
    # rm(.d1, .eta, .add_var)

    if (.vary_d2) {
        .d <- c(.d1, .d1)
    } else .d <- c(.d1, 0.1)
    .n <- 100

    .seed <- sample.int(2^31 - 1, 1)
    set.seed(.seed)

    Z <- quant_gen(q = 2, eta = .eta, d = .d, n = .n,
                   add_var = rep(.add_var, .n),
                   spp_gap_t = 500L, final_t = 50e3L, n_reps = 12,
                   save_every = 0L,
                   n_threads = .N_THREADS, show_progress = FALSE)

    Z$call[["eta"]] <- eval(.eta)
    Z$call[["d"]] <- eval(.d)
    Z$call[["n"]] <- eval(.n)
    Z$call[["add_var"]] <- eval(rep(.add_var, .n))

    jacs <- jacobians(Z)

    if (!is.null(.pb)) .pb$tick()

    return(list(NV = Z$nv %>%
                    mutate(trait = paste0("V", trait)) %>%
                    spread(trait, geno) %>%
                    mutate(d1 = .d1, eta = .eta, add_var = .add_var),
                J = jacs,
                seed = .seed))
}



# ------------------------------------*
# __ vary d1 & d2 ----
# ------------------------------------*

# Simulations varying both d values
if (.REDO_SIMS) {
    # Takes ~61 min w/ 3 threads
    coexist_sims <- crossing(.d1 = seq(-0.25, 2, length.out = 10),
                             .eta = c(-1,1) * etas[[2]],
                             .add_var = seq(0.01, 0.1, 0.01)) %>%
        # bc `seq` makes weird numbers that are ~1e-15 from what they should be:
        mutate(across(.fns = round, digits = 2))
    pb <- progress_bar$new(format = " [:bar] :percent in :elapsed | eta: :eta",
                           total = nrow(coexist_sims), clear = FALSE,
                           show_after = 0)
    set.seed(1558743552)
    coexist_sims <- coexist_sims %>%
        pmap(one_coexist_combo, .pb = pb, .vary_d2 = TRUE)
    # Because we're not only varying d1:
    for (i in 1:length(coexist_sims)) {
        coexist_sims[[i]][["NV"]] <- rename(coexist_sims[[i]][["NV"]], d = d1)
    }
    saveRDS(coexist_sims, rds("coexist_sims"))
} else {
    coexist_sims <- readRDS(rds("coexist_sims"))
}



coexist_spp_df <- map_dfr(coexist_sims, ~ .x[["NV"]]) %>%
    group_by(d, eta, add_var, rep) %>%
    summarize(n_spp = n()) %>%
    ungroup() %>%
    mutate(eta = factor(eta, levels = c(-1, 1) * etas[[2]],
                        labels = c("sub-additive", "super-additive")))

coexist_spp_p1 <- coexist_spp_df %>%
    mutate(n_spp = n_spp / 100) %>%
    ggplot(aes(add_var, n_spp, color = d)) +
    geom_jitter(width = 0.001, height = 0, shape = 1, size = 0.5) +
    facet_wrap(~ eta, ncol = 1) +
    scale_color_viridis_c(expression(d[1] ~ "and" ~ d[2]),
                          begin = 0.1, end = 0.85, option = "A",
                          breaks = c(0, 1, 2)) +
    scale_y_continuous("Proportion of species that survive",
                       breaks = c(0, 0.4, 0.8), limits = c(0, 1)) +
    xlab(expression(sigma^2)) +
    theme(legend.position = c(0.95, 0.55),
          legend.direction = "horizontal",
          legend.justification = c(1, 0),
          legend.key.height = unit(6, "pt"),
          legend.key.width = unit(12, "pt"),
          legend.title = element_text(face = "bold.italic")) +
    guides(color = guide_colorbar(title.position="top", title.hjust = 0.5,
                                  title.vjust = -1, label.vjust = 2))



coexist_spp_p2 <- coexist_spp_df %>%
    mutate(n_spp = n_spp / 100) %>%
    ggplot(aes(d, n_spp, color = add_var)) +
    geom_jitter(width = 0.02, height = 0, shape = 1, size = 0.5) +
    facet_wrap(~ eta, ncol = 1) +
    scale_color_viridis_c(expression(italic(sigma)^2), begin = 0.1, end = 0.85,
                          breaks = c(0.05, 0.1)) +
    scale_y_continuous("Proportion of species that survive",
                       breaks = c(0, 0.4, 0.8), limits = c(0, 1)) +
    xlab(expression(italic(d[1]) ~ "and" ~ italic(d[2]))) +
    theme(legend.position = c(0.95, 0.55),
          legend.direction = "horizontal",
          legend.justification = c(1, 0),
          legend.key.height = unit(6, "pt"),
          legend.key.width = unit(12, "pt"),
          legend.title = element_text(face = "italic")) +
    guides(color = guide_colorbar(title.position="top", title.hjust = 0.5,
                                  title.vjust = -1, label.vjust = 2))


#'
#' Code in `first-draft__stability.R` shows that...
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




# Simulations varying just d1
if (.REDO_SIMS) {
    # Takes ~2.2 min w/ 3 threads
    vary_d1_coexist_sims <- crossing(.d1 = seq(-0.15, 0.05, 0.05),
                                     .eta = c(-1,1) * etas[[2]],
                                     .add_var = 0.1) %>%
        # bc `seq` makes weird numbers that are ~1e-15 from what they should be:
        mutate(across(.fns = round, digits = 2))
    pb <- progress_bar$new(format = " [:bar] :percent in :elapsed | eta: :eta",
                           total = nrow(vary_d1_coexist_sims), clear = FALSE,
                           show_after = 0)
    set.seed(1472331523)
    vary_d1_coexist_sims <- vary_d1_coexist_sims %>%
        pmap(one_coexist_combo, .pb = pb)
    saveRDS(vary_d1_coexist_sims, rds("vary_d1_coexist_sims"))
} else {
    vary_d1_coexist_sims <- readRDS(rds("vary_d1_coexist_sims"))
}

vary_d1_coexist_spp_df <- map_dfr(vary_d1_coexist_sims, ~ .x[["NV"]]) %>%
    group_by(d1, eta, add_var, rep) %>%
    summarize(n_spp = n()) %>%
    ungroup() %>%
    mutate(eta = factor(eta, levels = c(-1, 1) * etas[[2]],
                        labels = c("sub-additive", "super-additive")))


coexist_spp_p3 <- vary_d1_coexist_spp_df %>%
    mutate(n_spp = n_spp / 100) %>%
    ggplot(aes(d1, n_spp)) +
    geom_vline(data = tibble(xint = c(-0.1, 0), col = factor(1:2)),
               aes(xintercept = xint, color = col), linetype = 2) +
    geom_jitter(width = 0.002, height = 0, shape = 1, size = 1) +
    facet_wrap(~ eta, ncol = 1) +
    scale_color_manual(NULL, values = c("firebrick", "dodgerblue"), guide = FALSE) +
    scale_y_continuous("Proportion of species that survive",
                       breaks = c(0, 0.4, 0.8), limits = c(0, 1)) +
    xlab(expression(italic(d[1]) ~ ("with" ~ italic(d[2]) == 0.1)))







coexist_spp_p <- plot_grid(coexist_spp_p1,
                           coexist_spp_p2 +
                               theme(axis.title.y = element_blank()),
                           coexist_spp_p3 +
                               theme(axis.title.y = element_blank()),
                           nrow = 1, labels = LETTERS[1:3], align = "vh",
                           label_fontface = "plain")

coexist_spp_p


if (.RESAVE_PLOTS) save_plot(coexist_spp_p, 6.5, 4, .prefix = "2-")








#' # ------------------------------------*
#' # __ vary d1 & d2; V0 ~ U ----
#' # ------------------------------------*
#'
#'
#'
#' #'
#' #' NOTE: this function always changes both `d` values
#' #'
#' one_coexist_U_combo <- function(.d, .eta, .add_var, .pb = NULL) {
#'
#'     if (!is.null(.pb)) .pb$tick(0)
#'
#'     # .d = 0.15; .eta = etas[[2]]; .add_var = 0.01
#'     # rm(.d, .eta, .add_var)
#'
#'     stopifnot(is.numeric(.d) && length(.d) == 1)
#'
#'     .d <- rep(.d, 2)
#'     .n <- 100
#'
#'     .V0 <- matrix(runif(.n * 2, 0, 3), nrow = .n)
#'
#'     .seed <- sample.int(2^31 - 1, 1)
#'     set.seed(.seed)
#'
#'     Z <- quant_gen(q = 2, eta = .eta, d = .d, n = .n,
#'                    V0 = .V0, sigma_V0 = 0,
#'                    add_var = rep(.add_var, .n),
#'                    spp_gap_t = 500L, final_t = 50e3L, n_reps = 12,
#'                    save_every = 0L,
#'                    n_threads = .N_THREADS, show_progress = FALSE)
#'
#'     Z$call[["eta"]] <- eval(.eta)
#'     Z$call[["d"]] <- eval(.d)
#'     Z$call[["n"]] <- eval(.n)
#'     Z$call[["add_var"]] <- eval(rep(.add_var, .n))
#'
#'     jacs <- jacobians(Z)
#'
#'     if (!is.null(.pb)) .pb$tick()
#'
#'     return(list(NV = Z$nv %>%
#'                     mutate(trait = paste0("V", trait)) %>%
#'                     spread(trait, geno) %>%
#'                     mutate(d = .d[1], eta = .eta, add_var = .add_var),
#'                 J = jacs,
#'                 V0 = .V0,
#'                 seed = .seed))
#' }
#'
#'
#' # Simulations varying both d values
#' if (.REDO_SIMS) {
#'     # Takes ~XX min w/ 3 threads
#'     t0 <- Sys.time()
#'     coexist_U_sims <- crossing(.d = seq(-0.25, 2, length.out = 10),
#'                              .eta = c(-1,1) * etas[[2]],
#'                              .add_var = seq(0.01, 0.1, 0.01)) %>%
#'         # bc `seq` makes weird numbers that are ~1e-15 from what they should be:
#'         mutate(across(.fns = round, digits = 2))
#'     pb <- progress_bar$new(format = " [:bar] :percent in :elapsed | eta: :eta",
#'                            total = nrow(coexist_U_sims), clear = FALSE,
#'                            show_after = 0)
#'     set.seed(2014151565)
#'     coexist_U_sims <- coexist_U_sims %>%
#'         pmap(one_coexist_U_combo, .pb = pb)
#'     saveRDS(coexist_U_sims, rds("coexist_U_sims"))
#'     t1 <- Sys.time()
#'     t1 - t0; rm(t0, t1)
#' } else {
#'     coexist_U_sims <- readRDS(rds("coexist_U_sims"))
#' }








# =============================================================================*
# =============================================================================*

# 3,S1. Conditional coexistence ----

# =============================================================================*
# =============================================================================*




cond_coexist_test <- function(.V0, .eta_sign, .d1) {

    # .V0 = "restricted"; .eta_sign = 1; .d1 = -0.1
    # rm(.V0, .eta_sign, .d1)

    .dist <- 0.6
    which_switch <- 3 # which spp to switch for unrestricted

    .n <- 5
    .q = 2
    stopifnot(is.numeric(.d1) && length(.d1) == 1)
    .ds <- c(.d1, abs(.d1))
    stopifnot(is.numeric(.eta_sign) && length(.eta_sign) == 1 &&
                  .eta_sign %in% c(-1,1))
    .V0 <- match.arg(.V0, c("restricted", "unrestricted"))
    .lab <- .V0

    .eta <- .eta_sign * etas[[2]]

    if (.lab == "restricted" && .eta < 0) {
        stop(paste("\n'restricted' with sub-additivity is not programmed bc",
                   "there's only one stable trait state.",
                   "Thus, there is no way to restrict starting traits to",
                   "be outside of that state's basin of attraction."))
    }

    if (.eta < 0) {
        .V0 <- rbind(seq(2, 3, length.out = 5),
                     seq(3, 2, length.out = 5))
    } else {
        .V0 <- rbind(seq(1.2, 0, length.out = 5),
                     seq(1.3, 2.5, length.out = 5))
        if (.lab == "unrestricted") {
            v2 <- .V0[2, which_switch]
            .V0[2, which_switch] <- .V0[1, which_switch]
            .V0[1, which_switch] <- v2
        }
    }


    qg <- quant_gen(q = .q, eta = .eta, d = .ds,
                   n_reps = 1, n = ncol(.V0),
                   V0 = .V0,
                   sigma_V0 = 0,
                   spp_gap_t = 500L,
                   final_t = 20e3L,
                   add_var = rep(0.05, .n),
                   show_progress = FALSE)

    qg$call[["q"]] <- eval(.q)
    qg$call[["n"]] <- eval(.n)
    qg$call[["d"]] <- eval(.ds)
    qg$call[["eta"]] <- eval(.eta)
    qg$call[["add_var"]] <- eval(rep(0.05, .n))

    out <- list(nv = qg %>%
                    .[["nv"]] %>%
                    mutate(trait = paste0("V", trait)) %>%
                    spread(trait, geno) %>%
                    select(-rep) %>%
                    mutate(V0 = .lab, eta = .eta, d1 = .d1),
                jacs = jacobians(qg))

    return(out)
}


if (.REDO_SIMS) {
    # Takes just a few seconds
    cond_coexist <- crossing(.V0 = c("restricted", "unrestricted"),
                             .eta_sign = c(-1,1),
                             .d1 = c(-0.1, 0.1)) %>%
        filter(!(.V0 == "restricted" & .eta_sign < 0)) %>%
        pmap(cond_coexist_test)
    saveRDS(cond_coexist, rds("cond_coexist"))
} else {
    cond_coexist <- readRDS(rds("cond_coexist"))
}




cond_coexist_df <- cond_coexist %>%
    map_dfr(~ .x$nv) %>%
    filter(time < 4e3L) %>%
    mutate(V0 = factor(V0, levels = c("restricted", "unrestricted")),
           eta = factor(eta, levels = c(-1, 1) * etas[[2]],
                        labels = c("sub-additive", "super-additive")),
           d1 = factor(d1, levels = c(-0.1, 0.1),
                       labels = c("conflicting",
                                  "non-conflicting"))) %>%
    # `trait_space` is a combination of starting trait values and eta:
    mutate(trait_space =
               case_when(V0 == "restricted" & eta == "super-additive" ~ "i",
                         V0 == "unrestricted" & eta == "super-additive" ~ "ii",
                         V0 == "unrestricted" & eta == "sub-additive" ~ "iii",
                         TRUE ~ NA_character_) %>%
               factor(levels = c("i", "ii", "iii")),
           trt = sprintf("'%s' ~ {%s}", paste(trait_space),
                         ifelse(d1 == "conflicting", "-{}", "+{}")) %>%
               factor(levels = sprintf("'%s' ~ {%s}",
                                       rep(levels(trait_space), 2),
                                       rep(c("-{}", "+{}"), each = 3))))

if (any(is.na(cond_coexist_df$trait_space))) {
    stop("\nERROR: unknown combination of V0 and eta")
}


# Starting conditions:

cond_coexist_sc_p <- cond_coexist_df %>%
    filter(d1 == "non-conflicting") %>%
    group_by(trait_space, spp) %>%
    filter(time == min(time)) %>%
    ungroup() %>%
    select(trait_space, spp, V1, V2) %>%
    ggplot(aes(V1, V2)) +
    geom_abline(data = tibble(trait_space = factor(c("i", "ii"),
                                                   levels = c("i", "ii",
                                                              "iii")),
                              slp = 1, int = 0),
                aes(slope = slp, intercept = int), linetype = 2, color = "gray70") +
    geom_point(aes(color = spp), size = 2) +
    geom_point(data = map_dfr(etas[[2]] * c(1, 1, -1), ~ stable_points(.x)) %>%
                   mutate(trait_space = map2(c("i", "ii", "iii"), c(2,2,1), rep) %>%
                              do.call(what = c) %>%
                              factor(levels = c("i", "ii", "iii")),
                          shp = case_when(V1 > 0 & V2 > 0 ~ 3L,
                                          V1 == 0 ~ 2L,
                                          TRUE ~ 1L) %>%
                              factor(levels = 1:3)),
               aes(shape = shp), size = 4, color = "black") +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    scale_shape_manual(values = 0:2, guide = FALSE) +
    scale_size_continuous(range = c(0.1, 1)) +
    facet_grid(~ trait_space) +
    coord_equal(xlim = c(0, 3), ylim = c(0, 3)) +
    ylab("Trait 2") +
    xlab("Trait 1") +
    theme(plot.margin = margin(0,0,0,b=6))


# Time series of abundances
cc_N_p_list <- map(c("non-conflicting", "conflicting"),
                   ~ cond_coexist_df %>%
                       filter(d1 == .x) %>%
                       ggplot(aes(time / 1000L, N, color = spp)) +
                       geom_line() +
                       facet_wrap(~ trait_space, nrow = 1) +
                       scale_color_brewer(palette = "Dark2",
                                          guide = FALSE) +
                       scale_y_continuous("Abundance", trans = "log",
                                          breaks = 10^(c(-3, 0, 3)),
                                          labels = parse(
                                              text = sprintf("10^{%i}",
                                                             c(-3, 0, 3)))) +
                       xlab("Time (× 1,000)") +
                       theme(plot.margin = margin(0,0,t=10,b=10),
                             axis.title.x = element_blank())
)


names(cc_N_p_list) <- c("non-conflicting", "conflicting")



stable_state_df <- map_dfr(1:3,
                           function(i) {
                               ts <- cond_coexist_df$trait_space %>%
                                   unique() %>%
                                   sort() %>%
                                   .[i]
                               stable_points((etas[[2]] * c(1,1,-1))[i]) %>%
                                   mutate(trait_space = ts)
                           }) %>%
    mutate(time = max(cond_coexist_df$time),
           shp = case_when(V1 > 0 & V2 > 0 ~ 3L,
                           V1 == 0 ~ 2L,
                           TRUE ~ 1L) %>%
               factor(levels = 1:3)) %>%
    gather(trait, value, V1:V2) %>%
    mutate(trait = gsub("^V", "trait ", trait))


cc_V_p_list <- map(c("non-conflicting", "conflicting"),
                   ~ cond_coexist_df %>%
                       filter(d1 == .x) %>%
                       gather(trait, value, V1:V2) %>%
                       mutate(trait = gsub("^V", "trait ", trait)) %>%
                       ggplot(aes(time / 1000L, value)) +
                       geom_hline(yintercept = 0, size = 0.5,
                                  linetype = 1, color = "gray70") +
                       geom_vline(xintercept = 0, size = 0.5,
                                  linetype = 1, color = "gray70") +
                       geom_point(data = stable_state_df,
                                  aes(shape = shp), size = 3) +
                       geom_line(aes(color = spp)) +
                       facet_grid(trait ~ trait_space) +
                       scale_color_brewer(palette = "Dark2", guide = FALSE) +
                       scale_shape_manual(values = 0:2, guide = FALSE) +
                       scale_y_continuous("Trait value", limits = c(-0.2, NA)) +
                       xlab("Time (× 1,000)") +
                       theme(plot.margin = margin(0,0,0,t=10))
)

names(cc_V_p_list) <- c("non-conflicting", "conflicting")



cond_coexist_ps <- map(c("non-conflicting", "conflicting"),
                       ~ plot_grid(cond_coexist_sc_p,
                                   plot_grid(cc_N_p_list[[.x]],
                                             cc_V_p_list[[.x]],
                                             ncol = 1, rel_heights = c(1, 2),
                                             align = "v", axis = "lr",
                                             labels = LETTERS[2:3],
                                             label_x = 0.06),
                                   ncol = 1, rel_heights = c(1, 2),
                                   labels = c(LETTERS[1], ""), label_x = 0.06))
names(cond_coexist_ps) <- c("non-conflicting", "conflicting")




if (.RESAVE_PLOTS) {
    save_plot(cond_coexist_ps[["non-conflicting"]], 6.5, 6.5,
              .name = "S1-cond_coexist_non-conflicting")
    save_plot(cond_coexist_ps[["conflicting"]], 6.5, 6.5,
              .name = "3-cond_coexist_conflicting")
}



#'
#' They're all stable except for sub-additive and conflicting trait 1,
#' which is neutrally stable (see `first-draft__stability.R`).
#'













# ===========================================================================*
# ===========================================================================*

# S2,S3. Invasibility ----

# ===========================================================================*
# ===========================================================================*

#'
#' This shouldn't be a figure in the main text, but it shows
#' the invasibility of resulting communities for these situations:
#'
#' - Conflicting evolution (d = -0.01)
#' - Neutral evolution (d = 0)
#' - Non-conflicting evolution (d = +0.01)
#'
#' ... for eta = -0.6 and +0.6
#'
#' For these simulations, q = 2
#'




one_d_invasion <- function(.x) {

    .d <- .x[[".d"]][[1]]
    .n <- .x[[".n"]][[1]]
    .inv_N0 <- .x[[".inv_N0"]][[1]]
    .eta <- .x[[".eta_sign"]][[1]] * etas[[2]]

    V0 <- stable_points(.eta) %>%
        as.matrix() %>%
        split(1:nrow(.)) %>%
        rep(ceiling(.n / length(.))) %>%
        .[1:.n] %>%
        do.call(what = cbind) %>%
        {colnames(.) <- NULL; .}
    N0 <- pop_sizes(.n, .eta, .d)

    v1h <- stable_points(.eta)[["V1"]][1]
    v2h <- stable_points(.eta)[["V2"]][1]

    inv_df <- map_dfr(1:8 * 0.1,
                      function(.dist) {
                          rbind(c(v1h, v2h + .dist),
                                c(v1h, v2h - .dist),
                                c(v1h + .dist, v2h),
                                c(v1h - .dist, v2h),
                                c(v1h + .dist / sqrt(2), v2h + .dist / sqrt(2)),
                                c(v1h - .dist / sqrt(2), v2h - .dist / sqrt(2)),
                                c(v1h + .dist / sqrt(2), v2h - .dist / sqrt(2)),
                                c(v1h - .dist / sqrt(2), v2h + .dist / sqrt(2))) %>%
                              {colnames(.) <- c(".v1", ".v2"); .} %>%
                              as_tibble()
                      }) %>%
        add_row(.v1 = v1h, .v2 = v2h) %>%
        filter(across(.fns = ~ . >= 0)) %>%
        pmap_dfr(function(.v1, .v2) {
            .inv <- quant_gen(q = 2, eta = .eta, d = .d, n = .n + 1,
                              save_every = 0L, final_t = 50e3L, n_reps = 1,
                              N0 = c(N0, N0[1] * .inv_N0),
                              V0 = cbind(V0, rbind(.v1, .v2)),
                              spp_gap_t = 0, sigma_V0 = 0,
                              show_progress = FALSE) %>%
                .[["nv"]] %>%
                filter(trait == 1, spp == .n+1) %>%
                nrow() %>%
                `>=`(1) %>%
                as.integer()
            tibble(dist = sqrt((.v1 - v1h)^2 + (.v2 - v2h)^2),
                   V1 = .v1, V2 = .v2,
                   invaded = .inv)
        }) %>%
        mutate(d = .d, res_n = .n, inv_N0 = .inv_N0, eta = .eta)

    return(inv_df)

}









if (.REDO_SIMS) {
    # Takes ~ 3 min w/ 3 threads
    invade_df <- crossing(.d = 0.01 * c(-1, 0, 1),
                          .n = c(2, 10),
                          .inv_N0 = c(0.01, 0.1, 1),
                          .eta_sign = c(-1, 1)) %>%
        split(1:nrow(.)) %>%
        mclapply(FUN = one_d_invasion, mc.cores = .N_THREADS)
    invade_df <- bind_rows(invade_df)
    saveRDS(invade_df, rds("invade_sims"))
} else {
    invade_df <- readRDS(rds("invade_sims"))
}



# invasion_caption <- "Successful invasion of an equilibrium community based on
#                      the invader's starting distance (in trait space)
#                      from the stable point.
#                      Columns of sub-panels separate whether evolution was
#                      conflicting, neutral, or non-conflicting.
#                      Rows of sub-panels separate the invaders' starting
#                      abundances ($N_{eq}$) in relation to the residents'
#                      abundances ($N_{res}$); all residents had the same
#                      abundance.
#                      The point color indicates the number of species present
#                      in the resident community."



invasion_p_fun <- function(.eta) {
    invade_df %>%
        filter(eta == .eta) %>%
        mutate(invaded = factor(invaded, levels = 0:1, labels = c("no", "yes")),
               d = factor(d, levels = 0.01 * c(-1, 0, 1),
                          labels = c("'conflicting'", "'neutral'",
                                     "'non-conflicting'")),
               id = sprintf("%s (%i spp)", invaded, res_n) %>%
                   factor(levels = rev(c("yes (2 spp)", "no (2 spp)",
                                         "yes (10 spp)", "no (10 spp)"))),
               inv_N0 = factor(inv_N0, levels = c(0.01, 0.1, 1),
                               labels = sprintf("N[inv] == %s", c("frac(N[res],100)",
                                                                  "frac(N[res],10)",
                                                                  "N[res]"))),
               res_n = factor(res_n, levels = c(2, 10),
                              labels = sprintf("%i species", c(2, 10)))) %>%
        ggplot(aes(id, dist)) +
        ggtitle(ifelse(.eta > 0, "super-additive", "sub-additive")) +
        geom_vline(dat = tibble(xint = c("yes (2 spp)", "no (2 spp)",
                                         "yes (10 spp)", "no (10 spp)")) %>%
                       mutate(xint = factor(xint, levels = rev(xint)),
                              lty = rep(1:2, 2) %>% factor()),
                   aes(xintercept = xint, linetype = lty),
                   size = 0.5, color = "gray70") +
        geom_jitter(aes(color = res_n, fill = res_n),
                    height = 0, width = 0.25, shape = 21) +
        facet_grid(inv_N0 ~ d, label = label_parsed) +
        xlab("Successful invasion") +
        ylab("Distance from stable point") +
        scale_color_viridis_d("starting community size:", begin = 0.3,
                              end = 0.7, option = "A", guide = FALSE) +
        scale_fill_viridis_d("starting community size:", begin = 0.3,
                             end = 0.7, option = "A", alpha = 0.5, guide = FALSE) +
        scale_linetype_manual(values = 1:2, guide = FALSE) +
        theme(panel.spacing = unit(2, "lines"),
              strip.text.y = element_text(angle = 0, margin = margin(0,0,0,l=6)),
              legend.position = "top",
              legend.title = element_text(size = 10),
              plot.title = element_text(hjust = 0.5,
                                        margin = margin(0,0,0,b=12))) +
        coord_flip() +
        NULL
}



invasion_ps <- map(c(-0.6, 0.6), invasion_p_fun)
names(invasion_ps) <- c("sub", "super")



if (.RESAVE_PLOTS) {
    save_plot(invasion_ps[["sub"]], 6, 5, .name = "S2-invasion_sub")
    save_plot(invasion_ps[["super"]], 6, 5, .name = "S3-invasion_super")
}








# ===========================================================================*
# ===========================================================================*

# Supp. info: "Filling in" of trait space ----

# ===========================================================================*
# ===========================================================================*

#'
#' This shouldn't be a figure or anything (or at least not one in
#' the main text), but it quickly shows that there is no tendency for
#' species to "disperse" throughout the trait space when tradeoffs
#' are additive.
#'

ring_points <- stable_points(0, line_n = 90)

fill_qg <- quant_gen(q = 2, eta = 0, d = 1e-2, max_t = 20e3L, n_reps = 1,
                     save_every = 0L, n = 2,
                     N0 = rep(5.568996746, 2),
                     V0 = ring_points %>%
                         .[c(30, 60),] %>%
                         as.matrix() %>%
                         split(1:nrow(.)),
                     start_t = 0, sigma_V0 = 0, show_progress = FALSE)

print_big_nums(ring_points[c(30, 60),])
print_big_nums(fill_qg$nv)


# Takes ~ 6 sec w/ 3 threads
filling_sims <- lapply(1:nrow(ring_points),
                         function(i) {
                             .new_V <- ring_points %>%
                                 .[i,] %>%
                                 as.matrix() %>%
                                 t() %>%
                                 list()
                             perturb(fill_qg, max_t = 20e3L, save_every = 0L,
                                     new_V = .new_V, new_N = 1) %>%
                                 .[["end"]] %>%
                                 .[["N"]]
                         })

map_dbl(filling_sims, ~ .x[3]) %>% range()


pop_sizes(3, 0, 1e-2)








# STOCHASTICITY AND COEXISTENCE ----


set.seed(4235668)
X <- quant_gen(eta = 0.6, d = c(-0.1, 0.1), q = 2, n = 10,
               spp_gap_t = 100L, final_t = 1e3L, save_every = 1L,
               # sigma_N = 0.1,
               sigma_V = 0.1,
               n_reps = 10)
set.seed(4235669)
Y <- quant_gen(eta = 0.6, d = -0.1, q = 2, n = 10,
               spp_gap_t = 100L, final_t = 1e3L, save_every = 1L,
               # sigma_N = 0.1,
               sigma_V = 0.1,
               n_reps = 10)
set.seed(4235670)
Z <- quant_gen(eta = 0.6, d = 0.1, q = 2, n = 10,
               spp_gap_t = 100L, final_t = 1e3L, save_every = 1L,
               # sigma_N = 0.1,
               sigma_V = 0.1,
               n_reps = 10)

Nts <- function(.QG, .title = NULL) {
    .QG$nv %>%
        filter(trait == 1) %>%
        ggplot(aes(time, N)) +
        geom_line(aes(color = spp), na.rm = TRUE) +
        ggtitle(.title) +
        facet_wrap(~ rep) +
        scale_color_viridis_d(begin = 0.1, end = 0.9, option = "C", guide = FALSE)
}

Vts <- function(.QG, .title = NULL) {
    .QG$nv %>%
        mutate(trait = paste0("V", trait)) %>%
        select(rep, time, spp, trait, geno) %>%
        spread(trait, geno) %>%
        arrange(rep, spp, time) %>%
        ggplot(aes(V1, V2)) +
        geom_path(aes(color = spp), na.rm = TRUE) +
        stable_points(.QG$call[["eta"]], return_geom = TRUE,
                      shape = 1, size = 4, color = "black") +
        geom_point(data = .QG$nv %>%
                       filter(time == max(time)) %>%
                       mutate(trait = paste0("V", trait)) %>%
                       select(rep, time, spp, trait, geno) %>%
                       spread(trait, geno) %>%
                       group_by(rep) %>%
                       filter(unq_spp_filter(V1, V2)) %>%
                       ungroup(),
                   shape = 4, size = 4, color = "black") +
        ggtitle(.title) +
        facet_wrap(~ rep) +
        scale_color_viridis_d(begin = 0.1, end = 0.9, option = "C", guide = FALSE)
}
Vpts <- function(.QG, .title = NULL) {
    stopifnot("pheno" %in% colnames(.QG$nv))
    .QG$nv %>%
        mutate(trait = paste0("V", trait)) %>%
        select(rep, time, spp, trait, pheno) %>%
        spread(trait, pheno) %>%
        arrange(rep, spp, time) %>%
        ggplot(aes(V1, V2)) +
        geom_path(aes(color = spp), na.rm = TRUE) +
        stable_points(.QG$call[["eta"]], return_geom = TRUE,
                      shape = 1, size = 4, color = "black") +
        geom_point(data = .QG$nv %>%
                       filter(time == max(time)) %>%
                       mutate(trait = paste0("V", trait)) %>%
                       select(rep, time, spp, trait, pheno) %>%
                       spread(trait, pheno) %>%
                       group_by(rep) %>%
                       filter(unq_spp_filter(V1, V2)) %>%
                       ungroup(),
                   shape = 4, size = 4, color = "black") +
        ggtitle(.title) +
        facet_wrap(~ rep) +
        scale_color_viridis_d(begin = 0.1, end = 0.9, option = "C", guide = FALSE)
}

Nts(X, expression(- {} + {}))
Nts(Y, expression(- {} - {}))
Nts(Z, expression(+ {} + {}))

Vts(X, expression(- {} + {}))
Vts(Y, expression(- {} - {}))
Vts(Z, expression(+ {} + {}))


Vpts(X, expression(- {} + {}))
Vpts(Y, expression(- {} - {}))
Vpts(Z, expression(+ {} + {}))


