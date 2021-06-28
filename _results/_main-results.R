

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

# Empty ggplot object, used in `plot_grid`
BLANK <- ggplot() + geom_blank() + theme_nothing()

# whether to re-do simulations (use rds files otherwise)
.REDO_SIMS <- FALSE
# whether to re-save plots
.RESAVE_PLOTS <- FALSE
# number of threads to use for simulations
.N_THREADS <- max(1, parallel::detectCores()-2)
# `mclapply` automatically uses this many threads
options(mc.cores = .N_THREADS)



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
    if (inherits(plot_obj, "function")) {
        cairo_pdf(fn, width = .width, height = .height)
        plot_obj()
        dev.off()
    } else {
        plot_fun <- ifelse(inherits(plot_obj, "egg"), print, plot)
        cairo_pdf(fn, width = .width, height = .height)
        plot_fun(plot_obj)
        dev.off()
    }
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

# 'Scaled community size'
comm_size_fun <- function(N, V1, V2, axes, d = NULL) {
    if (is.null(d)) .ds <- .d(.axes = axes[1]) else .ds = d
    sum(N * exp(- .ds[1] * V1^2 - .ds[2] * V2^2))
}


# Arguments for `egg::ggarrange`
label_args <- list(gp = grid::gpar(font = 1,
                                   fontsize = 16),
                   x = unit(0, "npc"),
                   y = unit(1, "npc"),
                   hjust = 0, vjust = 1)


pivot <- function(.df) {
    if ("pheno" %in% colnames(.df)) {
        .df <- rename(.df, V = geno, Vp = pheno) %>%
            pivot_wider(names_from = axis, values_from = c(V, Vp),
                        names_sep = "")
    } else {
        .df <- mutate(.df, axis = paste0("V", axis)) %>%
            pivot_wider(names_from = axis, values_from = geno)
    }
    return(.df)
}


#'
#' Make into `gtable` grob.
#'
make_gf <- function(x, ...) {
    x <- ggplotGrob(x)
    x <- gtable_frame(x, ...)
    return(x)
}



set.seed(1898348146)
etas <- map(1:6, ~ with(list(.q = .x), runif(.q * ((.q - 1) / 2), 0.1, 0.4)))
# With just one eta, it can be a simple number:
etas[[2]] <- 0.6








# =============================================================================*
# =============================================================================*

# Fig 2: Two-axis outcomes ----

# =============================================================================*
# =============================================================================*







one_eta_combo <- function(signs, .d, .n, ...) {


    # signs = 0; .d = -1
    # rm(signs, .d)

    stopifnot(is.numeric(.n) && .n > 0 && round(.n) == .n)
    stopifnot(is.numeric(signs))
    stopifnot(sum(is.na(signs)) == 0)
    stopifnot(all(signs %in% c(-1, 0, 1)))

    .q <- 1/2 * (1 + sqrt(1 + 8 * length(signs)))

    stopifnot(.q %in% 1:length(etas))

    .etas <- etas[[.q]]

    stopifnot(length(signs) == length(.etas))

    stopifnot(is.numeric(.d))
    stopifnot(all(!is.na(.d)))
    stopifnot(length(.d) %in% c(1L, .q))

    C <- matrix(0, .q, .q)
    C[lower.tri(C)] <- abs(.etas) * signs
    C <- C + t(C)
    diag(C) <- 1

    .V0_i <- 0.1 + 0.01 * (1:.q-1)

    .V0 <- unname(cbind(.V0_i, rev(.V0_i))[, rep(1:2, ceiling(.n / 2))[1:.n] ])

    args <- list(q = .q, eta = C, V0 = .V0,
                 n = .n,
                 d = .d, n_reps = 1,
                 spp_gap_t = 0L,
                 final_t = 50e3L,
                 save_every = 0L,
                 sigma_V0 = 0,
                 n_threads = .N_THREADS,
                 show_progress = FALSE)

    other_args <- list(...)
    if (length(other_args) > 0) {
        stopifnot(!is.null(names(other_args)))
        stopifnot(all(names(other_args) != ""))
        for (n in names(other_args)) args[[n]] <- other_args[[n]]
    }

    axis_to <- do.call(quant_gen, args)


    NV <- axis_to$nv %>%
        select(rep, spp, axis, geno) %>%
        pivot() %>%
        group_by(rep) %>%
        mutate(n_spp = length(unique(spp))) %>%
        ungroup()

    if (any(NV$n_spp != .n)) {
        message(sprintf(paste("\nsigns = %s, .d = %s, .n = %s had at least",
                              "one extinction"),
                        paste(signs, collapse = ", "),
                        paste(.d, collapse = ", "), .n))
    }

    for (i in 1:length(.etas)) {
        NV[,paste0("eta",i)] <- abs(.etas[i]) * signs[i]
    }
    NV[,"d"] <- .d

    output <- list(jacs = jacobians(axis_to),
                   ts = NV)


    return(output)
}



if (.REDO_SIMS) {
    # Takes ~9 sec with q=2 and 6 threads
    set.seed(145746935)
    eta_sims <- crossing(signs = c(-1, 0, 1), .d = c(-1, 0, 1), .n = 1:10) %>%
        pmap(one_eta_combo)
    saveRDS(eta_sims, rds("eta_sims"))
} else {
    eta_sims <- readRDS(rds("eta_sims"))
}

eta_sim_df <- map_dfr(eta_sims, ~.x$ts) %>%
    mutate(eta1 = factor(eta1, levels = sort(unique(eta1)),
                         labels = c("sub-additive", "additive",
                                    "super-additive")),
           d = factor(d, levels = c(-1, 0, 1),
                      labels = c("conflicting", "neutral", "ameliorative")),
           n_spp = factor(n_spp, levels = 1:max(n_spp)))

# eta_complex <- map_dfr(eta_sims, function(.z) {
#     zz <- map(.z$jacs, ~ eigen(.x, only.values = TRUE)$values)
#     zzz <- map_int(zz, function(.x) {
#         if (any(which(Im(.x) != 0))) {max(which(Im(.x) != 0))} else NA
#     })
#     zzz <- zzz[!is.na(zzz)]
#     if (length(zzz) == 0) {
#         zzz <- 0
#     } else zzz <- min(zzz)
#     distinct(.z$ts, eta1, d, n_spp) %>%
#         mutate(min_ind = zzz) %>%
#         filter(zzz > 0)
# })
# eta_stability <- map_dfr(eta_sims, function(.z) {
#     zz <- map(.z$jacs, ~ eigen(.x, only.values = TRUE)$values)
#     zzz <- map_dbl(zz, function(.x) max(Re(.x)))
#     distinct(.z$ts, eta1, d, n_spp) %>%
#         mutate(low = min(zzz), high = max(zzz))
# })





#'
#' Based on eigenvalues...
#'   * When the tradeoff is additive, the state is neutrally stable
#'   * Everything else is stable
#'   * Some eigenvalues are complex with d = -1, especially with n_spp = 2
#'
#'


# # If you want to show the relationship between radius and # species
# eta_sim_df %>%
#     mutate(eta1 = factor(eta1, levels = sort(unique(eta1)),
#                          labels = c("sub-additive", "additive",
#                                     "super-additive")),
#            d = factor(d, levels = c(-1, 0, 1),
#                       labels = c("conflicting", "neutral", "ameliorative")),
#            radius = sqrt(V1^2 + V2^2)) %>%
#     group_by(eta1, d, n_spp) %>%
#     filter(unq_spp_filter(V1, V2, .prec = 0.1)) %>%
#     ungroup() %>%
#     ggplot(aes(n_spp, radius)) +
#     geom_vline(xintercept = 1, linetype = 1, color = "gray70") +
#     geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
#     geom_point() +
#     facet_grid(eta1 ~ d)

outcomes_q2_p <- eta_sim_df %>%
    filter(!(eta1 == "additive" & n_spp != 1)) %>%
    filter(n_spp %in% 1:5) %>%
    filter(d == "neutral") %>%
    group_by(eta1, d, n_spp) %>%
    filter(unq_spp_filter(V1, V2, .prec = 0.001)) %>%
    ungroup() %>%
    ggplot(aes(V1, V2)) +
    geom_path(data = eta_sim_df %>%
                  filter(eta1 == "additive" & n_spp != 1) %>%
                  filter(n_spp %in% 1:5) %>%
                  filter(d == "neutral") %>%
                  group_by(eta1, d, n_spp) %>%
                  summarize(radius = median(sqrt(V1^2 + V2^2)),
                            .groups = "drop") %>%
                  mutate(V1 = map(radius,
                                  ~ .x * sin(seq(0, 0.5, length.out = 1001) *
                                                 pi)),
                         V2 = map(radius,
                                  ~ .x * cos(seq(0, 0.5, length.out = 1001) *
                                                 pi))) %>%
                  unnest(cols = c(V1, V2)),
              aes(color = n_spp), size = 0.75) +
    geom_point(aes(color = n_spp),
               alpha = 1, size = 2) +
    scale_x_continuous("Neutral axis 1", breaks = 0:2) +
    scale_y_continuous("Neutral axis 2", breaks = 0:2) +
    coord_equal(xlim = c(0, 2), ylim = c(0, 2)) +
    # facet_grid(d ~ eta1) +
    facet_grid( ~ eta1) +
    scale_color_viridis_d("Species in\ncommunity", end = 0.8) +
    theme(strip.text.x = element_text(size = 10, margin = margin(0,0,0,b=6)),
          strip.text.y = element_text(size = 10, margin = margin(0,0,0,l=6)),
          panel.border = element_rect(size = 0.5, fill = NA),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9, lineheight = 0.75),
          legend.key.size = unit(3, "pt"),
          legend.position = "right") +
    NULL



if (.RESAVE_PLOTS) save_plot(outcomes_q2_p, 5, 2, .prefix = "2-")









# =============================================================================*
# =============================================================================*

# Fig 3: History and axis strength ----

# =============================================================================*
# =============================================================================*




find_stable_comms <- function(.dat) {

    .eta <- .dat$.eta
    .d1 <- .dat$.d1
    .d2 <- .dat$.d2
    .state <- .dat$.state

    # .eta = 0.6
    # .d1 = -0.6
    # .d2 = 0.7
    # .state = "above"
    # rm(.eta, .d1, .d2, .state, .V0_0, sim0, L, out)

    .state <- match.arg(.state, c("above", "below", "both"))

    # This scenario takes a while:
    if (.eta == -0.6 && .d1 == -0.6 && .d2 > 2.38 && .d2 < 2.41 &&
        .state == "both") {
        .time <- 10e6L
    # Others only need this much time:
    } else .time <- 200e3L

    .V0_0 = cbind(c(0.3, 1), c(0.3, 1))
    if (.state == "below" || .state == "both") .V0_0 <- .V0_0[2:1,]
    if (.state == "both") .V0_0[,1] <- c(0.001, 2)

    sim0 <- quant_gen(q = 2, eta = .eta, V0 = .V0_0,
                      n = 2, d = c(.d1, .d2), n_reps = 1,
                      spp_gap_t = 0L,
                      final_t = .time,
                      save_every = 0L,
                      sigma_V0 = 0,
                      show_progress = FALSE)
    sim0$call[["eta"]] <- .eta
    sim0$call[["V0"]] <- .V0_0
    sim0$call[["d"]] <- c(.d1, .d2)

    if (all(!is.na(sim0$nv$spp))) {
        L <- eigen(jacobians(sim0)[[1]], only.values = TRUE)[["values"]] %>%
            Re() %>% max()
    } else L <- NA_real_

    out <- sim0$nv %>%
        pivot() %>%
        select(-rep) %>%
        mutate(O = comm_size_fun(N, V1, V2, d = c(.d1, .d2))) %>%
        mutate(eta = .eta, d1 = .d1, d2 = .d2, state = .state,
               lambda = L) %>%
        select(eta, d1, d2, state, everything())

    return(out)

}




if (.REDO_SIMS) {
    # Takes ~2 min
    find_comms_df <- bind_rows(
        crossing(.eta = c(-1,1)*0.6,
                 .d1 = seq(0.01, 2.7, 0.01) %>% round(2) %>% `*`(-1),
                 .d2 = 0.6,
                 .state = c("above", "below", "both")),
        crossing(.eta = c(-1,1)*0.6,
                 .d1 = -0.6,
                 .d2 = seq(0.01, 2.7, 0.01) %>% round(2),
                 .state = c("above", "below", "both"))) %>%
        distinct(.eta, .d1, .d2, .state) %>%
        filter(!(.eta < 0 & .state == "below")) %>%
        split(1:nrow(.)) %>%
        mclapply(find_stable_comms, mc.cores = .N_THREADS) %>%
        do.call(what = bind_rows)
    saveRDS(find_comms_df, rds("find_comms_df"))
} else {
    find_comms_df <- readRDS(rds("find_comms_df"))
}



#'
#' Returns a vector of which communities are duplicated
#' `X` is a data frame containing info for multiple communities,
#' in wide format (i.e., with 1 row per species, and
#' containing columns `V1` and `V2` instead of just `V`).
#' `.sep_col` is a string specifying the column that separates communities.
#' `.prec` is precision distinguishing separate communities.
#'
rm_dup_comm <- function(X, .sep_col = "state", .prec = 1e-4) {
    .V <- X %>%
        arrange(!!.sep_col, V1, V2) %>%
        select(V1, V2) %>%
        mutate(across(everything(), ~ round(., 20))) %>%
        split(X[[.sep_col]]) %>%
        lapply(t)
    stopifnot(length(unique(sapply(.V, ncol))) == 1)
    n <- length(.V)
    if (n == 1) return(X)
    dups <- rep(FALSE, n)
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            .z <- sum((.V[[i]] - .V[[j]])^2)
            if (.z < .prec) dups[j] <- TRUE
        }
    }
    unqs <- rep(!dups, each = ncol(.V[[1]]))
    return(X[unqs,])
}


# Unique stable communities:
# Takes ~6 sec
stable_comms_df <- find_comms_df %>%
    filter(lambda < 1) %>%
    select(-lambda) %>%
    split(interaction(.$eta, .$d1, .$d2, drop = TRUE)) %>%
    map_dfr(rm_dup_comm) %>%
    # group_by(eta, d1, d2, state) %>%
    # mutate(state = sum(V1 > 1e-6) %>% factor(levels = 0:2)) %>%
    # ungroup() %>%
    # rename(n_conf = state) %>%
    identity()






# ---------------------------------------------*
# Points (i.e., communities) to emphasize:
# ---------------------------------------------*

focal_stable_comms <- bind_rows(
    # First plot: effects of conflicting axis
    stable_comms_df %>%
        filter(d2 == 0.6, eta > 0, state == "below") %>%
        filter(d1 == round(min(d1) + 0.1, 2)),
    stable_comms_df %>%
        filter(d2 == 0.6, eta < 0, d1 > -1.7, state == "above") %>%
        filter(d1 == round(min(d1) + 0.1, 2)),
    # First plot: effects of ameliorative axis
    stable_comms_df %>%
        filter(d1 == -0.6, eta > 0, state == "both") %>%
        filter(d2 == 0.59),
    stable_comms_df %>%
        filter(d1 == -0.6, eta > 0, state == "above") %>%
        filter(d2 == round(max(d2) - 0.1, 2)),
    stable_comms_df %>%
        filter(d1 == -0.6, eta < 0, state == "above") %>%
        filter(d2 == round(max(d2) - 0.1, 2)),
    stable_comms_df %>%
        filter(d1 == -0.6, eta < 0, state == "both") %>%
        filter(d2 == round(min(d2) + 0.1, 2))) %>%
    mutate(comm = rep(1:(n()/2), each = 2) %>%
               as.roman() %>%
               tolower() %>%
               sprintf(fmt="(%s)"))




# Make fitness landscapes for these communities:

make_stable_comm_fit <- function(.dd) {

    # .dd <- stable_comms_df %>%
    #     filter(d1 == -0.6, eta < 0, state == "both") %>%
    #     filter(d2 == round(min(d2) + 0.1, 2))

    # rm(.eta, .strong, .barely, .N, .V, .q, .n, C, D, f, a0, r0, ..comm_fit,
    #    ..comm_exist, X)

    .eta <- .dd$eta[1]
    .d1 <- .dd$d1[1]
    .d2 <- .dd$d2[1]
    .N <- .dd$N
    .V <- rbind(.dd$V1, .dd$V2)
    .q <- 2
    .n <- ncol(.V)


    C <- diag(.q)
    C[1,2] <- C[2,1] <- .eta
    D <- diag(c(.d1, .d2))

    f <- formals(quant_gen)[["f"]]
    a0 <- formals(quant_gen)[["a0"]]
    r0 <- formals(quant_gen)[["r0"]]

    ..comm_fit <- function(.v1, .v2) {
        sauron:::F_it_cpp(V = cbind(c(.v1, .v2), .V), N = c(1, .N), i = 0,
                          f = f, a0 = a0, C = C, r0 = r0, D = D)
    }
    # Returns # species:
    ..comm_coexist <- function(.v1, .v2) {
        # .v1 = 1; .v2 = 0.5
        # rm(.v1, .v2)
        Z <- quant_gen(q = .q, eta = .eta,
                       V0 = cbind(c(.v1, .v2), .V), N0 = c(1, .N),
                       n = .n + 1, d = diag(D),
                       n_reps = 1,
                       spp_gap_t = 0L,
                       final_t = 50e3L,
                       save_every = 0L,
                       sigma_V0 = 0,
                       show_progress = FALSE) %>%
            .[["nv"]] %>%
            filter(axis == 1)
        nrow(Z)
    }

    X <- crossing(V1 = seq(0, 3, 0.1) %>% round(2), V2 = V1) %>%
        mutate(fit = map2_dbl(V1, V2, ..comm_fit)) %>%
        mutate(surv = fit >= 1, coexist_spp = .n,
               eta = .eta, d1 = .d1, d2 = .d2,
               state = .dd$state[1]) %>%
        select(eta, d1, d2, state, everything())

    X$coexist_spp[X$surv] <- map2_int(X$V1[X$surv], X$V2[X$surv],
                                      ..comm_coexist)

    lvls <- c("invader excluded",
              "2 residents excluded",
              "1 resident excluded",
              "3 species coexist")
    X <- X %>%
        mutate(outcome = case_when(!surv ~ lvls[1],
                                   coexist_spp == 1 ~ lvls[2],
                                   coexist_spp == 2 ~ lvls[3],
                                   coexist_spp == 3 ~ lvls[4]))
    stopifnot(!any(is.na(X$outcome)))
    X <- X %>%
        mutate(outcome = factor(outcome, levels = lvls))

    return(X)

}


if (.REDO_SIMS) {
    # Takes ~1 min
    stable_comm_fit <- focal_stable_comms %>%
        split(.$comm) %>%
        mclapply(make_stable_comm_fit) %>%
        do.call(what = bind_rows)
    saveRDS(stable_comm_fit, rds("stable_comm_fit"))
} else {
    stable_comm_fit <- readRDS(rds("stable_comm_fit"))
}



add_col_eta_cols <- function(.df) {
    .df %>%
        mutate(
            eta = factor(eta, levels = c(-1,1) * 0.6,
                         labels = c("Sub-additive", "Super-additive")),
            col = interaction(state, eta, drop = TRUE) %>%
                factor(levels = c(paste0(c("both", "above"), ".Sub-additive"),
                                  paste0(c("both", "above", "below"),
                                         ".Super-additive")) %>%
                           rev() %>%
                           .[c(4:1, 5)]))
}




#'
#' Plot for relationship between d and scaled community size.
#'
comms_d_p_fun <- function(.x_d, ...) {

    # .x_d = 1
    # rm(.x_d, .xlab, .other_d, .dd, .ddd, .p)

    .xlab <- paste(c("Conflicting", "Ameliorative")[.x_d], "axis strength")

    .other_d <- paste0("d", ifelse(.x_d == 1, 2, 1)) %>%
        as.name()
    .x_d <- paste0("d", .x_d) %>% as.name()

    .dd <- stable_comms_df %>%
        filter(d1 > -1.7) %>%
        filter(spp == 1, abs(!!.other_d) == 0.6) %>%
        select(eta, !!.x_d, state, O) %>%
        add_col_eta_cols()
    .ddd <- .dd %>%
        group_by(eta, state) %>%
        filter(!!.x_d == max(!!.x_d) | !!.x_d == min(!!.x_d)) %>%
        ungroup() %>%
        filter(!!.x_d != min(.dd[[paste(.x_d)]]),
               !!.x_d != max(.dd[[paste(.x_d)]]))
    .p <- .dd %>%
        ggplot(aes(abs(!!.x_d), O, color = col, linetype = eta)) +
        geom_line(size = 0.75) +
        geom_point(data = .ddd, shape = 4, size = 3) +
        scale_linetype_manual(values = c(2, 1)) +
        scale_color_viridis_d(begin = 0.2, end = 0.8, option = "A",
                              direction = -1, drop = FALSE) +
        scale_y_continuous("Scaled community size", labels = comma) +
        xlab(.xlab) +
        theme(...) +
        NULL
    return(.p)
}


#'
#' Heatmap for a focal community.
#'
focal_comm_hmp_fun <- function(.eta, .d1, .d2, .state, ...) {

    # .eta = -0.6
    # .d1 = -1.69
    # .d2 = 0.6
    # .state = "above"
    # rm(.eta, .d1, .d2, .state, .title, .p, .pal, .comm_df, .fit_df)

    .title <- focal_stable_comms %>%
        filter(eta == .eta, d1 == .d1, d2 == .d2, state == .state, spp == 1) %>%
        .[["comm"]]
    .pal <- magma(n = stable_comm_fit %>% .[["outcome"]] %>%
                      levels() %>% length() %>% `-`(1),
                  begin = 0.2, end = 0.8) %>%
        rev()

    .comm_df <- stable_comms_df %>%
        filter(eta == .eta, d1 == .d1, d2 == .d2, state == .state) %>%
        mutate(grp = group_spp(V1, V2, .prec = 0.001)) %>%
        group_by(eta, d1, d2, state, grp) %>%
        summarize(V1 = mean(V1), V2 = mean(V2), N = n(), O = O[1],
                  .groups = "drop") %>%
        select(-grp)
    .fit_df <- stable_comm_fit %>%
        filter(eta == .eta, d1 == .d1, d2 == .d2, state == .state) %>%
        select(V1, V2, fit, outcome)

    .p <- .comm_df %>%
        ggplot(aes(V1, V2)) +
        ggtitle(.title) +
        geom_raster(data = .fit_df,
                    aes(fill = outcome), interpolate = FALSE) +
        geom_abline(slope = 1, intercept = 0, linetype = 2,
                    color  = "gray70") +
        geom_point(size = 5, shape = 21, color = "black", fill = "white") +
        geom_text(aes(label = N), size = 10 / 2.83465,
                  fontface = "bold", color = "black") +
        coord_equal(xlim = c(-0.2, 3), ylim = c(-0.2, 3)) +
        scale_fill_manual(NULL, values = c("white", .pal), drop = FALSE,
                          aesthetics = c("color", "fill")) +
        xlab("Conflicting axis") +
        ylab("Ameliorative axis") +
        theme(plot.title = element_text(hjust = 0, size = 12,
                                        margin = margin(0,0,0,b=3)),
              legend.key = element_rect(colour = "black"),
              legend.text = element_text(size = 8),
              legend.key.size = unit(0.75, "lines"),
              plot.margin = margin(0,0,0,0)) +
        theme(...) +
        NULL

    return(.p)
}



comms_d_p_list <- map(1:2, comms_d_p_fun,
                      legend.position = "none",
                      axis.title.y = element_blank())

comms_d_p_list[[1]] <- comms_d_p_list[[1]] +
    geom_point(data = filter(focal_stable_comms, d2 == 0.6, spp == 1) %>%
                   add_col_eta_cols(),
         size = 2, shape = 1) +
    geom_text(data = filter(focal_stable_comms, d2 == 0.6, spp == 1) %>%
                  add_col_eta_cols(),
              aes(label = comm),
              size = 10 / 2.83465, nudge_x = -0.03, hjust = 1)
comms_d_p_list[[2]] <- comms_d_p_list[[2]] +
    geom_point(data = filter(focal_stable_comms, d1 == -0.6, spp == 1) %>%
                   add_col_eta_cols(),
                size = 2, shape = 1) +
    geom_text(data = filter(focal_stable_comms, d1 == -0.6, spp == 1) %>%
                  mutate(O = case_when(comm  != "(iv)" ~ O + 700,
                                       TRUE ~ O - 1600)) %>%
                  add_col_eta_cols(),
              aes(label = comm),
              size = 10 / 2.83465, vjust = 0)
comms_d_p <- ggarrange(plots = comms_d_p_list, draw = FALSE,
                       left = "Scaled community size")
# comms_d_p

focal_comm_p_list <- stable_comm_fit %>%
    distinct(eta, d1, d2, state) %>%
    rename_with(~ paste0(".", .x)) %>%
    pmap(focal_comm_hmp_fun, legend.position = "none",
         axis.title = element_blank())
focal_comm_p_list[c(2:3,5:6)]  <- focal_comm_p_list[c(2:3,5:6)] %>%
    map(~ .x + theme(axis.text.y = element_blank()))
focal_comm_p_list[1:3] <- focal_comm_p_list[1:3] %>%
    map(~ .x + theme(axis.text.x = element_blank()))



focal_comm_p <- arrangeGrob(textGrob("Ameliorative axis", rot = 90, vjust = 1),
                            gtable_rbind(focal_comm_p_list[1:3] %>%
                                             map(make_gf) %>%
                                             do.call(what = gtable_cbind),
                                         focal_comm_p_list[4:6] %>%
                                             map(make_gf) %>%
                                             do.call(what = gtable_cbind)),
                            make_gf(BLANK), textGrob("Conflicting axis"),
                            widths = c(0.05, 1), heights = c(1, 0.1),
                            nrow = 2, ncol = 2,
                            vp = viewport(height = unit(0.5, "npc"),
                                          y = 0, x = 0,
                                          just = c("left", "bottom")),
                            top = textGrob("(b)", x = 0, hjust = 0,
                                           gp = gpar(fontsize = 14)))




comm_p <- function() {
    grid.newpage()
    arrangeGrob(grobs = comms_d_p_list, left = "Scaled community size",
                vp = viewport(height = unit(0.5, "npc"),
                              y = 1, x = 0,
                              just = c("left", "top"))) %>%
        grid.draw()
    grid.draw(focal_comm_p)
    invisible(NULL)
}


if (.RESAVE_PLOTS) save_plot(comm_p, 4, 7, .prefix = "3-")








# =============================================================================*
# =============================================================================*

# ... old fig 3------

# =============================================================================*
# =============================================================================*


.d <- function(.strong, .barely, .axes = NULL) {
    if (!is.null(.axes)) {
        lvls <- c("C >> A", "C > A", "C < A", "C << A")
        stopifnot(.axes %in% lvls)
        ..df <- tibble(s = c(rep("c", 2), rep("a", 2)),
                       b = c(FALSE, TRUE, TRUE, FALSE))
        .strong <- ..df$s[lvls == .axes]
        .barely <- ..df$b[lvls == .axes]
    }
    .strong <- match.arg(.strong, c("conflicting", "ameliorative"))
    if (.strong == "conflicting") {
        z <- c(-2.7, 0.6)
        if (.barely) z[1] <- z[1] / 2
    } else {
        z <- c(-0.6, 2.7)
        if (.barely) z[2] <- z[2] / 2
    }
    return(z)
}


one_history_d_combo <- function(.dat) {

    .strong <- .dat$.strong
    .n <- .dat$.n
    .eta <- .dat$.eta
    .state <- .dat$.state
    .barely <- .dat$.barely


    # .strong = "conflicting"
    # .barely <- FALSE
    # .n = 2
    # .eta = 0.6
    # .state = "below"

    # rm(.dat, .dd, .strong, .n, .eta, .state, .V0_0, sim0, eig0, L, .N0, .V0, sims1)

    stopifnot(length(.strong) == 1 && .strong %in% c("conflicting", "ameliorative"))
    stopifnot(length(.n) == 1 && .n %in% 1:2)
    stopifnot(length(.eta) == 1 && .eta %in% (-1:1 * 0.6))
    stopifnot(length(.state) == 1 && .state %in% c("above", "below", "both"))
    stopifnot(length(.barely) == 1 && is.logical(.barely))

    .V0_0 = cbind(c(0.3, 1), c(0.3, 1))
    if (.n == 1) .V0_0 <- .V0_0[,1]
    if (.state == "below" || .state == "both") .V0_0 <- .V0_0[2:1,]
    if (.state == "both") .V0_0[,1] <- c(0.001, 2)

    .dd <- .d(.strong, .barely)

    sim0 <- quant_gen(q = 2, eta = .eta, V0 = .V0_0,
                      n = .n, d = .dd, n_reps = 1,
                      spp_gap_t = 0L,
                      final_t = 100e3L,
                      save_every = 0L,
                      sigma_V0 = 0,
                      show_progress = FALSE)


    if (any(is.na(sim0$nv$spp))) {
        msg <- sprintf(paste("Total extinction: strong = %s, barely = %s,",
                             "n = %i, eta = %.1f,",
                             "state = %s"), .strong, .barely, .n, .eta, .state)
        stop(msg)
    }

    sim0$call[["eta"]] <- .eta
    sim0$call[["V0"]] <- .V0_0
    sim0$call[["n"]] <- .n
    sim0$call[["d"]] <- .dd

    # You need to start with a stable community!
    J <- jacobians(sim0)[[1]]
    if (any(is.na(J))) {
        msg <- sprintf(paste("NA in Jacobian: strong = %s, barely = %s,",
                             "n = %i, eta = %.1f,",
                             "state = %s"), .strong, .barely, .n, .eta, .state)
        stop(msg)
    }
    eig0 <- eigen(J, only.values = TRUE)[["values"]]
    L <- max(Re(eig0))
    if ((L >= 1 && .eta != 0) || ((L - 1) > 1e-6 && .eta == 0)) {
        msg <- sprintf(paste("lambda = %.6f, strong = %s, barely = %s,",
                             "n = %i, eta = %.1f, state = %s"),
                       L, .strong, .barely, .n, .eta, .state)
        message(msg)
    }

    sim0 <- sim0 %>%
        .[["nv"]] %>%
        pivot()

    .N0 <- sim0$N
    .V0 <- sim0 %>%
        select(V1, V2) %>%
        t()

    sims1 <- tibble(.x = c(0.1, 0.1, 0.5),
                    .y = c(0.5, 1.5, 1)) %>%
        bind_rows(tibble(.x = .$.y, .y = .$.x)) %>%
        pmap_dfr(function(.x, .y) {
            quant_gen(q = 2, eta = .eta,
                      V0 = cbind(.V0, c(.x, .y)),
                      N0 = c(.N0, 1),
                      n = .n+1, d = .dd, n_reps = 1,
                      spp_gap_t = 0L,
                      final_t = 50e3L,
                      save_every = 10L,
                      sigma_V0 = 0,
                      show_progress = FALSE) %>%
                .[["nv"]] %>%
                pivot() %>%
                mutate(V1_0 = .x, V2_0 = .y) %>%
                select(-rep)
        }) %>%
        mutate(start = interaction(V1_0, V2_0, sep = "_"),
               eta = .eta, strong = .strong, barely = .barely,
               n_spp = .n, state = .state) %>%
        select(eta, strong, barely, n_spp, state, start, everything(), -V1_0, -V2_0)

    return(sims1)
}



if (.REDO_SIMS) {
    # Takes ~24 sec
    set.seed(199707675)
    hist_d_sims <- crossing(.strong = c("conflicting", "ameliorative"),
                            .barely = c(TRUE, FALSE),
                            .eta = -1:1 * 0.6,
                            .n = 1:2,
                            .state = c("above", "below", "both")) %>%
        # Remove combos that are boring or that can't produce equilibrium pops:
        filter(!(.state == "below" & (.eta < 0 | .n == 1)),
               !(.state == "both" & .n == 1),
               # !(.state == "both" & .eta < 0 & .strong == "ameliorative"),
               TRUE) %>%
        split(1:nrow(.)) %>%
        mclapply(one_history_d_combo, mc.cores = .N_THREADS) %>%
        do.call(what = rbind)
    saveRDS(hist_d_sims, rds("hist_d_sims"))
} else {
    hist_d_sims <- readRDS(rds("hist_d_sims"))
}





# ---------------*
# __unique stable communities ----
# ---------------*

check_stable <- function(.x) {

    if ("V" %in% colnames(.x) && is.list(.x$V)) {
        V <- .x$V[[1]]
        N <- .x$N[[1]]
    } else {
        V <- rbind(.x$V1, .x$V2)
        N <- .x$N
    }

    n <- ncol(V)
    q <- nrow(V)

    D <- matrix(0, q, q)
    diag(D) <- .d(.x$strong[1], .x$barely[1])
    C <- matrix(.x$eta[1], q, q)
    diag(C) <- 1

    J <- sauron:::jacobian_cpp(V = V, N = N,
                               formals(quant_gen)$f, formals(quant_gen)$a0,
                               formals(quant_gen)$r0, D, C,
                               eval(formals(quant_gen)$add_var), FALSE)
    J

    L <- max(Re(eigen(J, only.values = TRUE)[["values"]]))

    mutate(.x, lambda = L)

}


# Returns a vector of which communities are duplicated
dup_comm <- function(.V, .prec = 1e-4) {
    n <- length(.V)
    if (n == 1) return(FALSE)
    .V <- map(.V, ~ round(.x, 20))
    .V <- map(.V, ~ .x[,order(.x[1,], .x[2,])])
    dups <- rep(FALSE, n)
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            .z <- sum((.V[[i]] - .V[[j]])^2)
            if (.z < .prec) dups[j] <- TRUE
        }
    }
    return(dups)
}


# options(pillar.sigfig = 7)
# options(pillar.sigfig = NULL)

comms <- list()

# Two species communities
# --------------*
comms$two <- hist_d_sims %>%
    filter(time == 0 & eta != 0 & n_spp == 2 & spp %in% 1:2) %>%
    # filter(eta == 0.6, strong == "ameliorative", barely == TRUE, state == "both")
    # remove unstable communities:
    filter(!(barely == FALSE & strong == "ameliorative" & eta == -0.6 &
                     state == "above") &
               !(barely == FALSE & strong == "ameliorative" & eta == 0.6 &
                     state == "above") &
               !(barely == FALSE & strong == "conflicting" & eta == 0.6 &
                     state == "below") &
               !(barely == TRUE & strong == "ameliorative" & eta == 0.6 &
                     state == "above"),
               !(barely == TRUE & strong == "conflicting" & eta == 0.6 &
                     state == "below")) %>%
    group_by(eta, strong, barely, state, spp) %>%
    # Because these are repeated by invader starting points:
    summarize(V1 = V1[1], V2 = V2[1], N = N[1], .groups = "drop") %>%
    group_by(eta, strong, barely, state) %>%
    summarize(V = list(rbind(V1, V2)), N = list(N), .groups = "drop") %>%
    # split(1:nrow(.)) %>%
    # map_dfr(check_stable) %>%
    # arrange(desc(lambda))
    mutate(z1 = map_dbl(V, ~ sum(.x[1,])),
           z2 = map_dbl(V, ~ sum(.x[2,]))) %>%
    arrange(eta, strong, barely, (z1 - z2)^2) %>%
    group_by(eta, strong, barely) %>%
    mutate(dup = dup_comm(V)) %>%
    filter(!dup) %>%
    mutate(comm = 1:n()) %>%
    ungroup() %>%
    mutate(V1 = map(V, ~ .x[1,]),
           V2 = map(V, ~ .x[2,])) %>%
    select(-V, -z1, -z2, -state, -dup) %>%
    unnest(c(N, V1, V2)) %>%
    mutate(comm = factor(comm, levels = sort(unique(comm))),
           spp = factor(rep(1:2, n() / 2), levels = 1:2)) %>%
    select(eta, strong, barely, comm, spp, everything()) %>%
    identity()


#' These communities are not stable:
#'
#' lambda = 1.067329, strong = conflicting, barely = FALSE, eta = 0.6,
#' state = below
#'     - small perturbation and it goes to 1 species
#' lambda = 1.030829, strong = conflicting, barely = TRUE, eta = 0.6,
#' state = below
#'     - small perturbation and it goes to 1 species
#' lambda = 1.004679, strong = ameliorative, barely = FALSE, eta = 0.6,
#' state = above
#'     - small perturbation and it goes to one species investing in axis 2,
#'       the other not investing
#' lambda = 1.003003, strong = ameliorative, barely = TRUE, eta = 0.6,
#' state = above
#'     - small perturbation and it goes to one species investing 0 in V1
#'       and 1.26 in V2, the other investing ~0 in both
#' lambda = 1.001227, strong = ameliorative, barely = FALSE, eta = -0.6,
#' state = above
#'     - small perturbation and it goes to one species investing in both,
#'       the other not investing
#'


if (.REDO_SIMS) saveRDS(comms, rds("comms"))








# ---------------*
# __fitness landscapes ----
# ---------------*

# Make fitness landscapes for these communities:

comm_fit <- function(.dd) {

    # .dd <- comms$two %>%
    #     split(interaction(.$eta, .$strong, .$barely, .$comm, drop = TRUE)) %>%
    #     .[[1]]

    # rm(.eta, .strong, .barely, .N, .V, .q, .n, C, D, f, a0, r0, ..comm_fit,
    #    ..comm_exist, X)

    .eta <- .dd$eta[1]
    .strong <- .dd$strong[1]
    .barely <- .dd$barely[1]
    .N <- .dd$N
    .V <- rbind(.dd$V1, .dd$V2)
    .q <- nrow(.V)
    .n <- ncol(.V)


    C <- diag(.q)
    C[1,2] <- C[2,1] <- .eta
    D <- diag(.d(.strong, .barely))

    f <- formals(quant_gen)[["f"]]
    a0 <- formals(quant_gen)[["a0"]]
    r0 <- formals(quant_gen)[["r0"]]

    ..comm_fit <- function(.v1, .v2) {
        sauron:::F_it_cpp(V = cbind(c(.v1, .v2), .V), N = c(1, .N), i = 0,
                          f = f, a0 = a0, C = C, r0 = r0, D = D)
    }
    # Returns # species:
    ..comm_coexist <- function(.v1, .v2) {
        # .v1 = 1; .v2 = 0.5
        # rm(.v1, .v2)
        Z <- quant_gen(q = .q, eta = .eta,
                       V0 = cbind(c(.v1, .v2), .V), N0 = c(1, .N),
                       n = .n + 1, d = diag(D),
                       n_reps = 1,
                       spp_gap_t = 0L,
                       final_t = 50e3L,
                       save_every = 0L,
                       sigma_V0 = 0,
                       show_progress = FALSE) %>%
            .[["nv"]] %>%
            filter(axis == 1)
        nrow(Z)
    }

    X <- crossing(V1 = seq(0, 3, 0.1) %>% round(2), V2 = V1) %>%
        mutate(fit = map2_dbl(V1, V2, ..comm_fit)) %>%
        mutate(surv = fit >= 1, coexist_spp = .n,
               eta = .eta, strong = .strong, barely = .barely,
               comm = .dd$comm[1]) %>%
        select(eta, strong, comm, everything())

    X$coexist_spp[X$surv] <- map2_int(X$V1[X$surv], X$V2[X$surv],
                                      ..comm_coexist)

    return(X)

}





if (.REDO_SIMS) {
    hist_d_fit <- list()
    # Takes ~2 min
    hist_d_fit$two <- comms$two %>%
        split(interaction(.$eta, .$strong, .$barely, .$comm, drop = TRUE)) %>%
        mclapply(comm_fit, mc.cores = .N_THREADS) %>%
        do.call(what = rbind)
    saveRDS(hist_d_fit, rds("hist_d_fit"))
} else {
    hist_d_fit <- readRDS(rds("hist_d_fit"))
}


add_outcomes_col <- function(.df) {

    max_spp <- .df$coexist_spp %>% max()
    stopifnot(max_spp >= 3 && max_spp %% 1 == 0)

    lvls <- c("invader excluded",
              paste((max_spp-1):1, "residents excluded"),
              paste(max_spp, "species coexist"))
    .df_env <- as.environment(.df)
    parent.env(.df_env) <- environment()
    forms <- c(list(!surv ~ lvls[1]),
               map(2:length(lvls),
                   ~ sprintf("coexist_spp == %i ~ '%s'", .x-1, lvls[.x]) %>%
                       as.formula(env = .df_env)))
    environment(forms[[1]]) <- .df_env

    .df <- .df %>%
        mutate(outcome = do.call(case_when, forms))

    stopifnot(!any(is.na(.df$outcome)))

    .df <- .df %>%
        mutate(outcome = factor(outcome, levels = lvls))

    return(.df)

}

add_axes_col <- function(.df) {

    .df <- .df %>%
        mutate(axes = case_when(strong == "ameliorative" & !barely ~ "C << A",
                                strong == "ameliorative" & barely ~ "C < A",
                                strong == "conflicting" & !barely ~ "C >> A",
                                strong == "conflicting" & barely ~ "C > A"))

    stopifnot(!any(is.na(.df$axes)))

    .df <- .df %>%
        mutate(axes = factor(axes, levels = c("C >> A", "C > A",
                                              "C < A", "C << A")))
    return(.df)
}



hist_d_fit$two <- hist_d_fit$two %>%
    add_outcomes_col() %>%
    add_axes_col()




unq_comm_p_fun <- function(.e, .c, .a, ...) {

    # .c <- 1
    # .e <- -0.6
    # .a <- "C << A"
    # rm(.n, .c, .e, .xlab, .ylab, .title, .p, .pal)

    .n <- 2

    .nn <- c("one", "two", "three")[.n]

    .title <- case_when(.e > 0 ~ "Super-additive",
                        .e < 0 ~ "Sub-additive")
    .pal <- magma(n = hist_d_fit[[.nn]] %>% .[["outcome"]] %>%
                      levels() %>% length() %>% `-`(1),
                  begin = 0.2, end = 0.8) %>%
        rev()
    .pal <- magma(n = hist_d_fit[[.nn]] %>% .[["outcome"]] %>%
                        levels() %>% length() %>% `-`(1),
                    begin = 0.2, end = 0.8) %>%
        rev()

    .p <- comms %>%
        .[[.nn]] %>%
        filter(eta == .e, comm == .c) %>%
        add_axes_col() %>%
        filter(axes == .a) %>%
        select(spp, N, starts_with("V")) %>%
        mutate(grp = group_spp(V1, V2, .prec = 0.001)) %>%
        group_by(grp) %>%
        summarize(V1 = mean(V1), V2 = mean(V2), N = n(), .groups = "drop") %>%
        ggplot(aes(V1, V2)) +
        ggtitle(.title) +
        geom_raster(data = hist_d_fit[[.nn]] %>%
                        filter(eta == .e, comm == .c, axes == .a) %>%
                        select(V1, V2, fit, outcome),
                    aes(fill = outcome), interpolate = FALSE) +
        geom_abline(slope = 1, intercept = 0, linetype = 2,
                    color  = "gray70") +
        geom_point(size = 7, shape = 21, color = "black", fill = "white") +
        geom_text(aes(label = N), size = 12 / 2.83465,
                  fontface = "bold", color = "black") +
        coord_equal(xlim = c(-0.2, 3), ylim = c(-0.2, 3)) +
        scale_fill_manual(NULL, values = c("white", .pal), drop = FALSE,
                          aesthetics = c("color", "fill")) +
        xlab("Conflicting axis") +
        ylab("Ameliorative axis") +
        theme(plot.title = element_text(hjust = 0.5, size = 12,
                                        margin = margin(0,0,0,b=3)),
              legend.key = element_rect(colour = "black"),
              legend.text = element_text(size = 8),
              legend.key.size = unit(0.75, "lines"),
              plot.margin = margin(0,0,0,0)) +
        theme(...) +
        NULL

    return(tibble(eta = .e, comm = .c, axes = .a, p = list(.p)))
}


comms_p_df <- crossing(.e = c(-0.6, 0.6),
                         .c = 1:2,
                         .a = levels(hist_d_fit$two$axes)) %>%
    filter(!(.c == 2 & .e < 0)) %>%
    pmap_dfr(unq_comm_p_fun,
             legend.position = "none",
             plot.title = element_blank(),
             strip.text = element_blank(),
             axis.title = element_blank()) %>%
    mutate(p = ifelse(axes == "C << A", p,
                      map(p, ~ .x + theme(axis.text.x = element_blank()))),
           p = ifelse(eta <= 0, p,
                      map(p, ~ .x + theme(axis.text.y = element_blank()))),
           p = ifelse(axes != "C >> A", p,
                      map2(p, comm, ~ .x +
                              ggtitle(paste("comm.", .y)) +
                              theme(plot.title =
                                        element_text(size = 12, hjust = 0.5,
                                                     margin = margin(0,0,0,b=6)
                                        )))))


comms_p_comms_cols <- tibble(.c = c(1, 1, 2), .e = c(-1, 1, 1) * 0.6) %>%
    pmap(function(.c, .e) {
        # .c = 1; .e = -0.6
        z <- filter(comms_p_df, comm == .c, eta == .e)
        gtable_rbind(filter(z, axes == "C >> A")[["p"]][[1]] %>%
                         make_gf(height = unit(1, "null"),
                                 width = unit(1, "null")),
                     filter(z, axes == "C > A")[["p"]][[1]] %>%
                         make_gf(height = unit(1, "null"),
                                 width = unit(1, "null")),
                     filter(z, axes == "C < A")[["p"]][[1]] %>%
                         make_gf(height = unit(1, "null"),
                                 width = unit(1, "null")),
                     filter(z, axes == "C << A")[["p"]][[1]] %>%
                         make_gf(height = unit(1, "null"),
                                 width = unit(1, "null")))
    }) %>%
    {
        ..labs <- map(c(">>", ">", "<", "<<"),
                      ~ tibble(lab = paste0("conflicting ",.x,"\nameliorative"),
                               x = 0, y = 0) %>%
                          ggplot(aes(x, y, label = lab)) +
                          geom_text(hjust = 0.5, size = 12 / 2.83465) +
                          theme_nothing() +
                          coord_cartesian(xlim = c(NA, 0), expand = FALSE,
                                          clip = "off")) %>%
            map(make_gf, height = unit(1, "null"), width = unit(1, "null")) %>%
            do.call(what = gtable_rbind)
        ..blank1 <- replicate(4, make_gf(BLANK, width = unit(0.1, "null")))
        ..blank2 <- replicate(4, make_gf(BLANK, width = unit(0.5, "null")))
        gtable_cbind(..labs, do.call(gtable_rbind, ..blank1), .[[1]],
                     do.call(gtable_rbind, ..blank2), .[[2]], .[[3]])
    }

# grid.newpage(); grid.draw(comms_p_comms_cols)





# gtable_cbind
# gtable_rbind


comms_legend <- get_legend(comms_p_df$p[[1]] + theme(legend.position = "right"))

comm_lab <- function(i) textGrob(paste("comm.", i), gp = gpar(fontsize = 10))
comms_ylab <- textGrob("Ameliorative axis", just = c("center", "bottom"),
                       x = unit(0.7, "npc"), rot = 90)
axes_lab <- function(.comp) {
    textGrob(paste0("conflicting ", .comp, "\nameliorative"),
             gp = gpar(fontsize = 14))
}




main_block <- plot_grid(BLANK, BLANK, comm_lab(1), BLANK, comm_lab(1),
                            comm_lab(2),
                        comms_ylab, comms_p_list[[1]], BLANK, comms_p_list[[2]],
                            comms_p_list[[3]],
                        ncol = 6, rel_widths = c(0.5, 0.15, 1, 0.2, 1, 1),
                        nrow = 2, rel_heights = c(0.1, 1))
tradeoff_block <- plot_grid(BLANK,
                            textGrob("Sub-additive", gp = gpar(fontsize = 14)),
                            BLANK,
                            textGrob("Super-additive", gp = gpar(fontsize = 14)),
                            nrow = 1, rel_widths = c(0.15, 1, 0.2, 2))
xlab_block <- plot_grid(BLANK, textGrob("Conflicting axis"),
                        nrow = 1, rel_widths = c(0.15, 3.2))



# axes_plot <- crossing(strong = c("conflicting", "ameliorative"),
#          barely = c(TRUE, FALSE)) %>%
#     mutate(d = map2(strong, barely, .d)) %>%
#     unnest(d) %>%
#     mutate(var = rep(c("conflicting", "ameliorative"), n() / 2) %>%
#                factor(levels = sort(unique(.), decreasing = TRUE))) %>%
#     add_axes_col() %>%
#     ggplot(aes(var, abs(d))) +
#     geom_hline(yintercept = 0, linetype = 2, color  = "gray70") +
#     geom_segment(aes(xend = var, yend = 0)) +
#     geom_point(color = "black", fill = "dodgerblue", shape = 21, size = 3) +
#     scale_y_continuous("Axis strength", position = "right",
#                        limits = c(0, 1.1 *
#                                       max(map_dbl(c("c","a"),
#                                                   ~ max(abs(.d(.x, FALSE))))))) +
#     facet_grid(axes ~ .) +
#     theme(axis.title.x = element_blank(),
#           axis.text.x = element_text(size = 9, color = "black",
#                                      angle = 30, vjust = 0.6),
#           axis.ticks.x = element_blank(),
#           panel.spacing.y = unit(2, "lines"),
#           plot.margin = margin(0,0,0,0),
#           strip.text = element_blank())



comms_p <- plot_grid(plot_grid(tradeoff_block, main_block, xlab_block,
                                ncol = 1, rel_heights = c(0.04, 1.1, 0.09)),
                      plot_grid(comms_legend,
                                axes_plot,
                                BLANK,
                                # ncol = 1, rel_heights = c(0.1138211, 0.81300813,
                                #                           0.07317073)),
                                ncol = 1, rel_heights = c(0.1138211, 0.81300813,
                                                          0.07317073 - 0.041)),
                      nrow = 1, rel_widths = c(1, 0.3))




if (.RESAVE_PLOTS) save_plot(comms_p, 6.5, 8, .prefix = "3-")



#' The "scaled community size" parameter closely describes how conducive
#' the communities are to new invasion:
comms$two %>%
    add_axes_col() %>%
    mutate(d1 = map2_dbl(strong, barely, ~ .d(.x, .y)[1]),
           d2 = map2_dbl(strong, barely, ~ .d(.x, .y)[2])) %>%
    group_by(eta, d1, d2, comm) %>%
    summarize(comm_size = comm_size_fun(N, V1, V2, d = c(d1[1], d2[1])),
              .groups = "drop")










# =============================================================================*
# =============================================================================*

# Fig 4: Stabilizers ----

# =============================================================================*
# =============================================================================*


#' Also maybe show plot here about how it's not always most conducive to
#' coexistence when ameliorative axis is really strong.
#' (Because when it's really strong, all species investing in it isn't
#'  a stable community.)



stab_sims <- function(n_lower, d1, d2, .n = 10) {

    # n_lower = 4; d1 = -0.1; d2 = 0.5
    # rm(n_lower, d1, d2, .n, .V0_0, sim0, L, X)

    stopifnot(n_lower <= .n && n_lower >= 0)

    .V0_0 <- c(rep(list(c(1,0)), n_lower), rep(list(c(0,1)), .n - n_lower))
    .V0_0 <- do.call(cbind, .V0_0)

    sim0 <- quant_gen(q = 2, eta = 0.6, V0 = .V0_0,
                      n = .n, d = c(d1, d2), n_reps = 1,
                      spp_gap_t = 0L,
                      final_t = 100e3L,
                      save_every = 0L,
                      sigma_V0 = 0,
                      show_progress = FALSE)
    sim0$call$d <- c(d1, d2)
    sim0$call$n <- .n

    # stable?
    L <- sim0 %>%
        jacobians() %>%
        .[[1]] %>%
        eigen(only.values = TRUE) %>%
        .[["values"]] %>%
        Re() %>%
        max()
    if (L >= 1) return(rep(NA_real_, 3))

    X <- sim0 %>%
        .[["nv"]] %>%
        pivot()

    O <- map_dbl(1:nrow(X),
                 ~ X$N[.x] * exp(- d1 * X$V1[.x]^2 - d2 * X$V2[.x]^2)) %>%
        sum()


    V <- X %>%
        select(V1, V2) %>%
        map_dbl(max)

    return(c(V, Omega = O))
}

# Takes ~ 20 sec
stab_sim_df <- crossing(n_lower = 0:formals(stab_sims)$.n,
                        d1 = - c(0, 0.1, 0.5, 1),
                        d2 = abs(d1)) %>%
    mutate(Z = mcmapply(stab_sims, n_lower, d1, d2,
                         SIMPLIFY = FALSE, mc.cores = .N_THREADS)) %>%
    mutate(V1 = map_dbl(Z, ~ .x[1]),
           V2 = map_dbl(Z, ~ .x[2]),
           O = map_dbl(Z, ~ .x[3])) %>%
    select(-Z)

# Should all be zero:
stab_sim_df %>% filter(is.na(V1))
stab_sim_df %>% filter(is.na(V2))
stab_sim_df %>% filter(is.na(O))



#'
#' The only thing `d2` does  is make the `n_lower ~ V1` relationship steeper,
#' so the main-text plot below just uses `d2 == 0.5`.
#' Showing the effect of `d2` will go into the supplement.
#'
#' I'm also not using the case when `d1 == 0` because axis 1 is not actually
#' a conflicting axis in that case.
#' I simulated it just to understand the dynamics.
#'

stab_p_shared <-  list(
    geom_line(), geom_point(),
    scale_x_continuous("Species investing in conflicting axis",
                       breaks = 0:max(stab_sim_df$n_lower),
                       limits = c(0, max(stab_sim_df$n_lower))),
    scale_color_manual("Conflicting axis:",
                       values = c("black", "gray40", "gray70")),
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
          legend.position = "none"))


stab_p_list <- list()

stab_p_list[[1]] <- stab_sim_df %>%
    filter(d1 != 0, d2 == 0.5, n_lower > 0) %>%
    mutate(d1 = factor(d1, levels = sort(unique(d1)),
                        labels = c("strong", "moderate", "weak"))) %>%
    ggplot(aes(n_lower, V1, color = d1)) +
    scale_y_continuous("Conflicting\ninvestment\nper species",
                       limits = c(1.1, 2)) +
    stab_p_shared


stab_p_list[[2]] <- stab_sim_df %>%
    filter(d1 != 0, d2 == 0.5, n_lower < max(n_lower)) %>%
    mutate(d1 = factor(d1, levels = sort(unique(d1)),
                       labels = c("strong", "moderate", "weak"))) %>%
    ggplot(aes(n_lower, V2, color = d1)) +
    scale_y_continuous("Ameliorative\ninvestment\nper species",
                       limits = c(1.1, 2)) +
    stab_p_shared

stab_p_list[[3]] <- stab_sim_df %>%
    filter(d1 != 0, d2 == 0.5) %>%
    mutate(d1 = factor(d1, levels = sort(unique(d1)),
                       labels = c("strong", "moderate", "weak"))) %>%
    ggplot(aes(n_lower, O, color = d1)) +
    scale_y_continuous("Scaled\ncommunity\nsize", label = comma) +
    stab_p_shared


for (i in c(1,2)) {
    stab_p_list[[i]] <- stab_p_list[[i]] +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank())
}

stab_p_legend <- get_legend(stab_p_list[[1]] +
                                theme(legend.position = "top"))


stab_p <- plot_grid(plotlist = c(list(stab_p_legend), stab_p_list),
                    labels = c("", sprintf("(%s)", letters[1:length(stab_p_list)])),
                    label_fontface = "plain", label_size = 12,
                    align = "vh", axis = "lb",
                    rel_heights = c(0.2, 1, 1, 1), ncol = 1)

if (.RESAVE_PLOTS) save_plot(stab_p, 5, 5, .prefix = "4-")



#' Showing how ameliorative axis (i.e., `d2`) affects things.
#' Not as interesting, so it goes to the supplement.

# ... Fig S1: effects of ameliorative axis----

stab_p2_shared <- stab_p_shared
stab_p2_shared[[4]] <- scale_color_manual("Ameliorative axis:",
                                          values = c("black", "gray40",
                                                     "gray70"))

stab_p2_list <- list()

stab_p2_list[[1]] <- stab_sim_df %>%
    filter(d2 != 0, d1 == -0.5, n_lower > 0) %>%
    mutate(d2 = factor(d2, levels = sort(unique(d2), decreasing = TRUE),
                       labels = c("strong", "moderate", "weak"))) %>%
    ggplot(aes(n_lower, V1, color = d2)) +
    scale_y_continuous("Conflicting\ninvestment\nper species",
                       limits = c(0.96, 1.9)) +
    stab_p2_shared


stab_p2_list[[2]] <- stab_sim_df %>%
    filter(d2 != 0, d1 == -0.5, n_lower < max(n_lower)) %>%
    mutate(d2 = factor(d2, levels = sort(unique(d2), decreasing = TRUE),
                       labels = c("strong", "moderate", "weak"))) %>%
    ggplot(aes(n_lower, V2, color = d2)) +
    scale_y_continuous("Ameliorative\ninvestment\nper species",
                       limits = c(0.96, 1.9)) +
    stab_p2_shared

stab_p2_list[[3]] <- stab_sim_df %>%
    filter(d2 != 0, d1 == -0.5) %>%
    mutate(d2 = factor(d2, levels = sort(unique(d2), decreasing = TRUE),
                       labels = c("strong", "moderate", "weak"))) %>%
    ggplot(aes(n_lower, O, color = d2)) +
    scale_y_continuous("Scaled\ncommunity\nsize", label = comma) +
    stab_p2_shared

for (i in c(2,3)) {
    stab_p2_list[[i]] <- stab_p2_list[[i]] + theme(legend.position = "none")
}
for (i in c(1,2)) {
    stab_p2_list[[i]] <- stab_p2_list[[i]] +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank())
}


stab_p2 <- ggarrange(plots = stab_p2_list, ncol = 1, draw = FALSE)


if (.RESAVE_PLOTS) save_plot(stab_p2, 5, 5, .prefix = "S1-")




# =============================================================================*
# =============================================================================*

# Fig 5: Stochasticity ----

# =============================================================================*
# =============================================================================*



#'
#' .sigma_V should be wrapped into a list if you want it to vary by axis.
#'
one_stoch_sim <- function(.eta,
                          .sigma_V,
                          .d_mults = c(0.5, 1, 2),
                          .base_d = c(-0.25, 0.25),
                          ...) {

    # .eta = 1e-2; .sigma_V = 0.05

    stopifnot(length(.base_d) == 2)
    stopifnot(sign(.base_d) == c(-1,1))

    # Used for all calls to `quant_gen` below
    args <- list(eta = .eta,
                 q = 2,
                 n = 4,
                 V0 = cbind(c(0.2, 1), c(0.6, 0.7),
                            c(1, 0.2), c(0.7, 0.6)),
                 sigma_V0 = 0,
                 final_t = 200e3L,
                 spp_gap_t = 0L,
                 save_every = 100L,
                 sigma_V = unlist(.sigma_V),
                 n_reps = 12,
                 n_threads = .N_THREADS,
                 show_progress = FALSE)

    others <- list(...)
    for (n in names(others)) args[[n]] <- others[[n]]


    all_sims <- lapply(.d_mults, function(.m) {
        sim <- args %>%
            list_modify(d = .base_d * c(.m, 1)) %>%
            do.call(what = quant_gen) %>%
            .$nv %>%
            pivot() %>%
            mutate(d1 = .base_d[1] * .m, d2 = .base_d[2])
    })

    nv <- do.call(bind_rows, all_sims) %>%
        mutate(eta = .eta, sigma_V = .sigma_V) %>%
        select(eta, sigma_V, d1, d2, everything())

    return(nv)

}



#'
#' The function below produces a plot showing all simulation reps (from
#' simulation dataframes below).
#' I used it on all simulations to show that all reps conform to the
#' "representative" reps and that `d1` doesn't do anything useful here.
#' I've commented it out because it produces a very large, time-consuming plot.
#'
#'
#' stoch_big_plotter <- function(.sims) {
#'     .sims %>%
#'         ggplot(aes(V1, V2, color = spp)) +
#'         geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray70") +
#'         geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
#'         geom_vline(xintercept = 0, linetype = 1, color = "gray70") +
#'         geom_path() +
#'         geom_point(data = .sims %>% filter(time == 0),
#'                    shape = 1, size = 1) +
#'         geom_point(data = .sims %>% filter(time == max(time)),
#'                    shape = 19, color = "black", size = 3) +
#'         facet_grid(d1 + sigma_V ~ rep) +
#'         scale_color_viridis_d(end = 0.8, guide = FALSE) +
#'         coord_equal(xlim = c(0, 2.5), ylim = c(0, 2.5)) +
#'         NULL
#' }
#' # For example:
#' # stoch_big_plotter(stoch_comm_sims)
#'




#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#'
#' How does sigma_V affect relationship between conflicting axis strength
#' and scaled community size? (I predicted it makes it stronger.)
#'
#' ** Because this relationship depends on the specific scenario and community,
#' I'm not including this plot in the main text. **
#'
#'
# Takes ~20 sec
# set.seed(489421325)
# stoch_d1_sc_sims <- crossing(.eta = 1e-2,
#                        .sigma_V = c(0.05, 0.1, 0.2),
#                        .d_mults = round(seq(0, 0.5, 0.1), 1) / 0.25,
#                        save_every = 1000L) %>%
#     pmap_dfr(one_stoch_sim) %>%
#     select(-eta)
#
# stoch_d1_sc_sims %>%
#     filter(time > median(time)) %>%
#     group_by(sigma_V, d1, d2, rep, time) %>%
#     summarize(cs = comm_size_fun(N,
#                                  ifelse(is.na(Vp1), V1, Vp1),
#                                  ifelse(is.na(Vp2), V2, Vp2),
#                                  d = c(d1[1], d2[1])),
#               .groups = "drop") %>%
#     group_by(sigma_V, d1, rep) %>%
#     summarize(cs = median(cs)) %>%
#     group_by(sigma_V, d1) %>%
#     summarize(cs = median(cs), .groups = "drop") %>%
# mutate(sigma_V = factor(sigma_V, levels = sort(unique(sigma_V)),
#                         labels = sprintf("sigma[v] == %.2f",
#                                          sort(unique(sigma_V))))) %>%
#     ggplot(aes(abs(d1), cs)) +
#     geom_point() +
#     # geom_line() +
#     # geom_segment(aes(xend = d1, yend=0)) +
#     facet_grid(~ sigma_V, labeller = label_parsed) +
#     scale_y_continuous("Scaled community size", labels = comma,
#                        breaks = c(7,9,11)*1e3,
#                        limits = c(6500, 12500)) +
#     xlab("Conflicting axis strength") +
#     NULL







#'
#' Finds unique communities by `sigma_V`.
#' NOTE: This assumes that all reps are representative and that you're
#' not interested in `d1`.
#'
stoch_comms <- function(.sims) {
    .sims %>%
        filter(time == max(time), rep == 1, d1 == median(d1)) %>%
        group_by(sigma_V) %>%
        mutate(grp = group_spp(V1, V2, .prec = 1e-1)) %>%
        group_by(sigma_V, grp) %>%
        summarize(V1 = mean(V1), V2 = mean(V2), n_spp = n(), .groups = "drop")
}


#'
#' Plots how species evolve through axis space, for a given set of simulations
#' and `sigma_V`.
#'
stoch_comm_V_plotter <- function(.s, .sims) {
    .sims_comms <- stoch_comms(.sims) %>%
        filter(sigma_V == .s)
    .sims %>%
        filter(rep == 1, d1 == median(d1), sigma_V == .s) %>%
        ggplot(aes(V1, V2, color = spp)) +
        geom_abline(slope = 1, intercept = 0,
                    linetype = 2, color = "gray70") +
        geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
        geom_vline(xintercept = 0, linetype = 1, color = "gray70") +
        geom_path() +
        geom_point(data = .sims %>%
                       filter(time == 0, rep == 1, d1 == median(d1),
                              sigma_V == .s),
                   shape = 1, size = 1) +
        geom_point(data = .sims_comms, shape = 19, color = "black", size = 3) +
        geom_text(data = .sims_comms, aes(label = n_spp),
                  color = "white", size = 8 / 2.83465) +
        scale_color_viridis_d(end = 0.9, begin = 0.4, guide = FALSE) +
        scale_y_continuous("Ameliorative axis") +
        scale_x_continuous("Conflicting axis") +
        coord_equal(xlim = c(0, 1.5), ylim = c(0, 1.5)) +
        theme(plot.margin = margin(0,0,0,0)) +
        NULL
}
#'
#' Plots how Omega (scaled community size) changes through time,
#' for a given set of simulations and `sigma_V`.
#' NOTE: it filters out time < 1000 as we're most interested in
#' dynamics at "equilibrium".
#'
stoch_comm_O_plotter <- function(.s, .sims, .y_scale = NULL, .rollmean = FALSE) {
    .dd <- .sims %>%
        filter(time > 1000, rep == 1, d1 == median(d1), sigma_V == .s) %>%
        group_by(time) %>%
        summarize(cs = comm_size_fun(N,
                                     ifelse(is.na(Vp1), V1, Vp1),
                                     ifelse(is.na(Vp2), V2, Vp2),
                                     d = c(d1[1], d2[1])),
                  .groups = "drop")
    if (is.null(.y_scale)) {
        .y_scale <- scale_y_continuous("Scaled\ncommunity size",
                                       labels = comma,
                                       breaks = c(7,9,11)*1e3,
                                       limits = c(6500, 12500))
    }
    .p <- .dd %>%
        ggplot(aes(time, cs)) +
        geom_line(color = "black") +
        scale_color_viridis_d(end = 0.9, begin = 0.4, guide = FALSE) +
        .y_scale +
        scale_x_continuous("Time", breaks = c(0, 75e3, 150e3), labels = comma)
    if (.rollmean) {
        .p <- .p +
            geom_line(data = .dd %>%
                          mutate(cs = zoo::rollmean(cs, 50, fill = NA)) %>%
                          filter(!is.na(cs)), color = "dodgerblue")
    }
    return(.p)
}

#'
#' Plots PDF of normal distributions with varying `sigma_V`.
#' Set `.flip` to `TRUE` for plots on right side.
#'
sigma_distr_plotter <- function(.s,
                                .flip,
                                .max_y = 12,
                                .x_range = c(-0.6, 0.6),
                                .xlim_mult = 2) {
    .coord <- if (.flip) coord_flip else coord_cartesian
    tibble(x = seq(min(.x_range), max(.x_range), 0.001) %>% round(3)) %>%
        mutate(dens = map_dbl(x, ~ dnorm(.x, sd = .s))) %>%
        ggplot(aes(x, dens)) +
        geom_area(color = NA, fill = "gray60") +
        theme_nothing() +
        .coord(ylim = c(0, .max_y),
               xlim = .x_range * .xlim_mult,
               expand = FALSE)
}
#'
#' Combines plots of axis evolution, Omega, and sigmas into one plot.
#'
stoch_comm_p_combiner <- function(.V_plots,
                                  .O_plots,
                                  .sigmas) {

    # .V_plots = stoch_super_V_p_list
    # .O_plots = stoch_super_O_p_list
    # .sigmas =  c(0.05, 0.1, 0.2)
    # rm(.V_plots, .O_plots, .sigmas, i, combined_V, combined_O, combined)

    stopifnot(length(.V_plots) == length(.O_plots))

    stopifnot(inherits(.sigmas, "matrix") || inherits(.sigmas, "numeric"))

    if (inherits(.sigmas, "matrix")) {
        if (ncol(.sigmas)!=2) stop("\nIf a matrix, `.sigmas` must have 2 cols")
        if (nrow(.sigmas) != length(.V_plots)) {
            stop("\nIf `.sigmas` is a matrix, nrow(.sigmas) must equal ",
                 "`length(.V_plots)`")
            }
    } else {
        if (length(.sigmas) != length(.V_plots)) {
            stop("\nIf `.sigmas` is numeric, length(.sigmas) must equal ",
                 "`length(.V_plots)`")
        }
        .sigmas <- unname(cbind(.sigmas, .sigmas))
    }

    for (i in 1:length(.V_plots)) {
        if (i > 1) {
            .V_plots[[i]] <- .V_plots[[i]] +
                theme(axis.title.y = element_blank(),
                      axis.text.y = element_blank())
            .O_plots[[i]] <- .O_plots[[i]] +
                theme(axis.title.y = element_blank(),
                      axis.text.y = element_blank())
        }
        if (i != median(c(1, length(.V_plots)))) {
            .V_plots[[i]] <- .V_plots[[i]] +
                theme(axis.title.x = element_blank())
            .O_plots[[i]] <- .O_plots[[i]] +
                theme(axis.title.x = element_blank())
        }
    }

    combined_V <- lapply(1:nrow(.sigmas),
                         function(i){
                             main_plot <- .V_plots[[i]] %>%
                                 make_gf(width = unit(1, "null"))
                             top_sigma <- .sigmas[i,1] %>%
                                 sigma_distr_plotter(.flip = FALSE) %>%
                                 make_gf(width = unit(1, "null"))
                             right_sigma <- .sigmas[i,2] %>%
                                 sigma_distr_plotter(.flip = TRUE) %>%
                                 make_gf(width = unit(0.2, "null"))
                             blank <- BLANK %>%
                                 make_gf(width = unit(0.2, "null"))
                             top_row <- gtable_frame(
                                 gtable_cbind(blank, blank),
                                 height = unit(0.05, "null"))
                             mid_row <- gtable_frame(
                                 gtable_cbind(top_sigma, blank),
                                 height = unit(0.1, "null"))
                             bot_row <- gtable_frame(
                                 gtable_cbind(main_plot, right_sigma),
                                 height = unit(1, "null"))
                             combined <- gtable_rbind(top_row, mid_row, bot_row)
                             # grid.newpage(); grid.draw(combined)
                             return(combined)
                         }) %>%
        do.call(what = gtable_cbind)
    # grid.newpage(); grid.draw(combined_V)

    combined_O <- .O_plots %>%
        lapply(make_gf, width = unit(1, "null"), height = unit(0.2, "null")) %>%
        # lapply(function(.x) {
        #     gtable_cbind(.x,
        #                  make_gf(BLANK, width = unit(0.2, "null")))
        # }) %>%
        do.call(what = gtable_cbind) %>%
        identity()
    # grid.newpage(); grid.draw(combined_O)

    combined <- gtable_rbind(combined_V, combined_O)
    # grid.newpage(); grid.draw(combined)

    return(combined)

}



#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#'
#' For showing how sigma_V affects tradeoffs and resulting scaled
#' community size when tradeoffs are super-additive.
#'
# Takes ~20 sec
set.seed(1630767640)
stoch_super_sims <- crossing(.eta = 1e-2,
                            .sigma_V = c(0.05, 0.1, 0.2)) %>%
    pmap_dfr(one_stoch_sim) %>%
    select(-eta)


stoch_super_V_p_list <- map(sort(unique(stoch_super_sims$sigma_V)),
                            stoch_comm_V_plotter, .sims = stoch_super_sims)
# stoch_super_V_p_list[[3]]

stoch_super_O_p_list <- map(sort(unique(stoch_super_sims$sigma_V)),
                            stoch_comm_O_plotter, .sims = stoch_super_sims)
# stoch_super_O_p_list[[2]]

stoch_super_comm_p <- stoch_comm_p_combiner(stoch_super_V_p_list,
                                            stoch_super_O_p_list,
                                            c(0.05, 0.1, 0.2))
# grid.newpage(); grid.draw(stoch_super_comm_p)








#' ----------------------------------------------------------------------------
#' ----------------------------------------------------------------------------
#'
#' For showing how differences between sigma_V1 and sigma_V2 affect
#' tradeoffs and resulting scaled community size when tradeoffs are
#' sub-additive
#'
# Takes ~10 sec
set.seed(489421325)
stoch_sub_sims <- crossing(.eta = -1e-2,
                           .sigma_V = list(list(c(0.05, 0.05)),
                                           list(c(0.05, 0.20)),
                                           list(c(0.20, 0.05))),
                           .d_mults = 1) %>%
    pmap_dfr(one_stoch_sim) %>%
    select(-eta) %>%
    mutate(sigma_v1 = map_dbl(sigma_V, ~ .x[1]),
           sigma_v2 = map_dbl(sigma_V, ~ .x[2]),
           sigma_V = paste(sigma_v1, sigma_v2, sep = "__") %>%
               factor(levels = c("0.05__0.05", "0.05__0.2", "0.2__0.05")))


stoch_sub_V_p_list <- map(sort(unique(stoch_sub_sims$sigma_V)),
                          stoch_comm_V_plotter, .sims = stoch_sub_sims)
stoch_sub_O_p_list <- map(sort(unique(stoch_sub_sims$sigma_V)),
                          stoch_comm_O_plotter, .sims = stoch_sub_sims)

stoch_sub_comm_p <- stoch_comm_p_combiner(stoch_sub_V_p_list,
                                          stoch_sub_O_p_list,
                                          rbind(c(0.05, 0.05),
                                                c(0.05, 0.20),
                                                c(0.20, 0.05)))

# grid.newpage(); grid.draw(stoch_sub_comm_p)




stoch_p <- function() {
    p <- gtable_rbind(stoch_super_comm_p, stoch_sub_comm_p)
    grid.newpage()
    grid.draw(p)
    grid.text("(a)", x = 0, y = 1, just = c("left", "top"),
              gp = gpar(fontsize = 16))
    grid.text("Super-additive + equal variances",
              x = 0.1, y = 1, just = c("left", "top"),
              gp = gpar(fontsize = 14))
    grid.text("(b)", x = 0, y = 0.5, just = c("left", "top"),
              gp = gpar(fontsize = 16))
    grid.text("Sub-additive + unequal variances",
              x = 0.1, y = 0.5, just = c("left", "top"),
              gp = gpar(fontsize = 14))
    invisible(NULL)
}

if (.RESAVE_PLOTS) save_plot(stoch_p, 6, 6, .prefix = "5-")







# ----------------------------------------------------------------------------*
# ----------------------------------------------------------------------------*
# ... Fig S2: alt. scenario where sigma_V increases Omega ----


# Takes ~7 sec
set.seed(1412303795)
alt_stoch_super_sims <- crossing(.eta = 1e-2,
                                .sigma_V = c(0.05, 0.1, 0.2),
                                .d_mults = 1,
                                V0 = list(cbind(c(0.2, 1), c(0.6, 0.7),
                                                c(1, 0.2), c(0.5, 0.6)))) %>%
    pmap_dfr(one_stoch_sim) %>%
    select(-eta)


alt_stoch_super_V_p_list <- map(sort(unique(alt_stoch_super_sims$sigma_V)),
                                stoch_comm_V_plotter,
                                .sims = alt_stoch_super_sims)
alt_stoch_super_O_p_list <- map(sort(unique(alt_stoch_super_sims$sigma_V)),
                                stoch_comm_O_plotter,
                                .sims = alt_stoch_super_sims)

alt_stoch_super_comm_p <- stoch_comm_p_combiner(alt_stoch_super_V_p_list,
                                                alt_stoch_super_O_p_list,
                                                c(0.05, 0.1, 0.2))
# grid.newpage(); grid.draw(alt_stoch_super_comm_p)


alt_stoch_p <- function() {
    grid.newpage()
    grid.draw(alt_stoch_super_comm_p)
    grid.text("Super-additive + equal variances",
              x = 0.1, y = 1, just = c("left", "top"),
              gp = gpar(fontsize = 14))
}


if (.RESAVE_PLOTS) save_plot(alt_stoch_p, 6, 3, .prefix = "S2-")






#     # ===================================================================*
#     # ===================================================================*
#     # ===================================================================*
#     # ===================================================================*
#
#
#
#
#     zzz$nv %>%
#         # filter(rep == 1) %>%
#         # filter(spp == 1) %>%
#         pivot() %>%
#         ggplot(aes(time, V1, color = spp)) +
#         # geom_line() +
#         geom_line(aes(group = interaction(rep, spp)), alpha = 0.25) +
#         geom_point(data = zzz$nv %>% filter(time == 0, rep == 1) %>%
#                        # filter(spp == 1) %>%
#                        select(-rep) %>%
#                        pivot() %>%
#                        mutate(Vp1 = V1, Vp2 = V2),
#                    aes(fill = spp), size = 2, shape = 21, color = "black") +
#         facet_wrap(~ rep, ncol = 3) +
#         coord_cartesian(ylim = c(1.3, NA)) +
#         NULL
#
#
#
#
#     filter_comms <- function(.V1, .V2, .outcome, .comm) {
#         l1 <- .outcome != "invader excluded" &
#             (.V1 - mean(range(.V1)))^2 == min((.V1 - mean(range(.V1)))^2) &
#             (.V2 - mean(range(.V2)))^2 == min((.V2 - mean(range(.V2)))^2)
#
#         cmty <- communities[[.comm]]
#         .V1 - .V2
#         # l2 <-
#         .outcome == "invader excluded" &
#
#     }
#
#
#     hist_d_fit$two %>%
#         filter(strong == .strong, barely == .barely, eta == .eta) %>%
#         select(comm, V1, V2, outcome) %>%
#         group_by(comm, outcome) %>%
#         filter(filter_comms(V1, V2, outcome, comm)) %>%
#         summarize(V1 = mean(V1), V2 = mean(V2), .groups = "drop")
#
#
#         # arrange(comm) %>%
#         # split(.$comm)
#
#
#
#
#     if (.sigma_V == 0 && .sigma_N == 0) {
#         # These simulations are totally deterministic, so no need to run > once
#         .nreps <- 1
#         .N_THREADS <- 1
#     } else .nreps <- 96
#
#     .seed <- sample.int(2^31 - 1, 1)
#     set.seed(.seed)
#
#
#     one_V12_combo <- function(.v1, .v2) {
#         # .v1 = 1; .v2 = 1
#         # rm(.v1, .v2)
#         ..seed <- sample.int(2^31-1, 1)
#         set.seed(..seed)
#         X <- quant_gen(eta = .eta, d = .ds, q = 2, n = 2,
#                        V0 = cbind(t(stable_points(.eta)[1,]),
#                                   c(.v1, .v2)),
#                        sigma_V0 = 0,
#                        N0 = c(1, 1),
#                        spp_gap_t = 1e3L, final_t = 20e3L, save_every = 0L,
#                        sigma_V = .sigma_V,
#                        sigma_N = .sigma_N,
#                        add_var = rep(0.05, 2),
#                        n_reps = .nreps, n_threads = .N_THREADS,
#                        show_progress = FALSE)
#         # X$nv %>%
#         #     filter(trait == 1) %>%
#         #     group_by(rep) %>%
#         #     summarize(inv = any(spp == 2, na.rm = TRUE),
#         #               res = any(spp == 1, na.rm = TRUE)) %>%
#         #     select(-rep) %>%
#         #     summarise(across(.fns = mean)) %>%
#         #     mutate(V1 = .v1, V2 = .v2, seed = ..seed)
#         X$seed = ..seed
#         return(X)
#     }
#
#     Z <- crossing(V1 = seq(0, 4, 0.2),
#                   V2 = seq(0, 4, 0.2)) %>%
#         # bc `seq` makes weird numbers that are ~1e-15 from what they should be:
#         mutate(across(.fns = round, digits = 1)) %>%
#         mutate(sigma_V = .sigma_V, sigma_N = .sigma_N,
#                eta = .eta, d1 = .d1)
#     Z$qg = pmap(list(.v1 = Z$V1, .v2 = Z$V2), one_V12_combo)
#
#     return(Z)
# }
#
#
# crossing(.strong = c("c","a"),
#          .barely = c(TRUE, FALSE),
#          .eta = c(-1,1) * 0.6,
#          .sigma_N = c(0, 0.1, 0.2),
#          .sigma_V = c(0, 0.1, 0.2))





# =============================================================================*
# =============================================================================*

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# OBSOLETE Fig 3: "Global" coexistence ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----

# =============================================================================*
# =============================================================================*


#' # --------------*
#' # __3A - d --> # spp ----
#' # --------------*
#'
#' one_sim_combo <- function(.d, .eta, .add_var, .sigma_N, .sigma_V, .vary_d2) {
#'
#'     # .d = 0.5; .eta = 0.5; .add_var = 0.05; .sigma_N = 0.5; .sigma_V = 0; .vary_d2 = TRUE
#'
#'     .n <- 100
#'
#'     if (.vary_d2) {
#'         .ds <- rep(.d, 2)
#'     } else .ds <- c(.d, 0.1)
#'
#'     .seed <- sample.int(2^31 - 1, 1)
#'     set.seed(.seed)
#'
#'     Z <- quant_gen(q = 2, eta = .eta, d = .ds, n = .n,
#'                    add_var = rep(.add_var, .n),
#'                    # spp_gap_t = 500L,
#'                    spp_gap_t = 0L,
#'                    final_t = 50e3L, n_reps = 12,
#'                    sigma_N = .sigma_N, sigma_V = .sigma_V,
#'                    save_every = 0L,
#'                    n_threads = .N_THREADS, show_progress = FALSE)
#'
#'     if (.sigma_N == 0 && .sigma_V == 0) {
#'         Z$call[["eta"]] <- eval(.eta)
#'         Z$call[["d"]] <- eval(.ds)
#'         Z$call[["n"]] <- eval(.n)
#'         Z$call[["add_var"]] <- eval(rep(.add_var, .n))
#'         jacs <- jacobians(Z)
#'     } else jacs <- NULL
#'
#'
#'     if ("pheno" %in% colnames(Z$nv)) {
#'         .NV <- Z$nv %>%
#'             mutate(axis = paste0("V", axis)) %>%
#'             nest(value = c(geno, pheno)) %>%
#'             spread(axis, value) %>%
#'             unnest(c(V1, V2), names_sep = "_") %>%
#'             rename(V1 = V1_geno,
#'                    V2 = V2_geno,
#'                    Vp1 = V1_pheno,
#'                    Vp2 = V2_pheno)
#'     } else {
#'         .NV <- Z$nv %>%
#'             pivot()
#'         .NV$Vp1 <- NA_real_
#'         .NV$Vp2 <- NA_real_
#'     }
#'
#'     return(list(NV = .NV %>%
#'                     mutate(d = .d, eta = .eta, add_var = .add_var,
#'                            sigma_N = .sigma_N, sigma_V = .sigma_V,
#'                            vary_d2 = .vary_d2),
#'                 J = jacs,
#'                 seed = .seed))
#' }
#'
#' # Takes ~4 min w/ 6 threads and `.d` of length 7
#' coexist_d_spp_sims <- crossing(.d = -5:0,#seq(-0.15, 0, 0.025),
#'                             .eta = -1:1 * etas[[2]],
#'                             .add_var = 0.05,
#'                             .sigma_N = 0,
#'                             .sigma_V = 0,
#'                             .vary_d2 = TRUE) %>%
#'     pmap(one_sim_combo)
#'
#'
#' coexist_d_spp_df <- map_dfr(coexist_d_spp_sims, ~ .x[["NV"]]) %>%
#'     group_by(d, eta, rep) %>%
#'     summarize(n_spp = spp[N > 0] %>% unique() %>% length()) %>%
#'     ungroup() %>%
#'     mutate(eta = factor(eta, levels = -1:1 * etas[[2]],
#'                         labels = paste0(c("sub-", "", "super-"), "additive")))
#'
#'
#' # coexist_d_spp_p <-
#' coexist_d_spp_df %>%
#'     mutate(n_spp = n_spp / 100) %>%
#'     ggplot(aes(d, n_spp)) +
#'     geom_vline(xintercept = 0, color = "gray80", linetype = 1) +
#'     geom_hline(yintercept = 0, color = "gray80", linetype = 1) +
#'     geom_jitter(width = 0.002, height = 0, shape = 1, size = 1) +
#'     facet_wrap(~ eta, ncol = 1) +
#'     scale_y_continuous("Proportion of species that coexist",
#'                        breaks = c(0, 0.4, 0.8), limits = c(0, 1)) +
#'     # scale_x_continuous("Strength of\nconflicting axis",
#'     #                    breaks = c(-0.1, 0)) +
#'     theme(axis.title = element_text(size = 10),
#'           axis.text = element_text(size = 9),
#'           plot.margin = margin(0,l=8,r=8,t=8))
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' # --------------*
#' # __3B - sigma_V/N --> # spp ----
#' # --------------*
#'
#' # Takes ~3.5 min w/ 6 threads
#' coexist_stoch_spp_sims <- crossing(.d = 0,
#'                                   .eta = -1:1 * etas[[2]],
#'                                   .add_var = 0.05,
#'                                   .sigma_N = c(0, 0.05),
#'                                   .sigma_V = c(0, 0.05, 0.1),
#'                                   .vary_d2 = FALSE) %>%
#'     pmap(one_sim_combo)
#'
#'
#'
#' coexist_stoch_spp_df <- coexist_stoch_spp_sims %>%
#'     map_dfr(~ .x[["NV"]]) %>%
#'     group_by(eta, sigma_N, sigma_V, rep) %>%
#'     summarize(n_spp = spp[N > 0] %>% unique() %>% length()) %>%
#'     ungroup() %>%
#'     mutate(eta = factor(eta, levels = -1:1 * etas[[2]],
#'                         labels = c("sub-additive", "additive",
#'                                    "super-additive"))) %>%
#'     mutate(p_spp = n_spp / 100,
#'            sigma_N = factor(sigma_N, levels = sort(unique(sigma_N))))
#'
#'
#' # coexist_stoch_spp_p <-
#' coexist_stoch_spp_df %>%
#'     ggplot(aes(sigma_V, p_spp)) +
#'     geom_vline(xintercept = 0, color = "gray80", linetype = 1) +
#'     geom_hline(yintercept = 0, color = "gray80", linetype = 1) +
#'     geom_jitter(aes(color = sigma_N),
#'                 width = 0.002, height = 0,
#'                 size = 1, shape = 1) +
#'     geom_line(data = coexist_stoch_spp_df %>%
#'                   group_by(eta, sigma_N, sigma_V) %>%
#'                   summarize(p_spp = mean(p_spp), .groups = "drop"),
#'               aes(color = sigma_N)) +
#'     geom_text(data = coexist_stoch_spp_df %>%
#'                   distinct(eta, sigma_N) %>%
#'                   filter(eta == "additive") %>%
#'                   mutate(sigma_V = c(0.1, 0.1), p_spp = c(0.75, 0.01)) %>%
#'                   mutate(lab = paste("sigma[N] ==", sigma_N)),
#'               aes(color = sigma_N, label = lab),
#'               size = 8 / 2.83465, parse = TRUE, hjust = 1, vjust = 0) +
#'     facet_wrap(~ eta, ncol = 1) +
#'     scale_y_continuous("Proportion of species that coexist",
#'                        breaks = c(0, 0.4, 0.8), limits = c(0, 1)) +
#'     scale_x_continuous("Axis evolution SD",
#'                        breaks = c(0, 0.05, 0.1)) +
#'     scale_color_viridis_d("Population SD",
#'                           option = "A", end = 0.7) +
#'     theme(legend.position = "none",
#'           axis.title = element_text(size = 10),
#'           axis.text = element_text(size = 9),
#'           plot.margin = margin(0,l=8,r=8,t=8)) +
#'     NULL
#'
#' with(formals(quant_gen), r0 / a0)
#'
#'
#'
#'
#'
#'
#' coexist_spp_p <- ggarrange(coexist_d_spp_p,
#'                            coexist_stoch_spp_p +
#'                                theme(axis.title.y = element_blank(),
#'                                      axis.text.y = element_blank()),
#'                            nrow = 1, draw = FALSE,
#'                            labels = c("a", "b"),
#'                            label.args = label_args)
#'
#'
#' # coexist_spp_p
#'
#'
#'
#'
#'
#'
#'
#' if (.RESAVE_PLOTS) {
#'     save_plot(coexist_spp_p, 3, 4, .prefix = "3-")
#' }
#'
#'
#' # # (No effect at all, so no point in creating this.)
#' # inv_sims_one_p_fun("extinct",
#' #                    .fill_low = "#d01c8b",
#' #                    .fill_high = "#4dac26",
#' #                    .fill_limits = NULL)
#'
#'
#'
#'
#' # --------------*
#' # Fig S1 : sigma_i --> # spp ----
#' # --------------*
#'
#' add_var_effect_df <- grab_sims(.d = 0,
#'                                .eta = -1:1 * etas[[2]],
#'                                .add_var = c(0.01, 0.05, 0.1),
#'                                .sigma_N = 0,
#'                                .sigma_V = 0,
#'                                .vary_d2 = FALSE) %>%
#'     group_by(add_var, eta, rep) %>%
#'     summarize(n_spp = spp[N > 0] %>% unique() %>% length()) %>%
#'     ungroup() %>%
#'     mutate(eta = factor(eta, levels = -1:1 * etas[[2]],
#'                         labels = paste0(c("sub-", "", "super-"), "additive")))
#'
#' add_var_effect_p <- add_var_effect_df %>%
#'     mutate(n_spp = n_spp / 100) %>%
#'     ggplot(aes(add_var, n_spp)) +
#'     geom_hline(yintercept = 0, color = "gray80", linetype = 1) +
#'     geom_jitter(width = 0.002, height = 0, shape = 1, size = 1) +
#'     facet_wrap(~ eta, ncol = 1) +
#'     scale_color_viridis_d(expression(italic(sigma[i])^2),
#'                           begin = 0.1, end = 0.85, guide = FALSE) +
#'     scale_y_continuous("Proportion of species that coexist",
#'                        breaks = c(0, 0.4, 0.8), limits = c(0, 1)) +
#'     scale_x_continuous("Additive genetic variance")
#'
#'
#' if (.RESAVE_PLOTS) {
#'     save_plot(add_var_effect_p, 2.5, 4, .prefix = "S1-")
#' }
#'
#'
#'
#'
#'
#'
#'
#'
#' # --------------*
#' # Fig S2 - sigma_V - invader coexist ----
#' # --------------*
#'
#'
#' determ_inv_sims <- map_dfr(0:20,
#'                            function(i) {
#'                                fn <- rds(paste0("giant_inv_sims/giant_inv",
#'                                                 "_sims_outcomes_", i))
#'                                readRDS(fn) %>%
#'                                    filter(sigma_V == 0, sigma_N == 0,
#'                                           eta != 0) %>%
#'                                    select(-sigma_V, -sigma_N)
#'                            }) %>%
#'     bind_rows(map_dfr(0:6,
#'                       function(i) {
#'                           fn <- rds(paste0("giant_inv_sims_add-mid/",
#'                                            "giant_inv_sims_add-mid_",
#'                                            "outcomes_", i))
#'                           readRDS(fn) %>%
#'                               filter(sigma_V == 0, sigma_N == 0) %>%
#'                               select(-sigma_V, -sigma_N)
#'                       })) %>%
#'     mutate(d1 = factor(d1, levels = sort(unique(d1)),
#'                        labels = sprintf("d[1] == %.4f", sort(unique(d1)))),
#'            eta = factor(eta, levels = c(-0.6, 0, 0.6),
#'                         labels = paste0(c("'sub-", "'", "'super-"),
#'                                         "additive'"))) %>%
#'     select(-total)
#'
#' get_determ <- function(.par, .V1, .V2, .eta, .d1) {
#'     filter(determ_inv_sims, V1 == .V1[1], V2 == .V2[1],
#'            eta == .eta[1], d1 == .d1[1])[[.par]]
#' }
#'
#'
#' stoch_inv_sims <- map_dfr(0:20,
#'                           function(i) {
#'                               fn <- rds(paste0("giant_inv_sims/giant_inv",
#'                                                "_sims_outcomes_", i))
#'                               readRDS(fn) %>%
#'                                   filter((sigma_V == 0.1 & sigma_N == 0) |
#'                                              (sigma_V == 0 & sigma_N == 0.1)) %>%
#'                                   filter(eta != 0)
#'                           }) %>%
#'     bind_rows(map_dfr(0:6,
#'                       function(i) {
#'                           fn <- rds(paste0("giant_inv_sims_add-mid/",
#'                                            "giant_inv_sims_add-mid_",
#'                                            "outcomes_", i))
#'                           readRDS(fn) %>%
#'                               filter((sigma_V == 0.1 & sigma_N == 0) |
#'                                          (sigma_V == 0 & sigma_N == 0.1))
#'                       })) %>%
#'     mutate(d1 = factor(d1, levels = sort(unique(d1)),
#'                        labels = sprintf("d[1] == %.4f", sort(unique(d1)))),
#'            sigma_V = factor(sigma_V, levels = sort(unique(sigma_V)),
#'                             labels = sprintf("sigma[V] == %.4f",
#'                                              sort(unique(sigma_V)))),
#'            sigma_N = factor(sigma_N, levels = sort(unique(sigma_N)),
#'                             labels = sprintf("sigma[N] == %.4f",
#'                                              sort(unique(sigma_N)))),
#'            eta = factor(eta, levels = c(-0.6, 0, 0.6),
#'                         labels = paste0(c("'sub-", "'", "'super-"),
#'                                         "additive'"))) %>%
#'     mutate(across(coexist:extinct, ~ .x / total)) %>%
#'     select(-total) %>%
#'     group_by(V1, V2, eta, d1) %>%
#'     mutate(coexist_diff = coexist - get_determ("coexist", V1, V2, eta, d1),
#'            replace_diff = replace - get_determ("replace", V1, V2, eta, d1),
#'            reject_diff = reject - get_determ("reject", V1, V2, eta, d1),
#'            extinct_diff = extinct - get_determ("extinct", V1, V2, eta, d1)) %>%
#'     ungroup()
#'
#'
#'
#'
#'
#'
#' inv_sims_one_p_fun <- function(.par = "coexist",
#'                                .which_sigma = "V",
#'                                .fill_low = "#ca0020",
#'                                .fill_high = "#0571b0",
#'                                .fill_limits = c(-1, 1)) {
#'
#'     # .par = "coexist"; .which_sigma = "V"
#'     # .fill_low = "#ca0020"; .fill_high = "#0571b0"; .fill_limits = c(-1, 1)
#'     # rm(.par, .which_sigma, .fill_low, .fill_high, .fill_limits)
#'     # rm(fill_scale, make_eta_fct, midi, .sigma, ss_df)
#'
#'     .par <- match.arg(.par, c("coexist", "replace", "reject", "extinct"))
#'
#'     .which_sigma <- match.arg(.which_sigma, c("V", "N"))
#'
#'
#'     .fill_lab <- case_when(grepl("st$", .par) ~ paste0(.par, "ence"),
#'                            grepl("ct$", .par) ~ paste0(.par, "ion"),
#'                            grepl("ce$", .par) ~ paste0(.par, "ment"),
#'                            TRUE ~ "")
#'     .fill_lab <- bquote("Effect of" ~ sigma[.(.which_sigma)] ~
#'                             .(paste0("on ", .fill_lab, ":")))
#'
#'     fill_scale <- scale_fill_gradient2(.fill_lab,
#'                                        low = .fill_low,
#'                                        mid = "white",
#'                                        high = .fill_high,
#'                                        midpoint = 0,
#'                                        limits = .fill_limits,
#'                                        labels = function(x) {
#'                                            z <- case_when(
#'                                                x > 0 ~ sprintf("{}+%.1f",x),
#'                                                x == 0 ~ "{}~0.0",
#'                                                TRUE ~ sprintf("%.1f",x))
#'                                            parse(text = z)
#'                                        })
#'
#'     make_eta_fct <- function(..x) {
#'         factor(..x, levels = c(-0.6, 0, 0.6),
#'                labels = paste0(c("'sub-", "'", "'super-"),
#'                                "additive'"))
#'     }
#'     midi <- ceiling(formals(stable_points)[["line_n"]] / 2)
#'
#'     .sigma <- sprintf("sigma_%s", .which_sigma)
#'
#'     ss_df <- map_dfr(c(-0.6, 0, 0.6), ~ mutate(stable_points(.x),
#'                                                eta = make_eta_fct(.x)))
#'
#'     # p <-
#'     stoch_inv_sims %>%
#'         filter(!!as.name(.sigma) == sprintf("sigma[%s] == %.4f",
#'                                             .which_sigma, 0.1)) %>%
#'         ggplot(aes(V1, V2)) +
#'         geom_raster(aes_string(fill = paste0(.par, "_diff"))) +
#'         geom_tile(data = determ_inv_sims %>% filter(!!as.name(.par) > 0),
#'                   fill = NA, color = "black", size = 0.1) +
#'         facet_grid(d1 ~ eta, label = label_parsed) +
#'         geom_point(data = ss_df %>%
#'                        filter(eta == "'super-additive'", V1 == 0),
#'                    size = 2, shape = 16, color = "black") +
#'         geom_point(data = ss_df %>%
#'                        filter(eta == "'super-additive'", V2 == 0),
#'                    size = 2, shape = 21, color = "black", fill = "white") +
#'         geom_point(data = ss_df %>%
#'                        filter(eta == "'sub-additive'"),
#'                    size = 2, shape = 16, color = "black") +
#'         geom_path(data = ss_df %>%
#'                       filter(eta == "'additive'"),
#'                   color = "gray50", size = 1) +
#'         geom_point(data = tibble(V1 = sqrt(2), V2 = sqrt(2),
#'                                  eta = make_eta_fct(0)),
#'                    shape = 16, size = 2, color = "black") +
#'         scale_x_continuous("Conflicting axis", breaks = c(0, 2, 4)) +
#'         scale_y_continuous("Ameliorative axis", breaks = c(0, 2, 4)) +
#'         coord_equal() +
#'         NULL +
#'         theme(strip.text = element_text(size = 8),
#'               strip.text.y = element_text(angle = 0, hjust = 0,
#'                                           margin = margin(0,0,0,l=3)),
#'               strip.text.x = element_text(margin = margin(0,0,0,b=3)),
#'               panel.border = element_rect(size = 0.5, fill = NA),
#'               plot.title = element_text(size = 12, hjust = 0.5,
#'                                         margin = margin(0,0,0,b=6)),
#'               legend.title = element_text(size = 10),
#'               legend.text.align = 0.5,
#'               plot.margin = margin(0,0,0,0),
#'               legend.position = "bottom") +
#'         guides(fill = guide_colorbar(title.position = "top")) +
#'         fill_scale
#'
#' }
#'
#' inv_sims_V_coexist_p <- inv_sims_one_p_fun("coexist")
#'
#'
#'
#'
#' if (.RESAVE_PLOTS) {
#'     save_plot(inv_sims_V_coexist_p,
#'               4, 4, .name = "S2-sigmaV_coexist")
#' }
#'
#'
#' # --------------*
#' # Fig S3 - sigma_V - invader replace ----
#' # --------------*
#'
#' if (.RESAVE_PLOTS) {
#'     save_plot(inv_sims_one_p_fun("replace",
#'                                  .fill_low = "#e66101",
#'                                  .fill_high = "#5e3c99"),
#'               4, 4, .name = "S3-sigmaV_replace")
#' }
#'
#'
#' # --------------*
#' # Fig S4 - sigma_N - invader coexist ----
#' # --------------*
#'
#'
#' if (.RESAVE_PLOTS) {
#'     save_plot(inv_sims_one_p_fun("coexist", .which_sigma = "N"),
#'               4, 4, .name = "S4-sigmaN_coexist")
#' }
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' # =============================================================================*
#' # =============================================================================*
#'
#' # Fig 4: Conditional coexistence ----
#'
#' # =============================================================================*
#' # =============================================================================*
#'
#'
#'
#'
#' cond_coexist_test <- function(.V0, .eta_sign, .d1, .sigma_V = 0, .sigma_N = 0) {
#'
#'     # .V0 = "restricted"; .eta_sign = 0; .d1 = -0.1; .sigma_V = 0; .sigma_N = 0
#'     # rm(.V0, .eta_sign, .d1, .sigma_V, .sigma_N)
#'
#'     .dist <- 0.6
#'     which_switch <- 3 # which spp to switch for unrestricted
#'
#'     .n <- 5
#'     .q = 2
#'     stopifnot(is.numeric(.d1) && length(.d1) == 1)
#'     .ds <- c(.d1, abs(.d1))
#'     stopifnot(is.numeric(.eta_sign) && length(.eta_sign) == 1 &&
#'                   .eta_sign %in% -1:1)
#'     .V0 <- match.arg(.V0, c("restricted", "unrestricted"))
#'     .lab <- .V0
#'
#'     .eta <- .eta_sign * etas[[2]]
#'
#'     if (.lab == "restricted" && .eta < 0) {
#'         stop(paste("\n'restricted' with sub-additivity is not programmed bc",
#'                    "there's only one stable axis state.",
#'                    "Thus, there is no way to restrict starting axes to",
#'                    "be outside of that state's basin of attraction."))
#'     }
#'
#'     if (.eta < 0) {
#'         .V0 <- rbind(seq(3, 2, length.out = 5),
#'                      seq(2, 3, length.out = 5))
#'     } else if (.eta == 0) {
#'         .V0 <- rbind(seq(1.5, 0, length.out = 5),
#'                      seq(1.6, 3.1, length.out = 5))
#'         if (.lab == "unrestricted") {
#'             v2 <- .V0[2, which_switch]
#'             .V0[2, which_switch] <- .V0[1, which_switch]
#'             .V0[1, which_switch] <- v2
#'         }
#'     } else {
#'         .V0 <- rbind(seq(1.2, 0, length.out = 5),
#'                      seq(1.3, 2.5, length.out = 5))
#'         if (.lab == "unrestricted") {
#'             v2 <- .V0[2, which_switch]
#'             .V0[2, which_switch] <- .V0[1, which_switch]
#'             .V0[1, which_switch] <- v2
#'         }
#'     }
#'     .V0 <- round(.V0, 3)
#'
#'     if (.sigma_V == 0 && .sigma_N == 0) {
#'         .nreps <- 1
#'         .N_THREADS <- 1
#'     } else .nreps <- 12
#'
#'
#'     qg <- quant_gen(q = .q, eta = .eta, d = .ds,
#'                    n_reps = .nreps, n = ncol(.V0),
#'                    V0 = .V0,
#'                    sigma_V0 = 0,
#'                    sigma_V = .sigma_V,
#'                    sigma_N = .sigma_N,
#'                    spp_gap_t = 500L,
#'                    final_t = 20e3L,
#'                    add_var = rep(0.05, .n),
#'                    n_threads = .N_THREADS,
#'                    show_progress = FALSE)
#'
#'
#'     if (.sigma_V == 0) {
#'
#'         out <- list(nv = qg %>%
#'                         .[["nv"]] %>%
#'                         pivot())
#'         if (.sigma_N == 0) {
#'             out$nv <- select(out$nv, -rep)
#'
#'             qg$call[["q"]] <- eval(.q)
#'             qg$call[["n"]] <- eval(.n)
#'             qg$call[["d"]] <- eval(.ds)
#'             qg$call[["eta"]] <- eval(.eta)
#'             qg$call[["add_var"]] <- eval(rep(0.05, .n))
#'             out$jacs <- jacobians(qg)
#'         }
#'
#'
#'     } else {
#'
#'         out <- list(nv = qg %>%
#'                         .[["nv"]] %>%
#'                         mutate(axis = paste0("V", axis)) %>%
#'                         nest(value = c(geno, pheno)) %>%
#'                         spread(axis, value) %>%
#'                         unnest(c(V1, V2), names_sep = "_") %>%
#'                         rename(V1 = V1_geno,
#'                                V2 = V2_geno,
#'                                Vp1 = V1_pheno,
#'                                Vp2 = V2_pheno))
#'
#'     }
#'
#'     out$nv <- out$nv %>%
#'         mutate(V0 = .lab, eta = .eta, d1 = .d1,
#'                sigma_V = .sigma_V, sigma_N = .sigma_N)
#'
#'     return(out)
#' }
#'
#'
#' if (.REDO_SIMS) {
#'
#'     # Takes ~6 sec
#'     cond_coexist <- crossing(.V0 = c("restricted", "unrestricted"),
#'                              .eta_sign = -1:1,
#'                              .d1 = c(-0.1, 0.1)) %>%
#'         filter(!(.V0 == "restricted" & .eta_sign < 0)) %>%
#'         pmap(cond_coexist_test)
#'     saveRDS(cond_coexist, rds("cond_coexist"))
#'
#'     # Takes ~3.5 min
#'     set.seed(351367879)
#'     cond_coexist_sV <- crossing(.V0 = c("restricted", "unrestricted"),
#'                              .eta_sign = -1:1,
#'                              .d1 = c(-0.1, 0.1),
#'                              .sigma_V = c(0.05, 0.1, 0.2)) %>%
#'         filter(!(.V0 == "restricted" & .eta_sign < 0)) %>%
#'         pmap(cond_coexist_test)
#'     saveRDS(cond_coexist_sV, rds("cond_coexist_sV"))
#'
#'     # Takes ~1.7 min
#'     set.seed(602553504)
#'     cond_coexist_sN <- crossing(.V0 = c("restricted", "unrestricted"),
#'                              .eta_sign = -1:1,
#'                              .d1 = c(-0.1, 0.1),
#'                              .sigma_N = c(0.05, 0.1, 0.2)) %>%
#'         filter(!(.V0 == "restricted" & .eta_sign < 0)) %>%
#'         pmap(cond_coexist_test)
#'     saveRDS(cond_coexist_sN, rds("cond_coexist_sN"))
#'
#' } else {
#'     cond_coexist <- readRDS(rds("cond_coexist"))
#'     cond_coexist_sV <- readRDS(rds("cond_coexist_sV"))
#'     cond_coexist_sN <- readRDS(rds("cond_coexist_sN"))
#' }
#'
#'
#' cond_coexist_df_prep <- function(.dd) {
#'     .dd %>%
#'         map_dfr(~ .x$nv) %>%
#'         # filter(time < 4e3L) %>%
#'         mutate(V0 = factor(V0, levels = c("restricted", "unrestricted")),
#'                eta = factor(eta, levels = -1:1 * etas[[2]],
#'                             labels = c("sub-additive", "additive",
#'                                        "super-additive")),
#'                d1 = factor(d1, levels = c(-0.1, 0.1),
#'                            labels = c("conflicting",
#'                                       "ameliorative"))) %>%
#'         # `axis_space` is a combination of starting axis values and eta:
#'         mutate(axis_space =
#'                    case_when(V0 == "unrestricted" & eta == "sub-additive" ~
#'                                  "i",
#'                              V0 == "restricted" & eta == "super-additive" ~
#'                                  "ii",
#'                              V0 == "restricted" & eta == "additive" ~
#'                                  "iii",
#'                              V0 == "unrestricted" & eta == "super-additive" ~
#'                                  "iv",
#'                              V0 == "unrestricted" & eta == "additive" ~
#'                                  "v",
#'                              TRUE ~ NA_character_) %>%
#'                    factor(levels = c("i", "ii", "iii", "iv", "v")))
#' }
#'
#'
#'
#' cond_coexist_df <- cond_coexist %>%
#'     cond_coexist_df_prep() %>%
#'     select(-starts_with("sigma_"))
#'
#' # Reps with stochasticity:
#' cond_coexist_stoch_df <- map_dfr(list(cond_coexist_sV,
#'                                        cond_coexist_sN),
#'                                   cond_coexist_df_prep)
#'
#'
#' if (any(is.na(cond_coexist_df$axis_space))) {
#'     stop("\nERROR: unknown combination of V0 and eta in `cond_coexist_df`")
#' }
#' if (any(is.na(cond_coexist_stoch_df$axis_space))) {
#'     stop("\nERROR: unknown combination of V0 and eta in `cond_coexist_stoch_df`")
#' }
#'
#'
#'
#' # Starting conditions and trajectories:
#'
#'
#' cond_coexist_sc_p_fun <- function(.d1) {
#'     .dd <- cond_coexist_df %>%
#'         filter(d1 == .d1) %>%
#'         split(interaction(.$axis_space, .$spp, drop = TRUE)) %>%
#'         map_dfr(~ mutate(.x,
#'                          first = time == min(time),
#'                          last = time == max(time))) %>%
#'         select(axis_space, time, spp, V1, V2, first, last) %>%
#'         arrange(axis_space, time)
#'     .dd %>%
#'         ggplot(aes(V1, V2)) +
#'         geom_abline(data = tibble(axis_space = factor(c("ii", "iv")),
#'                                   slp = 1, int = 0),
#'                     aes(slope = slp, intercept = int), linetype = 3, color = "gray70") +
#'         geom_line(data = bind_rows(stable_points(0), stable_points(0)) %>%
#'                       mutate(axis_space = factor(rep(c("iii", "v"),
#'                                                       each = stable_points %>%
#'                                                           formals() %>%
#'                                                           .[["line_n"]]))),
#'                   linetype = 2, color = "black") +
#'         geom_point(data = .dd %>% filter(first),
#'                    aes(color = spp), size = 1.5, na.rm = TRUE) +
#'         geom_point(data = map_dfr(etas[[2]] * c(-1, 1, 1), ~ stable_points(.x)) %>%
#'                        mutate(axis_space = map2(c("i", "ii", "iv"), c(1,2,2), rep) %>%
#'                                   do.call(what = c) %>%
#'                                   factor(),
#'                               shp = case_when(V1 > 0 & V2 > 0 ~ 3L,
#'                                               V1 == 0 ~ 2L,
#'                                               TRUE ~ 1L) %>%
#'                                   factor(levels = 1:3)),
#'                    aes(shape = shp), size = 3, color = "black") +
#'         geom_point(data = .dd %>%
#'                        filter(time < max(time), last),
#'                    aes(color = spp), shape = 4, size = 1.5) +
#'         geom_path(aes(color = spp)) +
#'         scale_color_brewer(palette = "Dark2", guide = FALSE) +
#'         scale_shape_manual(values = c(5,1,2), guide = FALSE) +
#'         scale_size_continuous(range = c(0.1, 1)) +
#'         facet_grid(~ axis_space) +
#'         coord_equal(xlim = c(-0.1, 3.15), ylim = c(-0.1, 3.15)) +
#'         ylab("Axis 2") +
#'         xlab("Axis 1") +
#'         theme(plot.margin = margin(0,0,0,b=6),
#'               plot.title = element_text(hjust = 0.5,
#'                                         margin = margin(0,0,0,b=6)))
#' }
#'
#' cond_coexist_sc_p <- cond_coexist_sc_p_fun("conflicting")
#'
#'
#'
#' # Time series of abundances
#' cc_N_p <- cond_coexist_df %>%
#'     filter(d1 == "conflicting") %>%
#'     filter(time < 7e3) %>%
#'     ggplot(aes(time / 1000L, N, color = spp)) +
#'     geom_line() +
#'     geom_point(data = cond_coexist_df %>%
#'                    filter(d1 == "conflicting") %>%
#'                    group_by(axis_space, spp) %>%
#'                    filter(time == max(time)) %>%
#'                    ungroup() %>%
#'                    filter(time < max(time)),
#'                shape = 4, size = 1.5) +
#'     facet_wrap(~ axis_space, nrow = 1) +
#'     scale_color_brewer(palette = "Dark2",
#'                        guide = FALSE) +
#'     scale_y_continuous("Abundance", trans = "log",
#'                        breaks = 10^(c(-3, 0, 3)),
#'                        labels = parse(
#'                            text = sprintf("10^{%i}",
#'                                           c(-3, 0, 3)))) +
#'     xlab("Time (× 1,000)") +
#'     theme(plot.margin = margin(0,r=12,t=10,b=10),
#'           axis.title.x = element_blank(),
#'           plot.title = element_text(hjust = 0.5,
#'                                     margin = margin(0,0,0,b=6)))
#'
#'
#'
#'
#'
#' stable_state_df <- map_dfr(c(1:2, 4),
#'                            function(i) {
#'                                ts <- cond_coexist_df$axis_space %>%
#'                                    unique() %>%
#'                                    sort() %>%
#'                                    .[i]
#'                                stable_points((etas[[2]] * c(-1,1,0,1,0))[i]) %>%
#'                                    mutate(axis_space = ts)
#'                            }) %>%
#'     mutate(time = max(cond_coexist_df$time[cond_coexist_df$time < 7e3]),
#'            shp = case_when(V1 > 0 & V2 > 0 ~ 3L,
#'                            V1 == 0 ~ 2L,
#'                            TRUE ~ 1L) %>%
#'                factor(levels = 1:3)) %>%
#'     gather(axis, value, V1:V2) %>%
#'     mutate(axis = gsub("^V", "axis ", axis))
#'
#'
#'
#' dist_from_equil <- function(V1, V2, eta) {
#'     .eta <- case_when(eta[1] == "additive" ~ 0,
#'                       eta[1] == "sub-additive" ~ -etas[[2]],
#'                       TRUE ~ etas[[2]])
#'     if (.eta == 0) {
#'         eq_radius <- with(formals(quant_gen), sqrt((r0 / f) - 1))
#'         R <- sqrt(V1^2 + V2^2)
#'         dist <- abs(R - eq_radius)
#'     } else {
#'         equil <- stable_points(.eta)
#'         dist <- matrix(NA_real_, length(V1), nrow(equil))
#'         for (j in 1:nrow(equil)) {
#'             dist[,j] <- sqrt((V1 - equil$V1[j])^2 + (V2 - equil$V2[j])^2)
#'         }
#'         dist <- apply(dist, 1, min)
#'     }
#'     return(dist)
#' }
#'
#'
#'
#' cc_V_p <- cond_coexist_df %>%
#'     filter(d1 == "conflicting") %>%
#'     filter(time < 7e3) %>%
#'     split(.$eta) %>%
#'     map_dfr(function(.dd) {
#'         mutate(.dd, dist = dist_from_equil(V1, V2, eta))
#'     }) %>%
#'     mutate(dist = ifelse(eta == "additive", dist, mean(dist))) %>%
#'     gather(axis, value, V1:V2) %>%
#'     mutate(axis = gsub("^V", "axis ", axis)) %>%
#'     ggplot(aes(time / 1000L, value)) +
#'     geom_hline(yintercept = 0, size = 0.5,
#'                linetype = 1, color = "gray70") +
#'     geom_vline(xintercept = 0, size = 0.5,
#'                linetype = 1, color = "gray70") +
#'     geom_point(data = stable_state_df,
#'                aes(shape = shp), size = 3) +
#'     geom_line(aes(color = spp, size = dist)) +
#'     geom_line(aes(color = spp)) +
#'     geom_point(data = cond_coexist_df %>%
#'                    filter(d1 == "conflicting") %>%
#'                    gather(axis, value, V1:V2) %>%
#'                    mutate(axis = gsub("^V", "axis ", axis)) %>%
#'                    group_by(axis_space, spp) %>%
#'                    filter(time == max(time)) %>%
#'                    ungroup() %>%
#'                    filter(time < max(time)),
#'                aes(color = spp), shape = 4, size = 1.5) +
#'     facet_grid(axis ~ axis_space) +
#'     scale_color_brewer(palette = "Dark2", guide = FALSE) +
#'     scale_shape_manual(values = c(5,1,2), guide = FALSE) +
#'     scale_size_continuous(range = c(0.4, 2), guide = FALSE) +
#'     scale_y_continuous("Axis value", limits = c(-0.2, NA)) +
#'     scale_x_continuous("Time (× 1,000)",
#'                        limits = c(0, 7.2)) +
#'     theme(plot.margin = margin(0,0,0,t=10),
#'           plot.title = element_text(hjust = 0.5,
#'                                     margin = margin(0,0,0,b=6)))
#'
#'
#'
#'
#' cond_coexist_p <- plot_grid(cond_coexist_sc_p +
#'                                 xlab(sprintf("Axis 1\n(%s)", "conflicting")) +
#'                                 ylab("(ameliorative)\nAxis 2") +
#'                                 ggtitle(paste("Starting conditions and",
#'                                               "trajectories")) +
#'                                 theme(plot.margin = margin(0,0,0,r=12)),
#'                             cc_N_p +
#'                                 ggtitle("Abundance time series"),
#'                             cc_V_p +
#'                                 ggtitle("Axis time series"),
#'                             align = "v", axis = "l",
#'                             ncol = 1, rel_heights = c(3, 2, 4),
#'                             labels = LETTERS[1:3],
#'                             label_fontface = "plain", label_x = 0.06)
#'
#'
#'
#' if (.RESAVE_PLOTS) {
#'     save_plot(cond_coexist_p, 6, 7, .name = "4-cond_coexist")
#' }
#'
#'
#'
#'
#' #'
#' #' They're all stable except for sub-additive and conflicting axis 1,
#' #' which is neutrally stable (see `_main-results__stability.R`).
#' #'
#'
#'
#'
#'
#' # =============================================================================*
#' # =============================================================================*
#'
#' # <not updated> Figs 5,S2-S3: Conditional coexistence ----
#'
#' # =============================================================================*
#' # =============================================================================*
#'
#'
#'
#'
#' cc_N_stoch_plot_fun <- function(.V_stoch, .ts = FALSE) {
#'
#'     .d1 = "conflicting"
#'
#'     # .V_stoch = TRUE; .ts = FALSE
#'     # rm(.d1, .V_stoch, .ts, .dd, .dd2)
#'
#'     stopifnot(is.logical(.V_stoch) && length(.V_stoch) == 1)
#'
#'     if (.V_stoch) {
#'         .sigma_V <- 0.1
#'         .sigma_N <- 0
#'     } else {
#'         .sigma_V <- 0
#'         .sigma_N <- 0.1
#'     }
#'
#'     stopifnot((.sigma_V > 0 || .sigma_N) > 0 && !(.sigma_V > 0 && .sigma_N > 0))
#'
#'     .dd <- cond_coexist_stoch_df %>%
#'         filter(d1 == .d1,
#'                sigma_N == .sigma_N,
#'                sigma_V == .sigma_V)
#'
#'     if (.ts) {
#'         .dd2 <- cond_coexist_df %>%
#'             filter(d1 == .d1) %>%
#'             mutate(rep = 0L) %>%
#'             select(axis_space, rep, time, spp, N)
#'         .dd %>%
#'             mutate(rep = rep %>% paste() %>% as.integer()) %>%
#'             bind_rows(.dd2) %>%
#'             mutate(rep = factor(rep, levels = 0:12)) %>%
#'             mutate(id = interaction(spp, rep)) %>%
#'             ggplot(aes(time / 1000L, N, color = spp)) +
#'             geom_line(aes(group = id)) +
#'             geom_point(data = .dd %>%
#'                            group_by(axis_space, rep, spp) %>%
#'                            filter(time == max(time)) %>%
#'                            ungroup() %>%
#'                            filter(time < max(time)),
#'                        shape = 4, size = 1.5) +
#'             facet_grid(rep ~ axis_space) +
#'             scale_color_brewer(palette = "Dark2",
#'                                guide = FALSE) +
#'             scale_y_continuous("Abundance", trans = "log",
#'                                breaks = 10^(c(-3, 0, 3)),
#'                                labels = parse(
#'                                    text = sprintf("10^{%i}",
#'                                                   c(-3, 0, 3)))) +
#'             xlab("Time (× 1,000)") +
#'             theme(strip.text.y = element_blank())
#'     } else {
#'
#'         .dd <- .dd %>%
#'             filter(time == max(time)) %>%
#'             select(-time)
#'
#'         .dd2 <- cond_coexist_df %>%
#'             filter(d1 == .d1, time == max(time)) %>%
#'             group_by(axis_space) %>%
#'             mutate(tN = sqrt(N) / sum(sqrt(N))) %>%
#'             ungroup() %>%
#'             mutate(rep = 0L) %>%
#'             select(axis_space, rep, spp, N, tN)
#'
#'         .dd3 <- map_dfr(list(.dd, .dd2),
#'                         ~ .x %>%
#'                             group_by(axis_space, rep) %>%
#'                             summarize(n_spp = sum(N > 0)) %>%
#'                             ungroup() %>%
#'                             mutate(rep = factor(rep %>% paste() %>% as.integer(),
#'                                                 levels = 0:12)))
#'
#'         .dd %>%
#'             group_by(axis_space, rep) %>%
#'             mutate(tN = sqrt(N) / sum(sqrt(N))) %>%
#'             # mutate(tN = N / sum(N)) %>%
#'             ungroup() %>%
#'             select(axis_space, rep, spp, N, tN) %>%
#'             mutate(rep = rep %>% paste() %>% as.integer()) %>%
#'             bind_rows(.dd2) %>%
#'             mutate(rep = factor(rep, levels = 0:12)) %>%
#'             ggplot(aes(rep)) +
#'             geom_bar(aes(weight = tN, fill = spp), color = "white", size = 0.1) +
#'             geom_text(data = .dd3 %>% filter(n_spp > 1),
#'                       aes(label = n_spp, y = 1.05),
#'                       size = 8 / 2.83465) +
#'             geom_vline(xintercept = 1.5, size = 0.5) +
#'             facet_grid( ~ axis_space) +
#'             scale_fill_brewer(palette = "Dark2", guide = FALSE) +
#'             scale_y_continuous("Scaled relative abundance",
#'                                breaks = c(0, 0.5, 1)) +
#'             theme(axis.text.x = element_blank(),
#'                   axis.ticks.x = element_blank(),
#'                   axis.title.x = element_blank(),
#'                   plot.title = element_text(hjust = 0.5,
#'                                             margin = margin(0,0,0,b=6),
#'                                             size = 11)) +
#'             NULL
#'     }
#' }
#'
#'
#' cond_coexist_stoch_ps <- map(c(TRUE, FALSE), cc_N_stoch_plot_fun)
#' names(cond_coexist_stoch_ps) <- c("V_stoch", "N_stoch")
#'
#' cond_coexist_stoch_p <- ggarrange(cond_coexist_stoch_ps[["N_stoch"]] +
#'                                       ggtitle("With abundance stochasticity") +
#'                                       theme(plot.margin = margin(0,0,0,b=12),
#'                                             axis.title.y = element_blank()),
#'                                   cond_coexist_stoch_ps[["V_stoch"]] +
#'                                       ggtitle("With axis stochasticity") +
#'                                       theme(plot.margin = margin(0,0,0,0),
#'                                             axis.title.y = element_blank()),
#'                                   ncol = 1,
#'                                   padding = unit(1, "lines"),
#'                                   left = "Scaled relative abundance",
#'                                   draw = FALSE, labels = LETTERS[1:2],
#'                                   label.args = list(gp = grid::gpar(
#'                                       fontface = "plain", fontsize = 12),
#'                                       hjust = -2))
#'
#' # cond_coexist_stoch_p
#'
#' #'
#' #' This helps explain why you get totally different outcomes for situation v
#' #' when sigma_V > 0, compared to when sigma_V = 0.
#' #'
#' #' This is NOT explained by the distribution for V stochasticity
#' #' being lognormal and having a mean > 1.
#' #' I tried with simulations where the distribution was corrected for this
#' #' (and did indeed have a mean of 1).
#' #' The results were the same.
#' #'
#' cc_sigmaV_sit_v_p <- cond_coexist_stoch_df %>%
#'     filter(d1 == "conflicting",
#'            sigma_N == 0,
#'            sigma_V == 0.1) %>%
#'     filter(axis_space == "v") %>%
#'     mutate(id = interaction(spp, rep)) %>%
#'     arrange(time) %>%
#'     ggplot(aes(V1, V2, color = spp)) +
#'     geom_path(aes(group = id)) +
#'     geom_point(data = cond_coexist_stoch_df %>%
#'                    filter(d1 == "conflicting",
#'                           sigma_N == 0,
#'                           sigma_V == 0.1) %>%
#'                    filter(axis_space == "v") %>%
#'                    group_by(spp) %>%
#'                    summarize(V1 = V1[time == min(time)][1],
#'                              V2 = V2[time == min(time)][1])) +
#'     facet_wrap(~ rep, nrow = 3) +
#'     coord_equal(xlim = c(0, 3.11), ylim = c(0, 3.11)) +
#'     scale_color_brewer(palette = "Dark2", guide = FALSE) +
#'     ylab("(ameliorative)\nAxis 2") +
#'     xlab("Axis 1\n(conflicting)")
#'
#'
#'
#'
#'
#' qg <- quant_gen(q = 2, eta = 0, d = c(-0.1, 0.1),
#'                 n_reps = 24, n = 10,
#'                 sigma_V = 0.1,
#'                 spp_gap_t = 500L,
#'                 final_t = 20e3L,
#'                 add_var = rep(0.05, 10),
#'                 n_threads = 3,
#'                 show_progress = TRUE)
#'
#' nv <- qg %>%
#'     .[["nv"]] %>%
#'     mutate(axis = paste0("V", axis)) %>%
#'     nest(value = c(geno, pheno)) %>%
#'     spread(axis, value) %>%
#'     unnest(c(V1, V2), names_sep = "_") %>%
#'     rename(V1 = V1_geno,
#'            V2 = V2_geno,
#'            Vp1 = V1_pheno,
#'            Vp2 = V2_pheno)
#'
#' nv %>%
#'     mutate(id = interaction(spp, rep)) %>%
#'     arrange(time) %>%
#'     ggplot(aes(V1, V2, color = spp)) +
#'     # ggplot(aes(Vp1, Vp2, color = spp)) +
#'     geom_path(aes(group = id)) +
#'     geom_point(data = nv %>%
#'                    filter(!is.na(spp)) %>%
#'                    group_by(spp) %>%
#'                    summarize(V1 = V1[time == min(time)][1],
#'                              V2 = V2[time == min(time)][1])) +
#'                    # summarize(Vp1 = Vp1[time == min(time)][1],
#'                    #           Vp2 = Vp2[time == min(time)][1])) +
#'     facet_wrap(~ rep, nrow = 4) +
#'     coord_equal(xlim = c(0, 3.11), ylim = c(0, 3.11)) +
#'     scale_color_viridis_d(option = "D", guide = FALSE) +
#'     ylab("(ameliorative)\nAxis 2") +
#'     xlab("Axis 1\n(conflicting)")
#'
#' nv %>%
#'     mutate(id = interaction(spp, rep)) %>%
#'     ggplot(aes(time / 1000, N, color = spp)) +
#'     geom_path(aes(group = id)) +
#'     facet_wrap(~ rep, nrow = 4) +
#'     scale_color_viridis_d(option = "D", guide = FALSE) +
#'     ylab("Abundance") +
#'     xlab("Time (× 1000)")
#'
#'
#'
#' Z <- quant_gen(q = 2, eta = 0, d = c(-0.1, 0.1),
#'                n_reps = 24, n = 1,
#'                sigma_V = c(0.1, 0.025),
#'                final_t = 50e3L,
#'                # save_every = 10,
#'                save_every = 0,
#'                add_var = rep(0.05, 1),
#'                n_threads = 3)
#' znv <- Z$nv %>%
#'     mutate(axis = paste0("V", axis)) %>%
#'     nest(value = c(geno, pheno)) %>%
#'     spread(axis, value) %>%
#'     unnest(c(V1, V2), names_sep = "_") %>%
#'     rename(V1 = V1_geno,
#'            V2 = V2_geno,
#'            Vp1 = V1_pheno,
#'            Vp2 = V2_pheno)
#' znv %>%
#'     # filter(time == max(time)) %>%
#'     filter(!is.na(V1)) %>%
#'     ggplot(aes(V1, V2)) +
#'     geom_point(aes(color = spp)) +
#'     coord_equal(xlim = c(0, 3.11), ylim = c(0, 3.11)) +
#'     scale_color_viridis_d(option = "D", guide = FALSE) +
#'     ylab("(ameliorative)\nAxis 2") +
#'     xlab("Axis 1\n(conflicting)")
#'
#' # znv %>%
#' #     mutate(id = interaction(spp, rep)) %>%
#' #     arrange(time) %>%
#' #     ggplot(aes(V1, V2, color = spp)) +
#' #     # ggplot(aes(Vp1, Vp2, color = spp)) +
#' #     geom_path(aes(group = id)) +
#' #     geom_point(data = nv %>%
#' #                    filter(!is.na(spp)) %>%
#' #                    group_by(spp) %>%
#' #                    summarize(V1 = V1[time == min(time)][1],
#' #                              V2 = V2[time == min(time)][1])) +
#' #     # summarize(Vp1 = Vp1[time == min(time)][1],
#' #     #           Vp2 = Vp2[time == min(time)][1])) +
#' #     facet_wrap(~ rep, nrow = 4) +
#' #     coord_equal(xlim = c(0, 3.11), ylim = c(0, 3.11)) +
#' #     scale_color_viridis_d(option = "D", guide = FALSE) +
#' #     ylab("(ameliorative)\nAxis 2") +
#' #     xlab("Axis 1\n(conflicting)")
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' if (.RESAVE_PLOTS) {
#'     save_plot(cond_coexist_stoch_p, 6, 4, .prefix = "S1-")
#'     save_plot(cc_sigmaV_sit_v_p, 6, 5, .prefix = "S2-")
#' }
#'
#'
#'
#'
#'
#'
#'
#'
#' # =============================================================================*
#' # =============================================================================*
#'
#' # <not updated> Fig S3 Stoch. - # species ----
#'
#' # =============================================================================*
#' # =============================================================================*
#'
#'
#' # Simulations varying both d values
#'
#'
#' stoch_coexist_spp_df <- grab_sims(.d = seq(-0.25, 2, length.out = 10),
#'                                   .eta = -1:1 * etas[[2]],
#'                                   .add_var = c(0.01, 0.05, 0.1),
#'                                   .sigma_N = c(0, 0.05, 0.1, 0.2, 0.3),
#'                                   .sigma_V = c(0, 0.05, 0.1, 0.2, 0.3),
#'                                   .vary_d2 = TRUE) %>%
#'     group_by(d, eta, add_var, sigma_N, sigma_V, rep) %>%
#'     summarize(n_spp = spp[N > 0] %>% unique() %>% length()) %>%
#'     ungroup() %>%
#'     mutate(eta = factor(eta, levels = -1:1 * etas[[2]],
#'                         labels = c("sub-additive", "additive", "super-additive")))
#'
#'
#' stoch_coexist_p_fun <- function(.x) {
#'     stoch_coexist_spp_df %>%
#'         filter(eta == .x) %>%
#'         filter(sigma_N %in% c(0, 0.1, 0.2),
#'                sigma_V %in% c(0, 0.1, 0.2),
#'                add_var == 0.05) %>%
#'         mutate(sigma_N = factor(sigma_N, levels = sort(unique(sigma_N)),
#'                                 labels = sprintf("sigma[N] == %.2f",
#'                                                  sort(unique(sigma_N)))),
#'                sigma_V = factor(sigma_V, levels = sort(unique(sigma_V)),
#'                                 labels = sprintf("sigma[V] == %.2f",
#'                                                  sort(unique(sigma_V)))),
#'                add_var = factor(add_var, levels = sort(unique(add_var))),
#'                extinct = factor(n_spp == 0)) %>%
#'         mutate(n_spp = n_spp / 100) %>%
#'         # ggplot(aes(d, n_spp, color = add_var)) +
#'         ggplot(aes(d, n_spp)) +
#'         geom_hline(yintercept = 0, color = "gray80") +
#'         geom_vline(xintercept = 0, color = "gray80") +
#'         geom_jitter(aes(shape = extinct, size = extinct),
#'                     color = "dodgerblue",
#'                     width = 0.02, height = 0) +
#'         ggtitle(.x) +
#'         facet_grid(sigma_N ~ sigma_V, labeller = label_parsed) +
#'         # scale_color_viridis_d(expression(italic(sigma[i])^2),
#'         #                       begin = 0.1, end = 0.85) +
#'         scale_shape_manual(values = c(1, 4), guide = FALSE) +
#'         scale_size_manual(values = c(0.5, 2), guide = FALSE) +
#'         scale_y_continuous("Proportion of species that survive",
#'                            breaks = c(0, 0.4, 0.8), limits = c(0, 1)) +
#'         scale_x_continuous(expression(italic(d[1]) ~ "and" ~ italic(d[2])),
#'                            breaks = c(0, 1, 2)) +
#'         guides(color = guide_legend(override.aes = list(size = 2, shape = 19))) +
#'         theme(strip.text = element_text(size = 10),
#'               strip.text.y = element_text(angle = 0),
#'               plot.title = element_text(hjust = 0.5,
#'                                         margin = margin(0,0,0,b=12))) +
#'         NULL
#' }
#'
#' stoch_coexist_ps <- map(c("sub-additive", "additive", "super-additive"),
#'                         stoch_coexist_p_fun)
#' names(stoch_coexist_ps) <- c("sub", "add", "super")
#'
#' # stoch_coexist_ps[["sub"]]
#' # stoch_coexist_ps[["add"]]
#' # stoch_coexist_ps[["super"]]
#'
#'
#'
#' # Now looking at it near the boundaries, only varying d1:
#'
#'
#' stoch_vary_d1_coexist_spp_df <- grab_sims(.d = seq(-0.15, 0.05, 0.025),
#'                                           .eta = -1:1 * etas[[2]],
#'                                           .add_var = c(0.01, 0.05, 0.1),
#'                                           .sigma_N = c(0, 0.05, 0.1, 0.2, 0.3),
#'                                           .sigma_V = c(0, 0.05, 0.1, 0.2, 0.3),
#'                                           .vary_d2 = FALSE) %>%
#'     group_by(d, eta, add_var, sigma_N, sigma_V, vary_d2, rep) %>%
#'     summarize(n_spp = spp[N > 0] %>% unique() %>% length()) %>%
#'     ungroup() %>%
#'     mutate(eta = factor(eta, levels = -1:1 * etas[[2]],
#'                         labels = c("sub-additive", "additive", "super-additive")))
#'
#'
#' stoch_coexist_d1_p_fun <- function(.x) {
#'
#'     stoch_vary_d1_coexist_spp_df %>%
#'         filter(eta == .x) %>%
#'         filter(sigma_N %in% c(0, 0.1, 0.2),
#'                sigma_V %in% c(0, 0.1, 0.2),
#'                add_var == 0.05) %>%
#'         mutate(sigma_N = factor(sigma_N, levels = sort(unique(sigma_N)),
#'                                 labels = sprintf("sigma[N] == %.2f",
#'                                                  sort(unique(sigma_N)))),
#'                sigma_V = factor(sigma_V, levels = sort(unique(sigma_V)),
#'                                 labels = sprintf("sigma[V] == %.2f",
#'                                                  sort(unique(sigma_V)))),
#'                add_var = factor(add_var, levels = sort(unique(add_var))),
#'                extinct = factor(n_spp == 0)) %>%
#'         mutate(n_spp = n_spp / 100) %>%
#'         # ggplot(aes(d, n_spp, color = add_var)) +
#'         ggplot(aes(d, n_spp)) +
#'         ggtitle(.x) +
#'         geom_hline(yintercept = 0, color = "gray80") +
#'         geom_vline(xintercept = 0, color = "gray80") +
#'         geom_jitter(aes(shape = extinct, size = extinct), width = 0.002,
#'                     height = 0, color = "firebrick") +
#'         facet_grid(sigma_N ~ sigma_V, labeller = label_parsed) +
#'         # scale_color_viridis_d(expression(italic(sigma[i])^2),
#'         #                       begin = 0.1, end = 0.85) +
#'         scale_shape_manual(values = c(1, 4), guide = FALSE) +
#'         scale_size_manual(values = c(1, 3), guide = FALSE) +
#'         scale_y_continuous("Proportion of species that survive",
#'                            breaks = c(0, 0.4, 0.8), limits = c(-0.01, 1)) +
#'         scale_x_continuous(expression(italic(d[1]) ~
#'                                           ("with" ~ italic(d[2]) == 0.1)),
#'                            breaks = c(-0.1, 0)) +
#'         guides(color = guide_legend(override.aes = list(shape = 19, size = 2))) +
#'         theme(strip.text = element_text(size = 10),
#'               strip.text.y = element_text(angle = 0),
#'               legend.position = "top",
#'               plot.title = element_text(hjust = 0.5,
#'                                         margin = margin(0,0,0,b=12)))
#'
#' }
#'
#' stoch_coexist_d1_ps <- map(c("sub-additive", "additive", "super-additive"),
#'                            stoch_coexist_d1_p_fun)
#' names(stoch_coexist_d1_ps) <- c("sub", "add", "super")
#'
#' # stoch_coexist_d1_ps[["sub"]]
#' # stoch_coexist_d1_ps[["add"]]
#' # stoch_coexist_d1_ps[["super"]]
#'
#' stoch_coexist_ps <- map(stoch_coexist_ps,
#'                         ~.x + theme(axis.title.y = element_blank(),
#'                                     strip.text.y = element_blank()))
#' stoch_coexist_d1_ps <- map(stoch_coexist_d1_ps,
#'                            ~.x + theme(axis.title.y = element_blank(),
#'                                        axis.text.y = element_blank()))
#' for (f in c("sub", "add")) {
#'     stoch_coexist_ps[[f]] <- stoch_coexist_ps[[f]] +
#'         theme(axis.title.x = element_blank(),
#'               axis.text.x = element_blank())
#'     stoch_coexist_d1_ps[[f]] <- stoch_coexist_d1_ps[[f]] +
#'         theme(axis.title.x = element_blank(),
#'               axis.text.x = element_blank())
#' }
#' for (f in c("add", "super")) {
#'     stoch_coexist_ps[[f]] <- stoch_coexist_ps[[f]] +
#'         theme(strip.text.x = element_blank())
#'     stoch_coexist_d1_ps[[f]] <- stoch_coexist_d1_ps[[f]] +
#'         theme(strip.text.x = element_blank())
#' }
#'
#' stoch_coexist_p <- ggarrange(stoch_coexist_ps[["sub"]],
#'           stoch_coexist_d1_ps[["sub"]],
#'           stoch_coexist_ps[["add"]],
#'           stoch_coexist_d1_ps[["add"]],
#'           stoch_coexist_ps[["super"]],
#'           stoch_coexist_d1_ps[["super"]],
#'           ncol = 2, left = "Proportion of species that survive",
#'           draw = FALSE)
#'
#'
#' if (.RESAVE_PLOTS) {
#'     save_plot(stoch_coexist_p, 6.5, 7, .prefix = "S3-")
#' }
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' # =============================================================================*
#' # =============================================================================*
#'
#' # <not updated> Additive --> Sub-additive ?? ----
#'
#' # =============================================================================*
#' # =============================================================================*
#'
#'
#'
#' # What direction does selection move them?
#'
#' ss_t <- sauron:::sel_str_cpp
#' trnorm <- sauron:::trunc_rnorm_cpp
#'
#' D <- matrix(0, 2, 2)
#' diag(D) <- c(-0.1, 0.1)
#' C <- diag(2)
#'
#' V <- rbind(c(1.367653, 0.9889897),
#'            c(1.458830, 1.7362264))
#'
#' N <- c(9.6430886, 0.9798823, 0.3046547)
#'
#'
#' dir_df <- crossing(.v1 = seq(0, 4, 0.2),
#'                    .v2 = seq(0, 4, 0.2)) %>%
#'     mutate(across(.fns = round, digits = 3)) %>%
#'     pmap_dfr(function(.v1, .v2) {
#'         ss <- ss_t(cbind(V, c(.v1, .v2)),
#'                    N,
#'                    formals(quant_gen)$f,
#'                    formals(quant_gen)$a0,
#'                    C,
#'                    formals(quant_gen)$r0,
#'                    D)
#'         tibble(V1 = .v1, V2 = .v2, V1_delta = ss[1,3], V2_delta = ss[2,3])
#'     })
#'
#' dir_df %>%
#'     ggplot(aes(V1, V2)) +
#'     stable_points(0, return_geom = TRUE, color = "gray70") +
#'     geom_segment(aes(xend = V1 + 0.15 * V1_delta, yend = V2 + 0.15 * V2_delta),
#'                  arrow = arrow(length = unit(0.05, "inches"))) +
#'     scale_fill_viridis_c() +
#'     scale_x_continuous("Axis 1\n(conflicting)", breaks = c(0, 2, 4)) +
#'     scale_y_continuous("(ameliorative)\nAxis 2", breaks = c(0, 2, 4)) +
#'     coord_equal() +
#'     NULL
#'
#'
#' # set.seed(043506945)
#' maxt <- 2e4
#' .V_sigmas <- c(0.1, 0.2)
#'
#' V_test <- matrix(0, maxt+1, 4)
#' V_test[1,1:2] <- c(0.4, sqrt(4-0.4^2))
#' # V_test[1,1:2] <- c(sqrt(4-0.4^2), 0.4)
#'
#' for (k in 1:2) {
#'     if (.V_sigmas[k] > 0) {
#'         V_test[1,k+2] <- V_test[1,k] * rlnorm(1, - .V_sigmas[k]^2 / 2,
#'                                               .V_sigmas[k])
#'     } else V_test[1,k+2] <- V_test[1,k]
#' }
#'
#'
#' for (t in 1:maxt) {
#'
#'     deltaV <- 0.05 * ss_t(cbind(V_test[t,3:4]),
#'                           pop_sizes(1, 0, c(-0.1, 0.1)),
#'                           formals(quant_gen)$f,
#'                           formals(quant_gen)$a0,
#'                           C,
#'                           formals(quant_gen)$r0,
#'                           D)
#'
#'     V_test[t+1,1:2] <- pmax(0, V_test[t,1:2] + deltaV[,ncol(deltaV)])
#'
#'     for (k in 1:2) {
#'         if (.V_sigmas[k] > 0) {
#'             V_test[t+1,k+2] <- V_test[t+1,k] * rlnorm(1, - .V_sigmas[k]^2 / 2,
#'                                                       .V_sigmas[k])
#'         } else V_test[t+1,k+2] <- V_test[t+1,k]
#'     }
#'
#' }
#'
#' colnames(V_test) <- c("V1", "V2", "Vp1", "Vp2")
#'
#'
#' # Standard deviation of difference between genotype and phenotype:
#' sqrt(sum((V_test[,"Vp1"] - V_test[,"V1"])^2) / nrow(V_test))
#' #' Method 1:
#' #'   - 0.6783879 for .sigma_V = 0.4
#' #'   - 0.2892102 for .sigma_V = 0.2
#' #'   - 0.1341219 for .sigma_V = 0.1
#' #' Method 2:
#' #'   - 0.411026 for .sigma_V = 0.4
#' #'   - 0.1908039 for .sigma_V = 0.2
#' #'   - 0.09974028 for .sigma_V = 0.1
#' #' Method 3:
#' #'   - 0.6127865 for .sigma_V = 0.4
#' #'   - 0.2876067 for .sigma_V = 0.2
#' #'   - 0.1336868 for .sigma_V = 0.1
#'
#'
#'
#'
#' V_test %>%
#'     # .[1:20e3L,] %>%
#'     # .[seq(1, nrow(.), 100),] %>%
#'     as_tibble() %>%
#'     ggplot(aes(V1, V2)) +
#'     geom_path(aes(Vp1, Vp2), color = "dodgerblue", alpha = 0.2) +
#'     geom_path() +
#'     geom_point(data = V_test[1,,drop=FALSE] %>% as_tibble(),
#'                shape = 1, size = 4, color = "firebrick") +
#'     geom_point(data = V_test[nrow(V_test),,drop=FALSE] %>% as_tibble(),
#'                shape = 4, size = 4, color = "firebrick") +
#'     stable_points(0, return_geom = TRUE, color = "gray40") +
#'     geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray70") +
#'     coord_equal() #xlim = c(0, 3), ylim = c(0, 3))
#'
#'
#'
#' V_test %>%
#'     # .[1:100,] %>%
#'     as_tibble() %>%
#'     mutate(angle = atan(Vp1 / Vp2) * 180 / pi,
#'            radius = sqrt(Vp1^2 + Vp2^2),
#'            time = 1:n()) %>%
#'     gather(type, value, angle:radius) %>%
#'     ggplot(aes(time, value)) +
#'     geom_line(color = "dodgerblue", alpha = 0.5) +
#'     geom_line(data = V_test %>%
#'                   # .[1:100,] %>%
#'                   as_tibble() %>%
#'                   mutate(angle = atan(V1 / V2) * 180 / pi,
#'                          radius = sqrt(V1^2 + V2^2),
#'                          time = 1:n()) %>%
#'                   gather(type, value, angle:radius),
#'               color = "black") +
#'     geom_hline(data = tibble(type = c("angle", "radius"),
#'                              yint = c(45, 2)),
#'                aes(yintercept = yint), color = "firebrick") +
#'     geom_hline(data = tibble(yint = c(0, 0)),
#'                aes(yintercept = yint), color = "gray70") +
#'     facet_wrap(~ type, scales = "free", ncol = 1) +
#'     theme(strip.text = element_text(size = 16),
#'           panel.spacing = unit(2, "lines"))
#'
#'
#' #'
#' #' The simulations above show that when stochasticity for evolution is
#' #' a normal distribution instead of lognormal, it doesn't inextricably move
#' #' toward the center point on the ring (i.e., it doesn't become sub-additive).
#' #'
#'
#'
#'
#' V_test %>%
#'     as_tibble() %>%
#'     mutate(out = sqrt(Vp1^2 + Vp2^2) > 2,
#'            diff1 = V1 - lag(V1),
#'            diff2 = V2 - lag(V2)) %>%
#'     filter(!is.na(diff1)) %>%
#'     gather("axis", "value", diff1, diff2) %>%
#'     ggplot(aes(value, ..density..)) +
#'     geom_vline(xintercept = 0, linetype = 2, color = "gray70") +
#'     geom_freqpoly(aes(color = out), bins = 50) +
#'     facet_grid(~ axis, scales = "free") +
#'     scale_color_manual(values = c("firebrick", "dodgerblue"))
#'
#'
#'
#'
#' .V0 = "unrestricted"
#' .eta_sign = 0
#' .d1 = -0.1
#' .sigma_V = 0.1
#' .sigma_N = 0
#'
#'
#'
#' which_switch <- 3 # which spp to switch for unrestricted
#'
#' .n <- 5
#' .q = 2
#' stopifnot(is.numeric(.d1) && length(.d1) == 1)
#' .ds <- c(.d1, abs(.d1))
#' stopifnot(is.numeric(.eta_sign) && length(.eta_sign) == 1 &&
#'               .eta_sign %in% -1:1)
#' .V0 <- match.arg(.V0, c("restricted", "unrestricted"))
#' .lab <- .V0
#'
#' .eta <- .eta_sign * etas[[2]]
#'
#' if (.lab == "restricted" && .eta < 0) {
#'     stop(paste("\n'restricted' with sub-additivity is not programmed bc",
#'                "there's only one stable axis state.",
#'                "Thus, there is no way to restrict starting axes to",
#'                "be outside of that state's basin of attraction."))
#' }
#'
#' if (.eta < 0) {
#'     .V0 <- rbind(seq(3, 2, length.out = 5),
#'                  seq(2, 3, length.out = 5))
#' } else if (.eta == 0) {
#'     .V0 <- rbind(seq(1.5, 0, length.out = 5),
#'                  seq(1.6, 3.1, length.out = 5))
#'     if (.lab == "unrestricted") {
#'         v2 <- .V0[2, which_switch]
#'         .V0[2, which_switch] <- .V0[1, which_switch]
#'         .V0[1, which_switch] <- v2
#'     }
#' } else {
#'     .V0 <- rbind(seq(1.2, 0, length.out = 5),
#'                  seq(1.3, 2.5, length.out = 5))
#'     if (.lab == "unrestricted") {
#'         v2 <- .V0[2, which_switch]
#'         .V0[2, which_switch] <- .V0[1, which_switch]
#'         .V0[1, which_switch] <- v2
#'     }
#' }
#' .V0 <- round(.V0, 3)
#'
#' if (.sigma_V == 0 && .sigma_N == 0) {
#'     .nreps <- 1
#'     .N_THREADS <- 1
#' } else .nreps <- 12
#'
#'
#'
#' set.seed(1415272918)
#' qg <- quant_gen(q = .q, eta = .eta, d = .ds,
#'                 n_reps = .nreps, n = ncol(.V0),
#'                 V0 = .V0,
#'                 sigma_V0 = 0,
#'                 sigma_V = .sigma_V,
#'                 sigma_N = .sigma_N,
#'                 spp_gap_t = 500L,
#'                 final_t = 20e3L,
#'                 save_every = 1L,
#'                 add_var = rep(0.05, .n),
#'                 n_threads = .N_THREADS,
#'                 show_progress = FALSE)
#' qg0 <- quant_gen(q = .q, eta = .eta, d = .ds,
#'                 n_reps = 1, n = ncol(.V0),
#'                 V0 = .V0,
#'                 sigma_V0 = 0,
#'                 sigma_V = 0,
#'                 sigma_N = 0,
#'                 spp_gap_t = 500L,
#'                 final_t = 20e3L,
#'                 save_every = 1L,
#'                 add_var = rep(0.05, .n),
#'                 n_threads = 1,
#'                 show_progress = FALSE)
#'
#'
#' nv <- qg %>%
#'     .[["nv"]] %>%
#'     mutate(axis = paste0("V", axis)) %>%
#'     nest(value = c(geno, pheno)) %>%
#'     spread(axis, value) %>%
#'     unnest(c(V1, V2), names_sep = "_") %>%
#'     rename(V1 = V1_geno,
#'            V2 = V2_geno,
#'            Vp1 = V1_pheno,
#'            Vp2 = V2_pheno)
#'
#' nv %>%
#'     mutate(id = interaction(spp, rep)) %>%
#'     arrange(time) %>%
#'     # ggplot(aes(V1, V2, color = spp)) +
#'     ggplot(aes(Vp1, Vp2, color = spp)) +
#'     geom_path(aes(group = id)) +
#'     geom_point(data = nv %>%
#'                    group_by(spp) %>%
#'                    summarize(V1 = V1[time == min(time)][1],
#'                              V2 = V2[time == min(time)][1])) +
#'     facet_wrap(~ rep, nrow = 3) +
#'     coord_equal(xlim = c(0, 3.11), ylim = c(0, 3.11)) +
#'     scale_color_brewer(palette = "Dark2", guide = FALSE) +
#'     ylab("(ameliorative)\nAxis 2") +
#'     xlab("Axis 1\n(conflicting)")
#'
#' nv %>%
#'     mutate(id = interaction(spp, rep)) %>%
#'     ggplot(aes(time / 1000, N, color = spp)) +
#'     geom_path(aes(group = id)) +
#'     facet_wrap(~ rep, nrow = 3) +
#'     scale_color_brewer(palette = "Dark2", guide = FALSE) +
#'     ylab("Abundance") +
#'     xlab("Time (× 1000)")
