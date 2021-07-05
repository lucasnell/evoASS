

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
              size = 10 / 2.83465, nudge_x = -0.05, nudge_y = 500, hjust = 1)
comms_d_p_list[[2]] <- comms_d_p_list[[2]] +
    geom_point(data = filter(focal_stable_comms, d1 == -0.6, spp == 1) %>%
                   add_col_eta_cols(),
                size = 2, shape = 1) +
    geom_text(data = filter(focal_stable_comms, d1 == -0.6, spp == 1) %>%
                  mutate(O = case_when(!comm %in% c("(iv)", "(v)") ~ O + 1100,
                                       TRUE ~ O - 3000)) %>%
                  add_col_eta_cols(),
              aes(label = comm),
              size = 10 / 2.83465, vjust = 0) +
    coord_cartesian(ylim = c(0, NA))



focal_comm_p_list <- stable_comm_fit %>%
    distinct(eta, d1, d2, state) %>%
    rename_with(~ paste0(".", .x)) %>%
    pmap(focal_comm_hmp_fun, legend.position = "none",
         axis.title = element_blank())
focal_comm_p_list[c(2:3,5:6)]  <- focal_comm_p_list[c(2:3,5:6)] %>%
    map(~ .x + theme(axis.text.y = element_blank()))
focal_comm_p_list[1:3] <- focal_comm_p_list[1:3] %>%
    map(~ .x + theme(axis.text.x = element_blank()))



focal_comm_legend <- get_legend(focal_comm_p_list[[1]] +
                                    theme(legend.position = "top"))

focal_comm_p <- arrangeGrob(focal_comm_legend,
                            textGrob("Ameliorative axis", rot = 90, vjust = 1),
                            gtable_rbind(focal_comm_p_list[1:3] %>%
                                             map(make_gf) %>%
                                             do.call(what = gtable_cbind),
                                         focal_comm_p_list[4:6] %>%
                                             map(make_gf) %>%
                                             do.call(what = gtable_cbind)),
                            make_gf(BLANK), textGrob("Conflicting axis"),
                            widths = c(0.05, 1), heights = c(0.1, 1, 0.1),
                            nrow = 3, ncol = 2,
                            layout_matrix = rbind(c(1,1), c(2,3), c(4,5)),
                            vp = viewport(height = unit(0.6, "npc"),
                                          y = 0, x = 0,
                                          just = c("left", "bottom")),
                            top = textGrob("(b)", x = 0, hjust = 0,
                                           gp = gpar(fontsize = 14)))






comm_p <- function() {
    grid.newpage()
    arrangeGrob(ggarrange(plots = comms_d_p_list, draw = FALSE,
                          left = textGrob("Scaled community size", y = 0.55,
                                          vjust = 0.7, rot = 90), clip = "off"),
                top = textGrob("(a)", x = 0, hjust = 0,
                               gp = gpar(fontsize = 14)),
                vp = viewport(height = unit(0.4, "npc"),
                              y = 1, x = 0,
                              just = c("left", "top"))) %>%
        grid.draw()
    grid.draw(focal_comm_p)
    invisible(NULL)
}


if (.RESAVE_PLOTS) save_plot(comm_p, 5, 7, .prefix = "3-")

# ... Fig S1: `r` for tradeoffs across axis space----
#'
#' This explains why Omega doesn't explain trends across changes in eta:
#' The costs to the growth rate associated with super-additive tradeoffs
#' are overall greater across this range.
#'
tradeoffs_r_p <- crossing(v1 = round(seq(0, 3, 0.01), 2), v2 = v1,
                          e = c(-1,0,1)*0.6) %>%
    mutate(M = { formals(quant_gen)$r0 - formals(quant_gen)$f *
            (v1^2 + 2 * e * v1 * v2 + v2^2) } %>%
        exp() %>%
        identity()) %>%
    mutate(e = factor(e, levels = sort(unique(e)),
                      labels = paste0(c("sub-", "", "super-"), "additive"))) %>%
    # group_by(e) %>%
    # summarize(M = mean(M))
    ggplot(aes(v1, v2, fill = M)) +
    geom_raster(interpolate = FALSE) +
    # geom_contour(aes(z = M)) +
    facet_wrap(~ e) +
    xlab("Neutral axis 1") +
    ylab("Neutral axis 2") +
    scale_fill_viridis_c(expression(e^{r[0] - f ~ bold(v)[i]^T *
                                        bold(C) * bold(v)[i]})) +
    coord_equal()


if (.RESAVE_PLOTS) save_plot(tradeoffs_r_p, 6, 2.5, .prefix = "S1-")




# =============================================================================*
# =============================================================================*

# Fig 4: Stabilizers ----

# =============================================================================*
# =============================================================================*


#' Also maybe show plot here about how it's not always most conducive to
#' coexistence when ameliorative axis is really strong.
#' (Because when it's really strong, all species investing in it isn't
#'  a stable community.)



stab_sims <- function(.dd, .n = 10) {

    n_lower <- .dd$n_lower
    d1 <- .dd$d1
    d2 <- .dd$d2

    # n_lower = 9; d1 = -1; d2 = 0; .n = 10
    # rm(.dd, n_lower, d1, d2, .n, .V0_0, sim0, L, X)

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

    out <- mutate(.dd, V1 = NA_real_, V2 = NA_real_,
                  O = NA_real_, O_c = NA_real_, O_a = NA_real_,
                  V = list(NA),
                  N = list(NA))
    if (any(is.na(sim0$nv$geno)) || length(unique(sim0$nv$spp)) < .n) return(out)

    # stable?
    L <- sim0 %>%
        jacobians() %>%
        .[[1]] %>%
        eigen(only.values = TRUE) %>%
        .[["values"]] %>%
        Re() %>%
        max()
    if (L >= 1) return(out)

    X <- sim0 %>%
        .[["nv"]] %>%
        pivot()

    # Vector of `N_j * exp(v_i^T * D * v_i)`
    O_vec <- map_dbl(1:nrow(X),
                 ~ X$N[.x] * exp(- d1 * X$V1[.x]^2 - d2 * X$V2[.x]^2))

    if (n_lower == 0) O_conf <- NA_real_ else O_conf <- sum(O_vec[-1])
    if (n_lower == .n) O_amel <- NA_real_ else O_amel <- sum(O_vec[-.n])

    .V <- X %>%
        select(V1, V2) %>%
        map_dbl(max)

    mutate(out,
           V1 = .V[[1]], V2 = .V[[2]],
           O = sum(O_vec), O_c = O_conf, O_a = O_amel,
           V = list(t(X[,c("V1", "V2")])),
           N = list(X$N))
}

# Takes ~ 16 sec
stab_sim_df <- crossing(n_lower = 0:formals(stab_sims)$.n,
                        d1 = - c(0, 0.1, 0.5, 1),
                        d2 = abs(d1)) %>%
    split(1:nrow(.)) %>%
    mclapply(stab_sims) %>%
    do.call(what = bind_rows)


# Should all be zero:
stab_sim_df %>% filter(is.na(V1)) %>% nrow()
stab_sim_df %>% filter(is.na(V2)) %>% nrow()
stab_sim_df %>% filter(is.na(O)) %>% nrow()
stab_sim_df %>% filter(is.na(O_c), n_lower > 0) %>% nrow()
stab_sim_df %>% filter(is.na(O_a), n_lower < 10) %>% nrow()



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

# ... Fig S2: effects of ameliorative axis----

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


if (.RESAVE_PLOTS) save_plot(stab_p2, 5, 5, .prefix = "S2-")




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
        # geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
        # geom_vline(xintercept = 0, linetype = 1, color = "gray70") +
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
                           .sigma_V = list(list(c(0.1, 0.1)),
                                           list(c(0.1, 0.20)),
                                           list(c(0.20, 0.1))),
                           .d_mults = 1) %>%
    pmap_dfr(one_stoch_sim) %>%
    select(-eta) %>%
    mutate(sigma_v1 = map_dbl(sigma_V, ~ .x[1]),
           sigma_v2 = map_dbl(sigma_V, ~ .x[2]),
           sigma_V = paste(sigma_v1, sigma_v2, sep = "__") %>%
               factor(levels = c("0.1__0.1", "0.1__0.2", "0.2__0.1")))


stoch_sub_V_p_list <- map(sort(unique(stoch_sub_sims$sigma_V)),
                          stoch_comm_V_plotter, .sims = stoch_sub_sims)
stoch_sub_O_p_list <- map(sort(unique(stoch_sub_sims$sigma_V)),
                          stoch_comm_O_plotter, .sims = stoch_sub_sims)

stoch_sub_comm_p <- stoch_comm_p_combiner(stoch_sub_V_p_list,
                                          stoch_sub_O_p_list,
                                          rbind(c(0.10, 0.10),
                                                c(0.10, 0.20),
                                                c(0.20, 0.10)))

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

# stoch_p()

if (.RESAVE_PLOTS) save_plot(stoch_p, 6, 6, .prefix = "5-")







# ----------------------------------------------------------------------------*
# ----------------------------------------------------------------------------*
# ... Fig S3: alt. scenario where sigma_V increases Omega ----


# Takes ~7 sec
set.seed(1412303795)
alt_stoch_super_sims <- tibble(.eta = 1e-2,
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


if (.RESAVE_PLOTS) save_plot(alt_stoch_p, 6, 3, .prefix = "S3-")



