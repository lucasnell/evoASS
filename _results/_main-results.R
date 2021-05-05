
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
# if (file.exists(".Rprofile")) source(".Rprofile")


# Clean captions
cc <- function(.x) {
    .x <- gsub("\n", "", .x)
    .x <- gsub("\\s+", " ", .x)
    return(.x)
}

# whether to re-do simulations (use rds files otherwise)
.REDO_SIMS <- FALSE
# whether to re-save plots
.RESAVE_PLOTS <- FALSE
# number of threads to use for simulations
.N_THREADS <- max(1, parallel::detectCores()-2)



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

# Arguments for `egg::ggarrange`
label_args <- list(gp = grid::gpar(font = 1,
                                   fontsize = 16),
                   x = unit(0, "npc"),
                   y = unit(1, "npc"),
                   hjust = 0, vjust = 1)


pivot <- function(.df) {
    mutate(.df, axis = paste0("V", axis)) %>%
        pivot_wider(names_from = axis, values_from = geno)
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



.d <- function(.strong, .barely) {
    .strong <- match.arg(.strong, c("conflicting", "ameliorative"))
    if (.strong == "conflicting") {
        z <- c(-2.5, 0.6)
        if (.barely) z[1] <- z[1] / 2
    } else {
        z <- c(-0.6, 4)
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


    # .strong = "ameliorative"
    # .barely <- TRUE
    # .n = 2
    # .eta = 0.6
    # .state = "both"

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
    arrange(eta, strong, barely, z1) %>%
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
#' lambda = 1.069842, strong = conflicting, barely = FALSE, n = 2,
#' eta = 0.0, state = below
#'     - small perturbation and it goes to 1 species
#' lambda = 1.067329, strong = conflicting, barely = FALSE, n = 2,
#' eta = 0.6, state = below
#'     - small perturbation and it goes to 1 species
#' lambda = 1.001227, strong = ameliorative, barely = FALSE, n = 2,
#' eta = -0.6, state = above
#'     - small perturbation and it goes to one species investing in both,
#'       the other not investing
#' lambda = 1.004468, strong = ameliorative, barely = FALSE, n = 2,
#' eta = 0.0, state = above
#'     - small perturbation and it goes to one investing in both,
#'       another investing in neither
#' lambda = 1.004679, strong = ameliorative, barely = FALSE, n = 2,
#' eta = 0.6, state = above
#'     - small perturbation and it goes to one species investing in axis 2,
#'       the other not investing
#' lambda = 1.002689, strong = ameliorative, barely = TRUE, n = 2,
#' eta = 0.0, state = above
#'     - small perturbation and it goes to one species investing 0.361 in V1
#'       and 1.20 in V2, the other investing ~0 in both
#' lambda = 1.003003, strong = ameliorative, barely = TRUE, n = 2,
#' eta = 0.6, state = above
#'     - small perturbation and it goes to one species investing 0 in V1
#'       and 1.26 in V2, the other investing ~0 in both
#' lambda = 1.030829, strong = conflicting, barely = TRUE, n = 2,
#' eta = 0.6, state = below
#'     - small perturbation and it goes to 1 species




# Three species communities:
# --------------*
# This community isn't found in `hist_d_sims` (it's pretty hard to achieve)
Z <- quant_gen(q = 2, eta = -0.6,
               V0 = matrix(c(2, 2.1), 2, 3),
               N0 = rep(1, 3),
               n = 3, d = .d("conflicting", FALSE), n_reps = 1,
               spp_gap_t = 0L,
               final_t = 50e3L,
               save_every = 0L,
               sigma_V0 = 0,
               show_progress = FALSE) %>%
    # jacobians() %>%
    # .[[1]] %>%
    # eigen(only.values = TRUE) %>%
    # .[["values"]] %>%
    # Re() %>%
    # max()
    .[["nv"]] %>%
    pivot()



#' All these are stable, except for...
#'
#' lambda = 1.0044020, strong = ameliorative, barely = FALSE, n = 3, eta = 0.6,
#' all spp investing in V2, but at different magnitudes
#' (2 at 0.5253733, 1 at 0.950722)
#'     - small perturbation of either of two species at 0.5253733, and
#'       it goes to 2 species at V2 = 1.27, one species at V2 ≈ 0
#'
#' lambda = 1.001222, strong = ameliorative, barely = TRUE, n = 3, eta = 0.6,
#' all spp investing in V2, but at different magnitudes
#' (2 at 0.7520488, 1 at 0.8920994)
#'     - small perturbation of either of two species at 0.7520488, and
#'       it goes to 2 species at V2 = 1.26, one species at V2 ≈ 0
#'
#'



comms$three <- hist_d_sims %>%
    filter(time == max(time) & eta != 0, n_spp > 1) %>%
    group_by(eta, strong, barely, state, start) %>%
    mutate(spp3 = all(1:3 %in% spp)) %>%
    ungroup() %>%
    filter(spp3) %>%
    select(-n_spp, -spp3, -time) %>%
    arrange(eta, strong, barely, state, start, spp) %>%
    group_by(eta, strong, barely, state, start) %>%
    summarize(N = list(N), V = list(rbind(V1, V2)),
              .groups = "drop") %>%
    add_row(eta = -0.6, strong = "conflicting", barely = FALSE,
            N = list(Z$N), V = list(t(Z[,c("V1", "V2")]))) %>%
    mutate(z1 = map_dbl(V, ~ sum(.x[1,])),
           z2 = map_dbl(V, ~ sum(.x[2,]))) %>%
    arrange(eta, strong, barely, z1) %>%
    group_by(eta, strong, barely) %>%
    mutate(dup = dup_comm(V, .prec = 1e-3)) %>%
    filter(!dup) %>%
    mutate(comm = 1:n()) %>%
    ungroup() %>%
    select(-state, -start, -z1, -z2, -dup) %>%
    # Remove unstable communities:
    filter(barely |
               !map_lgl(V, ~ all(.x[1,] == 0) &
                            sum((.x[2,] - 0.5253733)^2 < 1e-10) == 2 &
                            sum((.x[2,] - 0.9507220)^2 < 1e-10) == 1)) %>%
    filter(!barely |
               !map_lgl(V, ~ all(.x[1,] == 0) &
                            sum((.x[2,] - 0.7520488)^2 < 1e-10) == 2 &
                            sum((.x[2,] - 0.8920994)^2 < 1e-10) == 1)) %>%
    arrange(eta, strong, barely, comm) %>%
    # split(1:nrow(.)) %>%
    # map_dfr(check_stable) %>%
    # arrange(desc(lambda))
    mutate(V1 = map(V, ~ .x[1,]),
           V2 = map(V, ~ .x[2,])) %>%
    select(-V) %>%
    unnest(c(N, V1, V2)) %>%
    mutate(spp = factor(rep(1:3, n() / 3), levels = 1:3)) %>%
    mutate(comm = factor(comm, levels = sort(unique(comm)))) %>%
    select(eta, strong, barely, comm, spp, everything())

# rm(Z)


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
    # Takes ~2.5 min
    hist_d_fit$three <- comms$three %>%
        split(interaction(.$eta, .$strong, .$barely, .$comm, drop = TRUE)) %>%
        mclapply(comm_fit, mc.cores = .N_THREADS) %>%
        do.call(what = rbind)
    saveRDS(hist_d_fit, rds("hist_d_fit"))
} else {
    hist_d_fit <- readRDS(rds("hist_d_fit"))
}


hist_d_fit$two <- hist_d_fit$two %>%
    mutate(outcome = case_when(!surv ~ "invader excluded",
                               coexist_spp == 3 ~ "3 species coexist",
                               coexist_spp == 2 ~ "1 resident excluded",
                               coexist_spp == 1 ~ "2 residents excluded") %>%
               factor(levels = c("invader excluded", "2 residents excluded",
                                 "1 resident excluded", "3 species coexist")))
hist_d_fit$three <- hist_d_fit$three %>%
    mutate(outcome = case_when(!surv ~ "invader excluded",
                               coexist_spp == 4 ~ "4 species coexist",
                               coexist_spp == 3 ~ "1 resident excluded",
                               coexist_spp == 2 ~ "2 residents excluded",
                               coexist_spp == 1 ~ "3 residents excluded") %>%
    factor(levels = c("invader excluded", "3 residents excluded",
                      "2 residents excluded", "1 resident excluded",
                      "4 species coexist")))

stopifnot(!any(is.na(hist_d_fit$two$outcome)))
stopifnot(!any(is.na(hist_d_fit$three$outcome)))



unq_comm_p_fun <- function(.n, .s, .b, .e, ...) {

    # .n <- 2
    # .s <- "conflicting"
    # .b <- FALSE
    # .e <- -0.6
    # rm(.n, .s, .e, .xlab, .ylab, .title, .p, .pal)

    .nn <- c("one", "two", "three")[.n]

    .xlab <- ifelse(.s == "ameliorative", "Weak conflicting axis",
                    expression(bold("Strong") ~ "conflicting axis"))
    .ylab <- ifelse(.s == "conflicting", "Weak ameliorative axis",
                    expression(bold("Strong") ~ "ameliorative axis"))
    .title <- paste0(c("Sub-a", "A", "Super-a"),
                     "dditive")[which(-1:1 * 0.6 == .e)]
    .pal <- magma(n = hist_d_fit[[.nn]] %>% .[["outcome"]] %>%
                        levels() %>% length() %>% `-`(1),
                    begin = 0.2, end = 0.8) %>%
        rev()

    .p <- comms %>%
        .[[.nn]] %>%
        filter(strong == .s, barely == .b, eta == .e) %>%
        group_by(comm) %>%
        mutate(grp = group_spp(V1, V2, .prec = 0.001)) %>%
        group_by(comm, grp) %>%
        summarize(V1 = mean(V1), V2 = mean(V2), N = n(), .groups = "drop") %>%
        ggplot(aes(V1, V2)) +
        ggtitle(.title) +
        geom_raster(data = hist_d_fit[[.nn]] %>%
                        filter(strong == .s, barely == .b, eta == .e),
                    aes(fill = outcome), interpolate = FALSE) +
        geom_abline(slope = 1, intercept = 0, linetype = 2,
                    color  = "gray70") +
        geom_point(size = 7, shape = 21, color = "black", fill = "white") +
        geom_text(aes(label = N), size = 12 / 2.83465,
                  fontface = "bold", color = "black") +
        facet_wrap(~ comm, nrow = 1) +
        coord_equal(xlim = c(-0.2, 3), ylim = c(-0.2, 3)) +
        scale_fill_manual(NULL, values = c("white", .pal), drop = FALSE,
                          aesthetics = c("color", "fill")) +
        xlab(.xlab) +
        ylab(.ylab) +
        theme(strip.text = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 12,
                                        margin = margin(0,0,0,b=3)),
              legend.key = element_rect(colour = "black")) +
        theme(...) +
        NULL

    tibble(n_spp = .n, strong = .s, barely = .b, eta = .e, plot = list(.p))
}


comms_p_df <- crossing(.n = 2,
                       .b = FALSE,
                       .s = c("conflicting", "ameliorative"),
                       .e = c(-0.6, 0.6)) %>%
    pmap_dfr(unq_comm_p_fun)

comms_p_legend <- get_legend(comms_p_df[["plot"]][[1]])


comms_p_list <- comms_p_df %>%
    mutate(plot = ifelse(eta > 0,
                         map(plot, ~ .x + theme(axis.title.y = element_blank(),
                                                axis.text.y = element_blank())),
                         plot),
           plot = map(plot, ~ .x + theme(legend.position = "none",
                                         axis.title.x = element_blank())))


comms_p <- plot_grid(plot_grid(plot_grid(plotlist = comms_p_list[["plot"]][1:2],
                               nrow = 1, rel_widths = c(1, 1.9),
                               axis = "bl", align = "vh"),
                     textGrob("Weak conflicting axis", y = 1, vjust = 1),
                     plot_grid(plotlist = comms_p_list[["plot"]][3:4],
                               nrow = 1, rel_widths = c(1, 1.9),
                               axis = "bl", align = "vh"),
                     textGrob(expression(bold("Strong") ~
                                             "conflicting axis"),
                              y = 1, vjust = 1),
                     ncol = 1, rel_heights = c(1, 0.08, 1, 0.08)),
                     comms_p_legend,
                     nrow = 1, rel_widths = c(1, 0.3))


if (.RESAVE_PLOTS) save_plot(comms_p, 6.5, 5, .prefix = "3-")



comms_p_df3 <- crossing(.n = 3,
                        .b = FALSE,
                       .s = c("conflicting", "ameliorative"),
                       .e = c(-0.6, 0.6)) %>%
    pmap_dfr(unq_comm_p_fun)

comms_p_legend3 <- get_legend(comms_p_df3[["plot"]][[1]])


comms_p_list3 <- comms_p_df3 %>%
    mutate(plot = map(plot, ~ .x + theme(legend.position = "none")))

comms_p3 <- plot_grid(plot_grid(plotlist = comms_p_list3[["plot"]],
                                align = "hv", axis = "tl", ncol = 1, hjust = 0),
                      comms_p_legend3, nrow = 1, rel_widths = c(1, 0.3))



if (.RESAVE_PLOTS) save_plot(comms_p3, 7, 9, .prefix = "S1-")



# LEFT OFF ----
# Add plots (probably supplementary) for when barely == TRUE







# =============================================================================*
# =============================================================================*

# Fig 4: Stabilizers ----

# =============================================================================*
# =============================================================================*


z <- comms$three %>%
    filter(strong == "ameliorative", eta > 0) %>%
    group_by(comm) %>%
    mutate(z = sum(V2 > 1e-9)) %>%
    ungroup() %>%
    filter(z == 2) %>%
    select(-z)

# z[3,"V2"] <- z[3,"V1"]
# z[3,"V1"] <- 0

zz <- quant_gen(q = 2, eta = z$eta[1],
          V0 = cbind(rbind(z$V1, z$V2), c(0, 0.1)),
          N0 = c(z$N, 1),
          n = nrow(z)+1, d = .d(z$strong[1]), n_reps = 1,
          spp_gap_t = 0L,
          final_t = 50e3L,
          save_every = 10L,
          sigma_V0 = 0,
          show_progress = FALSE) %>%
    .[["nv"]] %>%
    pivot()

zz %>% filter(time == max(time))







zz %>%
    filter(time < 2000) %>%
    ggplot(aes(time, N, color = spp)) +
    geom_line()









# =============================================================================*
# =============================================================================*

# OBSOLETE Fig 3: "Global" coexistence ----

# =============================================================================*
# =============================================================================*


# --------------*
# __3A - d --> # spp ----
# --------------*

one_sim_combo <- function(.d, .eta, .add_var, .sigma_N, .sigma_V, .vary_d2) {

    # .d = 0.5; .eta = 0.5; .add_var = 0.05; .sigma_N = 0.5; .sigma_V = 0; .vary_d2 = TRUE

    .n <- 100

    if (.vary_d2) {
        .ds <- rep(.d, 2)
    } else .ds <- c(.d, 0.1)

    .seed <- sample.int(2^31 - 1, 1)
    set.seed(.seed)

    Z <- quant_gen(q = 2, eta = .eta, d = .ds, n = .n,
                   add_var = rep(.add_var, .n),
                   # spp_gap_t = 500L,
                   spp_gap_t = 0L,
                   final_t = 50e3L, n_reps = 12,
                   sigma_N = .sigma_N, sigma_V = .sigma_V,
                   save_every = 0L,
                   n_threads = .N_THREADS, show_progress = FALSE)

    if (.sigma_N == 0 && .sigma_V == 0) {
        Z$call[["eta"]] <- eval(.eta)
        Z$call[["d"]] <- eval(.ds)
        Z$call[["n"]] <- eval(.n)
        Z$call[["add_var"]] <- eval(rep(.add_var, .n))
        jacs <- jacobians(Z)
    } else jacs <- NULL


    if ("pheno" %in% colnames(Z$nv)) {
        .NV <- Z$nv %>%
            mutate(axis = paste0("V", axis)) %>%
            nest(value = c(geno, pheno)) %>%
            spread(axis, value) %>%
            unnest(c(V1, V2), names_sep = "_") %>%
            rename(V1 = V1_geno,
                   V2 = V2_geno,
                   Vp1 = V1_pheno,
                   Vp2 = V2_pheno)
    } else {
        .NV <- Z$nv %>%
            pivot()
        .NV$Vp1 <- NA_real_
        .NV$Vp2 <- NA_real_
    }

    return(list(NV = .NV %>%
                    mutate(d = .d, eta = .eta, add_var = .add_var,
                           sigma_N = .sigma_N, sigma_V = .sigma_V,
                           vary_d2 = .vary_d2),
                J = jacs,
                seed = .seed))
}

# Takes ~4 min w/ 6 threads and `.d` of length 7
coexist_d_spp_sims <- crossing(.d = -5:0,#seq(-0.15, 0, 0.025),
                            .eta = -1:1 * etas[[2]],
                            .add_var = 0.05,
                            .sigma_N = 0,
                            .sigma_V = 0,
                            .vary_d2 = TRUE) %>%
    pmap(one_sim_combo)


coexist_d_spp_df <- map_dfr(coexist_d_spp_sims, ~ .x[["NV"]]) %>%
    group_by(d, eta, rep) %>%
    summarize(n_spp = spp[N > 0] %>% unique() %>% length()) %>%
    ungroup() %>%
    mutate(eta = factor(eta, levels = -1:1 * etas[[2]],
                        labels = paste0(c("sub-", "", "super-"), "additive")))


# coexist_d_spp_p <-
coexist_d_spp_df %>%
    mutate(n_spp = n_spp / 100) %>%
    ggplot(aes(d, n_spp)) +
    geom_vline(xintercept = 0, color = "gray80", linetype = 1) +
    geom_hline(yintercept = 0, color = "gray80", linetype = 1) +
    geom_jitter(width = 0.002, height = 0, shape = 1, size = 1) +
    facet_wrap(~ eta, ncol = 1) +
    scale_y_continuous("Proportion of species that coexist",
                       breaks = c(0, 0.4, 0.8), limits = c(0, 1)) +
    # scale_x_continuous("Strength of\nconflicting axis",
    #                    breaks = c(-0.1, 0)) +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          plot.margin = margin(0,l=8,r=8,t=8))









# --------------*
# __3B - sigma_V/N --> # spp ----
# --------------*

# Takes ~3.5 min w/ 6 threads
coexist_stoch_spp_sims <- crossing(.d = 0,
                                  .eta = -1:1 * etas[[2]],
                                  .add_var = 0.05,
                                  .sigma_N = c(0, 0.05),
                                  .sigma_V = c(0, 0.05, 0.1),
                                  .vary_d2 = FALSE) %>%
    pmap(one_sim_combo)



coexist_stoch_spp_df <- coexist_stoch_spp_sims %>%
    map_dfr(~ .x[["NV"]]) %>%
    group_by(eta, sigma_N, sigma_V, rep) %>%
    summarize(n_spp = spp[N > 0] %>% unique() %>% length()) %>%
    ungroup() %>%
    mutate(eta = factor(eta, levels = -1:1 * etas[[2]],
                        labels = c("sub-additive", "additive",
                                   "super-additive"))) %>%
    mutate(p_spp = n_spp / 100,
           sigma_N = factor(sigma_N, levels = sort(unique(sigma_N))))


# coexist_stoch_spp_p <-
coexist_stoch_spp_df %>%
    ggplot(aes(sigma_V, p_spp)) +
    geom_vline(xintercept = 0, color = "gray80", linetype = 1) +
    geom_hline(yintercept = 0, color = "gray80", linetype = 1) +
    geom_jitter(aes(color = sigma_N),
                width = 0.002, height = 0,
                size = 1, shape = 1) +
    geom_line(data = coexist_stoch_spp_df %>%
                  group_by(eta, sigma_N, sigma_V) %>%
                  summarize(p_spp = mean(p_spp), .groups = "drop"),
              aes(color = sigma_N)) +
    geom_text(data = coexist_stoch_spp_df %>%
                  distinct(eta, sigma_N) %>%
                  filter(eta == "additive") %>%
                  mutate(sigma_V = c(0.1, 0.1), p_spp = c(0.75, 0.01)) %>%
                  mutate(lab = paste("sigma[N] ==", sigma_N)),
              aes(color = sigma_N, label = lab),
              size = 8 / 2.83465, parse = TRUE, hjust = 1, vjust = 0) +
    facet_wrap(~ eta, ncol = 1) +
    scale_y_continuous("Proportion of species that coexist",
                       breaks = c(0, 0.4, 0.8), limits = c(0, 1)) +
    scale_x_continuous("Axis evolution SD",
                       breaks = c(0, 0.05, 0.1)) +
    scale_color_viridis_d("Population SD",
                          option = "A", end = 0.7) +
    theme(legend.position = "none",
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          plot.margin = margin(0,l=8,r=8,t=8)) +
    NULL

with(formals(quant_gen), r0 / a0)






coexist_spp_p <- ggarrange(coexist_d_spp_p,
                           coexist_stoch_spp_p +
                               theme(axis.title.y = element_blank(),
                                     axis.text.y = element_blank()),
                           nrow = 1, draw = FALSE,
                           labels = c("a", "b"),
                           label.args = label_args)


# coexist_spp_p







if (.RESAVE_PLOTS) {
    save_plot(coexist_spp_p, 3, 4, .prefix = "3-")
}


# # (No effect at all, so no point in creating this.)
# inv_sims_one_p_fun("extinct",
#                    .fill_low = "#d01c8b",
#                    .fill_high = "#4dac26",
#                    .fill_limits = NULL)




# --------------*
# Fig S1 : sigma_i --> # spp ----
# --------------*

add_var_effect_df <- grab_sims(.d = 0,
                               .eta = -1:1 * etas[[2]],
                               .add_var = c(0.01, 0.05, 0.1),
                               .sigma_N = 0,
                               .sigma_V = 0,
                               .vary_d2 = FALSE) %>%
    group_by(add_var, eta, rep) %>%
    summarize(n_spp = spp[N > 0] %>% unique() %>% length()) %>%
    ungroup() %>%
    mutate(eta = factor(eta, levels = -1:1 * etas[[2]],
                        labels = paste0(c("sub-", "", "super-"), "additive")))

add_var_effect_p <- add_var_effect_df %>%
    mutate(n_spp = n_spp / 100) %>%
    ggplot(aes(add_var, n_spp)) +
    geom_hline(yintercept = 0, color = "gray80", linetype = 1) +
    geom_jitter(width = 0.002, height = 0, shape = 1, size = 1) +
    facet_wrap(~ eta, ncol = 1) +
    scale_color_viridis_d(expression(italic(sigma[i])^2),
                          begin = 0.1, end = 0.85, guide = FALSE) +
    scale_y_continuous("Proportion of species that coexist",
                       breaks = c(0, 0.4, 0.8), limits = c(0, 1)) +
    scale_x_continuous("Additive genetic variance")


if (.RESAVE_PLOTS) {
    save_plot(add_var_effect_p, 2.5, 4, .prefix = "S1-")
}








# --------------*
# Fig S2 - sigma_V - invader coexist ----
# --------------*


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
                        labels = paste0(c("'sub-", "'", "'super-"),
                                        "additive'"))) %>%
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
                        labels = paste0(c("'sub-", "'", "'super-"),
                                        "additive'"))) %>%
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
    # rm(fill_scale, make_eta_fct, midi, .sigma, ss_df)

    .par <- match.arg(.par, c("coexist", "replace", "reject", "extinct"))

    .which_sigma <- match.arg(.which_sigma, c("V", "N"))


    .fill_lab <- case_when(grepl("st$", .par) ~ paste0(.par, "ence"),
                           grepl("ct$", .par) ~ paste0(.par, "ion"),
                           grepl("ce$", .par) ~ paste0(.par, "ment"),
                           TRUE ~ "")
    .fill_lab <- bquote("Effect of" ~ sigma[.(.which_sigma)] ~
                            .(paste0("on ", .fill_lab, ":")))

    fill_scale <- scale_fill_gradient2(.fill_lab,
                                       low = .fill_low,
                                       mid = "white",
                                       high = .fill_high,
                                       midpoint = 0,
                                       limits = .fill_limits,
                                       labels = function(x) {
                                           z <- case_when(
                                               x > 0 ~ sprintf("{}+%.1f",x),
                                               x == 0 ~ "{}~0.0",
                                               TRUE ~ sprintf("%.1f",x))
                                           parse(text = z)
                                       })

    make_eta_fct <- function(..x) {
        factor(..x, levels = c(-0.6, 0, 0.6),
               labels = paste0(c("'sub-", "'", "'super-"),
                               "additive'"))
    }
    midi <- ceiling(formals(stable_points)[["line_n"]] / 2)

    .sigma <- sprintf("sigma_%s", .which_sigma)

    ss_df <- map_dfr(c(-0.6, 0, 0.6), ~ mutate(stable_points(.x),
                                               eta = make_eta_fct(.x)))

    # p <-
    stoch_inv_sims %>%
        filter(!!as.name(.sigma) == sprintf("sigma[%s] == %.4f",
                                            .which_sigma, 0.1)) %>%
        ggplot(aes(V1, V2)) +
        geom_raster(aes_string(fill = paste0(.par, "_diff"))) +
        geom_tile(data = determ_inv_sims %>% filter(!!as.name(.par) > 0),
                  fill = NA, color = "black", size = 0.1) +
        facet_grid(d1 ~ eta, label = label_parsed) +
        geom_point(data = ss_df %>%
                       filter(eta == "'super-additive'", V1 == 0),
                   size = 2, shape = 16, color = "black") +
        geom_point(data = ss_df %>%
                       filter(eta == "'super-additive'", V2 == 0),
                   size = 2, shape = 21, color = "black", fill = "white") +
        geom_point(data = ss_df %>%
                       filter(eta == "'sub-additive'"),
                   size = 2, shape = 16, color = "black") +
        geom_path(data = ss_df %>%
                      filter(eta == "'additive'"),
                  color = "gray50", size = 1) +
        geom_point(data = tibble(V1 = sqrt(2), V2 = sqrt(2),
                                 eta = make_eta_fct(0)),
                   shape = 16, size = 2, color = "black") +
        scale_x_continuous("Conflicting axis", breaks = c(0, 2, 4)) +
        scale_y_continuous("Ameliorative axis", breaks = c(0, 2, 4)) +
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
              legend.text.align = 0.5,
              plot.margin = margin(0,0,0,0),
              legend.position = "bottom") +
        guides(fill = guide_colorbar(title.position = "top")) +
        fill_scale

}

inv_sims_V_coexist_p <- inv_sims_one_p_fun("coexist")




if (.RESAVE_PLOTS) {
    save_plot(inv_sims_V_coexist_p,
              4, 4, .name = "S2-sigmaV_coexist")
}


# --------------*
# Fig S3 - sigma_V - invader replace ----
# --------------*

if (.RESAVE_PLOTS) {
    save_plot(inv_sims_one_p_fun("replace",
                                 .fill_low = "#e66101",
                                 .fill_high = "#5e3c99"),
              4, 4, .name = "S3-sigmaV_replace")
}


# --------------*
# Fig S4 - sigma_N - invader coexist ----
# --------------*


if (.RESAVE_PLOTS) {
    save_plot(inv_sims_one_p_fun("coexist", .which_sigma = "N"),
              4, 4, .name = "S4-sigmaN_coexist")
}









# =============================================================================*
# =============================================================================*

# Fig 4: Conditional coexistence ----

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
                   "there's only one stable axis state.",
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
                        pivot())
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
                        mutate(axis = paste0("V", axis)) %>%
                        nest(value = c(geno, pheno)) %>%
                        spread(axis, value) %>%
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
        # `axis_space` is a combination of starting axis values and eta:
        mutate(axis_space =
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


if (any(is.na(cond_coexist_df$axis_space))) {
    stop("\nERROR: unknown combination of V0 and eta in `cond_coexist_df`")
}
if (any(is.na(cond_coexist_stoch_df$axis_space))) {
    stop("\nERROR: unknown combination of V0 and eta in `cond_coexist_stoch_df`")
}



# Starting conditions and trajectories:


cond_coexist_sc_p_fun <- function(.d1) {
    .dd <- cond_coexist_df %>%
        filter(d1 == .d1) %>%
        split(interaction(.$axis_space, .$spp, drop = TRUE)) %>%
        map_dfr(~ mutate(.x,
                         first = time == min(time),
                         last = time == max(time))) %>%
        select(axis_space, time, spp, V1, V2, first, last) %>%
        arrange(axis_space, time)
    .dd %>%
        ggplot(aes(V1, V2)) +
        geom_abline(data = tibble(axis_space = factor(c("ii", "iv")),
                                  slp = 1, int = 0),
                    aes(slope = slp, intercept = int), linetype = 3, color = "gray70") +
        geom_line(data = bind_rows(stable_points(0), stable_points(0)) %>%
                      mutate(axis_space = factor(rep(c("iii", "v"),
                                                      each = stable_points %>%
                                                          formals() %>%
                                                          .[["line_n"]]))),
                  linetype = 2, color = "black") +
        geom_point(data = .dd %>% filter(first),
                   aes(color = spp), size = 1.5, na.rm = TRUE) +
        geom_point(data = map_dfr(etas[[2]] * c(-1, 1, 1), ~ stable_points(.x)) %>%
                       mutate(axis_space = map2(c("i", "ii", "iv"), c(1,2,2), rep) %>%
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
        facet_grid(~ axis_space) +
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
                   group_by(axis_space, spp) %>%
                   filter(time == max(time)) %>%
                   ungroup() %>%
                   filter(time < max(time)),
               shape = 4, size = 1.5) +
    facet_wrap(~ axis_space, nrow = 1) +
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
                               ts <- cond_coexist_df$axis_space %>%
                                   unique() %>%
                                   sort() %>%
                                   .[i]
                               stable_points((etas[[2]] * c(-1,1,0,1,0))[i]) %>%
                                   mutate(axis_space = ts)
                           }) %>%
    mutate(time = max(cond_coexist_df$time[cond_coexist_df$time < 7e3]),
           shp = case_when(V1 > 0 & V2 > 0 ~ 3L,
                           V1 == 0 ~ 2L,
                           TRUE ~ 1L) %>%
               factor(levels = 1:3)) %>%
    gather(axis, value, V1:V2) %>%
    mutate(axis = gsub("^V", "axis ", axis))



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
    gather(axis, value, V1:V2) %>%
    mutate(axis = gsub("^V", "axis ", axis)) %>%
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
                   gather(axis, value, V1:V2) %>%
                   mutate(axis = gsub("^V", "axis ", axis)) %>%
                   group_by(axis_space, spp) %>%
                   filter(time == max(time)) %>%
                   ungroup() %>%
                   filter(time < max(time)),
               aes(color = spp), shape = 4, size = 1.5) +
    facet_grid(axis ~ axis_space) +
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
#' They're all stable except for sub-additive and conflicting axis 1,
#' which is neutrally stable (see `_main-results__stability.R`).
#'




# =============================================================================*
# =============================================================================*

# <not updated> Figs 5,S2-S3: Conditional coexistence ----

# =============================================================================*
# =============================================================================*




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
            select(axis_space, rep, time, spp, N)
        .dd %>%
            mutate(rep = rep %>% paste() %>% as.integer()) %>%
            bind_rows(.dd2) %>%
            mutate(rep = factor(rep, levels = 0:12)) %>%
            mutate(id = interaction(spp, rep)) %>%
            ggplot(aes(time / 1000L, N, color = spp)) +
            geom_line(aes(group = id)) +
            geom_point(data = .dd %>%
                           group_by(axis_space, rep, spp) %>%
                           filter(time == max(time)) %>%
                           ungroup() %>%
                           filter(time < max(time)),
                       shape = 4, size = 1.5) +
            facet_grid(rep ~ axis_space) +
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
            group_by(axis_space) %>%
            mutate(tN = sqrt(N) / sum(sqrt(N))) %>%
            ungroup() %>%
            mutate(rep = 0L) %>%
            select(axis_space, rep, spp, N, tN)

        .dd3 <- map_dfr(list(.dd, .dd2),
                        ~ .x %>%
                            group_by(axis_space, rep) %>%
                            summarize(n_spp = sum(N > 0)) %>%
                            ungroup() %>%
                            mutate(rep = factor(rep %>% paste() %>% as.integer(),
                                                levels = 0:12)))

        .dd %>%
            group_by(axis_space, rep) %>%
            mutate(tN = sqrt(N) / sum(sqrt(N))) %>%
            # mutate(tN = N / sum(N)) %>%
            ungroup() %>%
            select(axis_space, rep, spp, N, tN) %>%
            mutate(rep = rep %>% paste() %>% as.integer()) %>%
            bind_rows(.dd2) %>%
            mutate(rep = factor(rep, levels = 0:12)) %>%
            ggplot(aes(rep)) +
            geom_bar(aes(weight = tN, fill = spp), color = "white", size = 0.1) +
            geom_text(data = .dd3 %>% filter(n_spp > 1),
                      aes(label = n_spp, y = 1.05),
                      size = 8 / 2.83465) +
            geom_vline(xintercept = 1.5, size = 0.5) +
            facet_grid( ~ axis_space) +
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
                                      ggtitle("With axis stochasticity") +
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
    filter(axis_space == "v") %>%
    mutate(id = interaction(spp, rep)) %>%
    arrange(time) %>%
    ggplot(aes(V1, V2, color = spp)) +
    geom_path(aes(group = id)) +
    geom_point(data = cond_coexist_stoch_df %>%
                   filter(d1 == "conflicting",
                          sigma_N == 0,
                          sigma_V == 0.1) %>%
                   filter(axis_space == "v") %>%
                   group_by(spp) %>%
                   summarize(V1 = V1[time == min(time)][1],
                             V2 = V2[time == min(time)][1])) +
    facet_wrap(~ rep, nrow = 3) +
    coord_equal(xlim = c(0, 3.11), ylim = c(0, 3.11)) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    ylab("(ameliorative)\nAxis 2") +
    xlab("Axis 1\n(conflicting)")





qg <- quant_gen(q = 2, eta = 0, d = c(-0.1, 0.1),
                n_reps = 24, n = 10,
                sigma_V = 0.1,
                spp_gap_t = 500L,
                final_t = 20e3L,
                add_var = rep(0.05, 10),
                n_threads = 3,
                show_progress = TRUE)

nv <- qg %>%
    .[["nv"]] %>%
    mutate(axis = paste0("V", axis)) %>%
    nest(value = c(geno, pheno)) %>%
    spread(axis, value) %>%
    unnest(c(V1, V2), names_sep = "_") %>%
    rename(V1 = V1_geno,
           V2 = V2_geno,
           Vp1 = V1_pheno,
           Vp2 = V2_pheno)

nv %>%
    mutate(id = interaction(spp, rep)) %>%
    arrange(time) %>%
    ggplot(aes(V1, V2, color = spp)) +
    # ggplot(aes(Vp1, Vp2, color = spp)) +
    geom_path(aes(group = id)) +
    geom_point(data = nv %>%
                   filter(!is.na(spp)) %>%
                   group_by(spp) %>%
                   summarize(V1 = V1[time == min(time)][1],
                             V2 = V2[time == min(time)][1])) +
                   # summarize(Vp1 = Vp1[time == min(time)][1],
                   #           Vp2 = Vp2[time == min(time)][1])) +
    facet_wrap(~ rep, nrow = 4) +
    coord_equal(xlim = c(0, 3.11), ylim = c(0, 3.11)) +
    scale_color_viridis_d(option = "D", guide = FALSE) +
    ylab("(ameliorative)\nAxis 2") +
    xlab("Axis 1\n(conflicting)")

nv %>%
    mutate(id = interaction(spp, rep)) %>%
    ggplot(aes(time / 1000, N, color = spp)) +
    geom_path(aes(group = id)) +
    facet_wrap(~ rep, nrow = 4) +
    scale_color_viridis_d(option = "D", guide = FALSE) +
    ylab("Abundance") +
    xlab("Time (× 1000)")



Z <- quant_gen(q = 2, eta = 0, d = c(-0.1, 0.1),
               n_reps = 24, n = 1,
               sigma_V = c(0.1, 0.025),
               final_t = 50e3L,
               # save_every = 10,
               save_every = 0,
               add_var = rep(0.05, 1),
               n_threads = 3)
znv <- Z$nv %>%
    mutate(axis = paste0("V", axis)) %>%
    nest(value = c(geno, pheno)) %>%
    spread(axis, value) %>%
    unnest(c(V1, V2), names_sep = "_") %>%
    rename(V1 = V1_geno,
           V2 = V2_geno,
           Vp1 = V1_pheno,
           Vp2 = V2_pheno)
znv %>%
    # filter(time == max(time)) %>%
    filter(!is.na(V1)) %>%
    ggplot(aes(V1, V2)) +
    geom_point(aes(color = spp)) +
    coord_equal(xlim = c(0, 3.11), ylim = c(0, 3.11)) +
    scale_color_viridis_d(option = "D", guide = FALSE) +
    ylab("(ameliorative)\nAxis 2") +
    xlab("Axis 1\n(conflicting)")

# znv %>%
#     mutate(id = interaction(spp, rep)) %>%
#     arrange(time) %>%
#     ggplot(aes(V1, V2, color = spp)) +
#     # ggplot(aes(Vp1, Vp2, color = spp)) +
#     geom_path(aes(group = id)) +
#     geom_point(data = nv %>%
#                    filter(!is.na(spp)) %>%
#                    group_by(spp) %>%
#                    summarize(V1 = V1[time == min(time)][1],
#                              V2 = V2[time == min(time)][1])) +
#     # summarize(Vp1 = Vp1[time == min(time)][1],
#     #           Vp2 = Vp2[time == min(time)][1])) +
#     facet_wrap(~ rep, nrow = 4) +
#     coord_equal(xlim = c(0, 3.11), ylim = c(0, 3.11)) +
#     scale_color_viridis_d(option = "D", guide = FALSE) +
#     ylab("(ameliorative)\nAxis 2") +
#     xlab("Axis 1\n(conflicting)")









if (.RESAVE_PLOTS) {
    save_plot(cond_coexist_stoch_p, 6, 4, .prefix = "S1-")
    save_plot(cc_sigmaV_sit_v_p, 6, 5, .prefix = "S2-")
}








# =============================================================================*
# =============================================================================*

# <not updated> Fig S3 Stoch. - # species ----

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

# <not updated> Additive --> Sub-additive ?? ----

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
               "there's only one stable axis state.",
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
    mutate(axis = paste0("V", axis)) %>%
    nest(value = c(geno, pheno)) %>%
    spread(axis, value) %>%
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
