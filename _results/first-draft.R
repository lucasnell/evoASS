
# load packages
suppressPackageStartupMessages({
    library(sauron)
    library(tidyverse)
    library(grid)
    library(gridExtra)
    library(pbmcapply)
    library(parallel)
    library(egg)
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



save_plot <- function(plot_obj, .width, .height, .prefix = NULL, .suffix = NULL) {
    fn <- gsub("_p$", "", paste(substitute(plot_obj)))
    if (!is.null(.prefix)) fn <- paste0(.prefix, fn)
    if (!is.null(.suffix)) fn <- paste0(fn, .suffix)
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
                 spp_gap_t = 1000L, final_t = 10e3L, save_every = 0L,
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




# Takes ~5 sec with q=2 and 3 threads
set.seed(145746935)
eta_sims <- list(-1, 0, 1) %>%
    map(one_eta_combo)
eta_sim_df <- map_dfr(eta_sims, ~.x$ts)




#' #'
#' #' Based on these eigenvalues...
#' #'   * When the tradeoff is additive, the state is neutrally stable
#' #'   * Everything else is stable
#' #'
#' #'
#' eta_sim_eigs <- map_dfr(1:length(eta_sims),
#'                         function(i) {
#'                             eigs <- map_dbl(eta_sims[[i]][["jacs"]],
#'                                             function(.x){
#'                                                 if (any(is.na(.x))) return(NA)
#'                                                 max(eigen(.x)$values)
#'                                             })
#'                             eta_sims[[i]]$ts %>%
#'                                 distinct(eta1) %>%
#'                                 mutate_all(sign) %>%
#'                                 rename_all(~ gsub("^eta", "sign", .x)) %>%
#'                                 mutate(e_min = min(eigs),
#'                                        e_max = max(eigs))
#'                         })
#'
#'
#' print_big_nums(eta_sim_eigs)




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



if (.RESAVE_PLOTS) save_plot(outcomes_q2_p, 5, 2, "1-")





# =============================================================================*
# =============================================================================*

# 2. Conditional coexistence ----

# =============================================================================*
# =============================================================================*


.eta <- 0.6

D <- matrix(0, 2, 2)
diag(D) <- 0.1
C <- matrix(.eta, 2, 2)
diag(C) <- 1
f <- formals(quant_gen)[["f"]]
a0 <- formals(quant_gen)[["a0"]]
r0 <- formals(quant_gen)[["r0"]]

.dist <- 0.6

V <- stable_points(.eta) %>%
    filter(V2 > 0) %>%
    add_row(V1 = .$V1 + c(0, .dist / sqrt(2), .dist,
                          .dist / sqrt(2), 0),
            V2 = .$V2 + c(-.dist, -.dist / sqrt(2), 0,
                          .dist / sqrt(2), .dist)) %>%
    .[-1,] %>%
    arrange(V2) %>%
    as.matrix() %>%
    t()

for (i in 2:5) {
    V0 <- do.call(cbind, c(rep(list(stable_points(.eta)[1,] %>%
                                        as.matrix() %>% t()),
                               i - 1),
                           list(V[,i])))
    if (i > 4) {
        V0 <- V0[,-3:-4]
    } else if (i > 3) V0 <- V0[,-3]
    F_ <- sauron:::F_t_cpp(V0,
                           # rep(1, ncol(V)),
                           c(pop_sizes(ncol(V0) - 1, .eta, 0.1), 1),
                           # 1,
                           f, a0, C, r0, D)
    print(F_)
    print(F_[length(F_)] - F_[length(F_)-1])
}; rm(i, F_)


sauron:::F_t_cpp(cbind(c(0, 2), c(0.424, 1.58)),
                       c(pop_sizes(1, 0.6, 0.1), 1),
                       f, a0, C, r0, D)



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

    if (.eta < 0) {
        if (.V0 == "unrestricted") {
            .V0 <- stable_points(.eta) %>%
                add_row(V1 = .$V1 + c(-.dist / sqrt(2), .dist / sqrt(2),
                                      -.dist / sqrt(2),
                                      0, .dist / sqrt(2)),
                        V2 = .$V2 + c(-.dist / sqrt(2), -.dist / sqrt(2),
                                      .dist / sqrt(2),
                                      .dist, .dist / sqrt(2))) %>%
                .[-1,] %>%
                as.matrix() %>%
                t()
        } else {
            .V0 <- stable_points(.eta) %>%
                add_row(V1 = c(.$V1 + c(-.dist / sqrt(2), -.dist,
                                        -.dist / sqrt(2),
                                        0, .dist / sqrt(2))),
                        V2 = .$V2 + c(-.dist / sqrt(2), 0, .dist / sqrt(2),
                                      .dist, .dist / sqrt(2))) %>%
                .[-1,] %>%
                as.matrix() %>%
                t()
        }
    } else {
        # if (.V0 == "unrestricted") {
        #     .V0 <- stable_points(.eta) %>%
        #         add_row(V1 = c(.$V1[.$V2 > 0] + c(0, .dist, .dist / sqrt(2), 0),
        #                        .$V1[.$V1 > 0]),
        #                 V2 = c(.$V2[.$V2 > 0] + c(-.dist, 0, .dist / sqrt(2),
        #                                           .dist),
        #                        .$V2[.$V1 > 0] + .dist)) %>%
        #         .[-1:-2,] %>%
        #         arrange(V2) %>%
        #         as.matrix() %>%
        #         t()
        # } else {
        #     .V0 <- stable_points(.eta) %>%
        #         filter(V2 > 0) %>%
        #         add_row(V1 = .$V1 + c(0, .dist / sqrt(2), .dist,
        #                               .dist / sqrt(2), 0),
        #                 V2 = .$V2 + c(-.dist, -.dist / sqrt(2), 0,
        #                               .dist / sqrt(2), .dist)) %>%
        #         .[-1,] %>%
        #         arrange(V2) %>%
        #         as.matrix() %>%
        #         t()
        # }
        .V0 <- rbind(seq(1.2, 0, length.out = 5),
                     seq(1.3, 2.5, length.out = 5))
        if (.lab == "unrestricted") {
            v2 <- .V0[2, which_switch]
            .V0[2, which_switch] <- .V0[1, which_switch]
            .V0[1, which_switch] <- v2
        }
    }

    Z <- quant_gen(q = .q, eta = .eta, d = .ds,
                   n_reps = 1, n = ncol(.V0),
                   V0 = .V0,
                   sigma_V0 = 0,
                   spp_gap_t = 5e3L,
                   add_var = rep(0.05, ncol(.V0)),
                   final_t = 20e3L,
                   show_progress = FALSE) %>%
        .[["nv"]] %>%
        mutate(trait = paste0("V", trait)) %>%
        spread(trait, geno) %>%
        select(-rep) %>%
        mutate(V0 = .lab, eta = .eta, d1 = .d1)

    return(Z)
}




Z <- quant_gen(q = 2, eta = .eta, d = 0.1,
               n_reps = 1, n = 1,
               V0 = rbind(2, 2 + 1e-6),
               N0 = 1, add_var = 0.1,
               sigma_V0 = 0,
               spp_gap_t = 0L,
               final_t = 5e3L,
               save_every = 1L,
               show_progress = FALSE) %>%
    .[["nv"]] %>%
    mutate(trait = paste0("V", trait)) %>%
    spread(trait, geno) %>%
    select(-rep)

Z %>% filter(time == 650)


Z %>%
    filter(time < 650, time >= 500) %>%
    gather(trait, value, V1:V2) %>%
    {
        ggplot(., aes(time, value, color = spp)) +
        geom_line() +
        geom_point(data = filter(., time %in% (min(time) + 25 * 1:5))) +
        ggtitle(Sys.time()) +
        facet_wrap(~ trait) +
        scale_color_brewer(palette = "Dark2", guide = FALSE) +
        ylab("Trait value") +
        xlab("Time")
    }


Z %>%
    filter(time %in% round(750 * (0.9^(0:4)))) %>%
    select(starts_with("V", FALSE)) %>%
    as.matrix() %>%
    t()


Z %>%
    filter(time %in% (500 + 25 * 1:5)) %>%
    select(starts_with("V", FALSE)) %>%
    lm(formula = V2 ~ V1)



ggarrange(Z %>%
              arrange(time, spp) %>%
              ggplot(aes(V1, V2, color = spp)) +
              geom_polygon(data = tibble(x = c(-1, 5, -1), y = c(-1, 5, 5)), aes(x,y),
                           color = NA, fill = "gray80") +
              ggtitle(Sys.time()) +
              geom_path() +
              geom_point(data = Z %>%
                             filter(time %in% (500 + 25 * 1:5))) +
              geom_point(data = Z %>%
                             filter(time == min(time)),
                         shape = 1) +
              geom_point(data = Z %>%
                             filter(time == max(time)),
                         shape = 4, size = 4, color = "black") +
              coord_equal(xlim = c(0, 4.25), ylim = c(0, 4.25)) +
              scale_color_brewer(palette = "Dark2", guide = FALSE) +
              ylab("Trait 2") +
              xlab("Trait 1"),
          Z %>%
              filter(time < 1e3L) %>%
              gather(trait, value, V1:V2) %>%
              ggplot(aes(time, value, color = spp)) +
              geom_line() +
              geom_line(data = Z %>%
                            filter(time < 1e3L),
                        aes(time, N / 5), linetype = 2) +
              ggtitle(Sys.time()) +
              facet_wrap(~ trait) +
              scale_color_brewer(palette = "Dark2", guide = FALSE) +
              ylab("Trait value") +
              xlab("Time"),
          nrow = 1)






invader_fitness <- function(.V1, .V2, .res_n_spp) {
    V0_ <- do.call(cbind, c(rep(list(c(0, 2)), .res_n_spp), list(c(.V1, .V2))))
    F_ <- sauron:::F_it_cpp(i = .res_n_spp,
                            V = V0_,
                            N = c(pop_sizes(.res_n_spp, C[1,2], D[1,1]), 1),
                            f, a0, C, r0, D)
    return(tibble(res_n_spp = .res_n_spp, V1 = .V1, V2 = .V2, fitness = F_))
}


invader_heatmap <- crossing(.V1 = seq(0, 3, 0.1),
         .V2 = seq(0, 3, 0.1),
         .res_n_spp = 1:4) %>%
    pmap_dfr(.f = invader_fitness)


invader_heatmap %>%
    mutate(res_n_spp = factor(res_n_spp)) %>%
    ggplot(aes(V1, V2, fill = fitness)) +
    geom_tile() +
    # geom_contour(aes(z = fitness), color = "black", size = 0.25) +
    geom_point(data = invader_heatmap %>%
                   mutate(res_n_spp = factor(res_n_spp)) %>%
                   filter(fitness > 0.8), shape = 1) +
    scale_fill_gradient2(midpoint = 1) +
    facet_wrap(~ res_n_spp, nrow = 2) +
    coord_equal()


invader_heatmap %>%
    group_by(res_n_spp) %>%
    summarize(V1 = V1[fitness == max(fitness)],
              V2 = V2[fitness == max(fitness)])

invader_heatmap %>%
    filter(res_n_spp == 3, fitness > 0.8) %>%
    # filter(abs(V1 - V2) == min(abs(V1 - V2)))
    filter(V2 == max(V2))







cond_coexist_df <- crossing(.V0 = c("restricted", "unrestricted"),
         # .eta_sign = c(-1,1),
         .eta_sign = 1,
         # .d1 = 0.1) %>%
         .d1 = c(0.1, -0.1)) %>%
    pmap_dfr(cond_coexist_test) %>%
    mutate(V0 = factor(V0, levels = c("restricted", "unrestricted")),
           eta = factor(eta, levels = c(-1, 1) * etas[[2]],
                        labels = c("sub-additive", "super-additive")))


cond_coexist_df %>%
    filter(time == max(time)) %>%
    group_by(d1, eta, V0) %>%
    summarize(n_spp = n()) %>%
    ungroup()


# cond_coexist_df %>%
#     filter(d1 > 0, eta == "super-additive", V0 == "restricted", time > 5e3L-10)



# Time series for coexistence and exclusion
# cond_coexist_ts_p <-
cond_coexist_df %>%
    mutate(d1 = factor(d1)) %>%
    # filter(d1 > 0) %>%
    ggplot(aes(time / 1000L, N, color = spp)) +
    geom_line() +
    # facet_grid(eta ~ V0, scales = "free_y") +
    facet_grid(d1 ~ V0, scales = "free_y") +
    # scale_color_viridis_d(begin = 0.1, end = 0.9, option = "C", guide = FALSE) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    ylab("Abundance") +
    scale_x_continuous(expression("Time (" %*% "1,000" * ")")) +
    theme(#strip.text.x = element_blank(),
          plot.margin = margin(0,0,0,t=20))



# Movement through trait space for coexistence and exclusion:
# cond_coexist_sp_p <-
    cond_coexist_df %>%
    mutate(d1 = factor(d1)) %>%
    # filter(d1 > 0) %>%
    # filter(time < 2.5e3) %>%
    arrange(time, spp) %>%
    ggplot(aes(V1, V2, color = spp)) +
    geom_polygon(data = tibble(x = c(-1, 5, -1), y = c(-1, 5, 5)), aes(x,y),
                 color = NA, fill = "gray80") +
    geom_path() +
    geom_point(data = cond_coexist_df %>%
                   mutate(d1 = factor(d1)) %>%
                   group_by(d1, eta, V0, spp) %>%
                   filter(time == min(time)) %>%
                   ungroup(),
               shape = 1) +
    geom_point(data = cond_coexist_df %>%
                   mutate(d1 = factor(d1)) %>%
                   filter(time == max(time)) %>%
                   group_by(d1, eta, V0) %>%
                   filter(unq_spp_filter(V1, V2)) %>%
                   ungroup(),
               shape = 4, size = 4, color = "black") +
    # scale_color_viridis_d(begin = 0.1, end = 0.9, option = "C", guide = FALSE) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    facet_grid(d1 ~ V0) +
    coord_equal(xlim = c(0, 4.25), ylim = c(0, 4.25)) +
    ylab("Trait 2") +
    xlab("Trait 1") +
    theme(# strip.text = element_blank(),
          plot.margin = margin(0,0,0,t=20))



cond_coexist_p_label <-
    tableGrob(cbind("exclusion", "coexistence"),
              theme = ttheme_minimal(core = list(fg_params =
                                                     list(hjust=0.5, x=0.5,
                                                          vjust=0.5, y=0.5))),
              widths = unit(rep(1/2,2), "npc"),
              heights = unit(2,"npc"))


cond_coexist_p <- ggarrange(cond_coexist_ts_p,
                            cond_coexist_sp_p,
                            heights = c(0.4, 0.6),
                            ncol = 1, labels = LETTERS[1:2],
                            top = cond_coexist_p_label,
                            label.args = list(gp = gpar(fontsize = 14,
                                                        fontface =  "plain"),
                                              vjust = 1.7, hjust = -3))


if (.RESAVE_PLOTS) save_plot(cond_coexist_p, 4.5, 5, "2-")



# # Are these stable?
#
# cond_coexist_test_stability <- function(.V0, .lab, .ds) {
#     .q <- 2
#     stopifnot(length(.ds) == .q)
#     .ds <- abs(.ds)*c(-1,rep(1,.q-1))
#     Z <- quant_gen(q = .q, eta = etas[[2]], d = .ds, max_t = 20e3L, n_reps = 1,
#                    save_every = 0L, n = 100, N0 = rep(1, 100),
#                    V0 = .V0,
#                    start_t = 0, sigma_V0 = 0,
#                    show_progress = FALSE)
#     Z$call[["q"]] <- .q
#     Z$call[["d"]] <- .ds
#     Z$call[["etas"]] <- etas[[2]]
#     return(jacobians(Z)[[1]])
# }
#
# cond_coexist_J <- tibble(.V0 = list(V0_exclude, V0_coexist),
#                          .lab = c("exclusion", "coexistence"),
#                          .ds = list(c(0.1, 0.1))) %>%
#     pmap(cond_coexist_test_stability)
#
# map(cond_coexist_J, ~ max(eigen(.x, only.values = TRUE)[["values"]]))













#' # =============================================================================*
#' # =============================================================================*
#'
#' # 3. 3-trait outcomes ----
#'
#'
#' # =============================================================================*
#' # =============================================================================*
#'
#' to_cap <- "Unique trait values for surviving species, for all 27 combinations
#'            of pairwise combinations of traits being
#'            sub-additive ($\\eta < 0$ or \"$-$\"), neutral ($\\eta = 0$),
#'            or super-additive ($\\eta > 0$ or \"$+$\").
#'            The size of points indicates the value for trait 3.
#'            Sub-panels separate the number of pairwise combinations (out of 3)
#'            that are sub-additive, neutral, or super-additive.
#'            All but the \"all super\" and \"all neutral\" sub-panels
#'            contain multiple pairwise combinations, so each specific combination
#'            is separated by color and indicated by the colored label.
#'            The symbol $0^{+}$ inside labels indicates when that $\\eta$ can be
#'            $\\ge 0$.
#'            Therefore, all 3 points in the \"1 sub\" sub-panel represent
#'            4 combinations each.
#'            When one sub-panel contains multiple points of the same color,
#'            this indicates that multiple alternative states are possible
#'            for that combination."
#'
#'
#'
#' if (.REDO_SIMS) {
#'     # Takes ~2.5 min with q=3 and 3 threads
#'     set.seed(1611377415)
#'     eta_sims_q3 <- crossing(eta1 = -1:1, eta2 = -1:1,
#'                             eta3 = -1:1) %>%
#'         split(1:nrow(.)) %>%
#'         map(as.numeric) %>%
#'         map(one_eta_combo, .d = 1, max_t = 50e3L)
#'     saveRDS(eta_sims_q3, rds("eta_sim-q3"))
#'     eta_sim_q3_df <- map_dfr(eta_sims_q3, ~.x$ts)
#'
#'     # # Takes >>18 min (didn't finish) with q=4 and 3 threads
#'     # t0 <- Sys.time()
#'     # set.seed(266613920)
#'     # eta_sims_q4 <- crossing(eta1 = -1:1, eta2 = -1:1,
#'     #                         eta3 = -1:1, eta4 = -1:1,
#'     #                         eta5 = -1:1, eta6 = -1:1) %>%
#'     #     split(1:nrow(.)) %>%
#'     #     map(as.numeric) %>%
#'     #     map(one_eta_combo, .d = 1, max_t = 50e3L)
#'     # saveRDS(eta_sims_q4, rds("eta_sim-q4"))
#'     # t1 <- Sys.time()
#'     # t1 - t0; # rm(t0, t1)
#'     # eta_sim_q4_df <- map_dfr(eta_sims_q4, ~.x$ts)
#' } else {
#'     eta_sims_q3 <- readRDS(rds("eta_sim-q3"))
#'     eta_sim_q3_df <- map_dfr(eta_sims_q3, ~.x$ts)
#'     # eta_sims_q4 <- readRDS(rds("eta_sim-q4"))
#'     # eta_sim_q4_df <- map_dfr(eta_sims_q4, ~.x$ts)
#' }
#'
#'
#'
#' #'
#' #' Based on these eigenvalues...
#' #'   *
#' #'
#' #'
#' eta_sims_q3_eigs <- map_dfr(1:length(eta_sims_q3),
#'                         function(i) {
#'                             eigs <- map_dbl(eta_sims_q3[[i]][["jacs"]],
#'                                             function(.x){
#'                                                 if (any(is.na(.x))) return(NA)
#'                                                 max(eigen(.x)$values)
#'                                             })
#'                             eta_sims_q3[[i]]$ts %>%
#'                                 distinct(eta1, eta2, eta3) %>%
#'                                 mutate_all(sign) %>%
#'                                 rename_all(~ gsub("^eta", "sign", .x)) %>%
#'                                 mutate(e_min = min(eigs),
#'                                        e_max = max(eigs))
#'                         })
#'
#'
#' print_big_nums(eta_sims_q3_eigs, n = 30)
#'
#' eta_sims_q3_eigs %>%
#'     mutate_at(vars(starts_with("sign")),
#'               ~ case_when(.x < 0 ~ "-",
#'                           .x > 0 ~ "+",
#'                           TRUE ~ "0")) %>%
#'     mutate_at(vars(starts_with("e_m")),
#'               ~ case_when(.x < 1 ~ "stable",
#'                           .x > 1 ~ paste(.x - 1),
#'                           TRUE ~ "neutrally stable")) %>%
#'     mutate(id = paste0(sign1, sign2, sign3) %>% factor()) %>%
#'     select(id, e_min, e_max) %>%
#'     print(n = 30)
#'
#'
#' eta_sim_q3_df %>%
#'     group_by(eta1, eta2, eta3) %>%
#'     summarize(N = n()) %>%
#'     ungroup() %>%
#'     print(n = 30)
#'
#' parseable <- function(.id) {
#'     .id = paste(.id)
#'     IDD1 = str_sub(.id, 1, 1) %>%
#'         {case_when(. == "+" ~ paste0("{}", .),
#'                    . == "-" ~ paste0(., "''"),
#'                    . == "0" ~ "0",
#'                    . == "X" ~ "0^{'+'}",
#'                    TRUE ~ NA_character_)}
#'     IDD2 = str_sub(.id, 2, 2) %>%
#'         {case_when(. == "-" ~ paste0(., "''"),
#'                    . == "X" ~ "0^{'+'}",
#'                    TRUE ~ .)}
#'     IDD3 = str_sub(.id, 3, 3) %>%
#'         {case_when(. == "+" ~ paste0(., "{}"),
#'                    . == "-" ~ paste0(., "''"),
#'                    . == "0" ~ "0",
#'                    . == "X" ~ "0^{'+'}",
#'                    TRUE ~ NA_character_)}
#'     if(any(is.na(IDD1)) || any(is.na(IDD3))) {
#'         print(.id[is.na(IDD1) | is.na(IDD3)])
#'         stop("Non-programmed option in `parseable`")
#'     }
#'     paste0(IDD1, IDD2, IDD3) %>%
#'         str_replace_all("''0", "0") %>%
#'         str_replace_all("\\{\\}0", "0") %>%
#'         str_replace_all("\\}0", "\\} * 0") %>%
#'         str_replace_all("00", "0 * 0") %>%
#'         str_replace_all("00", "0 * 0")
#' }
#'
#'
#' n_states_q3_df <- eta_sim_q3_df %>%
#'     group_by(eta1, eta2, eta3) %>%
#'     summarize(N = n()) %>%
#'     ungroup() %>%
#'     mutate_at(vars(starts_with("eta")),
#'               ~ case_when(.x < 0 ~ "-",
#'                           .x > 0 ~ "+",
#'                           TRUE ~ "0")) %>%
#'     mutate(id = paste0(eta1, eta2, eta3) %>% parseable() %>% factor(),
#'            id = reorder(id, N))
#'
#' outcomes_q3_p <- eta_sim_q3_df %>%
#'     mutate_at(vars(starts_with("eta")),
#'               ~ case_when(.x < 0 ~ "-",
#'                           .x > 0 ~ "+",
#'                           TRUE ~ "0")) %>%
#'     mutate(id = paste0(eta1, eta2, eta3) %>% parseable() %>%
#'                factor(levels = n_states_q3_df$id %>% levels())) %>%
#'     group_by(id) %>%
#'     filter(unq_spp_filter(V1, V2, V3, .prec = 0.01)) %>%
#'     ungroup() %>%
#'     ggplot() +
#'     geom_point(aes(V1, V2, size = V3), shape = 1, color = "firebrick") +
#'     geom_text(data = n_states_q3_df,
#'               aes(x = 2.5, y = 2.5, label = N),
#'               size = 10 / 2.83465, fontface = "bold",
#'               hjust = 1, vjust = 1) +
#'     scale_size_continuous("Trait 3", range = c(1, 5), breaks = 0.5 * 1:3) +
#'     scale_x_continuous("Trait 1", breaks = 0:2) +
#'     scale_y_continuous("Trait 2", breaks = 0:2) +
#'     coord_equal(xlim = c(0, 2.5), ylim = c(0, 2.5)) +
#'     facet_wrap(~ id, nrow = 5, labeller = label_parsed) +
#'     theme(strip.text = element_text(size = 9),
#'           panel.border = element_rect(size = 0.5, fill = NA)) +
#'     NULL
#'
#'
#'
#' if (.RESAVE_PLOTS) save_plot(outcomes_q3_p, 6.5, 6, "3-")
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
#' # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
#' # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
#'
#' # Interactions between d and eta ----
#'
#' # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
#' # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*
#'
#' ## The code below explores what happens when d >= 0.9
#'
#' # eta_sim_df %>%
#' #     filter((eta1 < 0 & eta2 == 0 & eta3 > 0) |
#' #                (eta1 < 0 & eta2 > 0 & eta3 == 0) |
#' #                (eta1 == 0 & eta2 < 0 & eta3 > 0) |
#' #                (eta1 == 0 & eta2 > 0 & eta3 < 0) |
#' #                (eta1 > 0 & eta2 == 0 & eta3 < 0) |
#' #                (eta1 > 0 & eta2 < 0 & eta3 == 0)) %>%
#' #     group_by(eta1, eta2, eta3) %>%
#' #     summarize(N = n()) %>%
#' #     ungroup()
#'
#' eta_sim_df %>%
#'     filter(eta1 > 0)
#'
#'
#' sign1 = 1; sign2 = -1; sign3 = 1
#' C <- matrix(0, .q, .q)
#' C[lower.tri(C)] <- abs(etas) * c(sign1, sign2, sign3)[1:length(etas)]
#' C <- C + t(C)
#' diag(C) <- 1
#'
#' trait_to <- quant_gen(q = .q, eta = C, d = 10, max_t = 20e3L, n_reps = 24,
#'                       save_every = 100L, n = 100, N0 = rep(1, 100),
#'                       start_t = 0, sigma_V0 = 2, n_threads = .N_THREADS,
#'                       show_progress = TRUE)
#'
#' trait_to$nv %>%
#'     filter(time == max(time)) %>%
#'     mutate(trait = paste0("V", trait)) %>%
#'     spread(trait, value) %>%
#'     filter(unq_spp_filter(V1, V2)) %>%
#'     select(starts_with("V")) %>%
#'     # mutate(eta1 = C[2,1], eta2 = C[3,1], eta3 = C[3,2]) %>%
#'     identity()
#'
#' trait_to$nv %>%
#'     filter(time == max(time)) %>%
#'     mutate(trait = paste0("V", trait)) %>%
#'     spread(trait, value) %>%
#'     # filter(unq_spp_filter(V1, V2, V3)) %>%
#'     # select(starts_with("V")) %>%
#'     filter(!((V1 - 1.63)^2 < 0.001 & (V2 - 1.63)^2 < 0.001 & V3 == 0)) %>%
#'     .[["N"]] %>%
#'     range()
#'
#'
#'
#'
#' # # A tibble: 1 x 6
#' #      V1    V2    V3   eta1  eta2  eta3
#' #   <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl>
#' # 1  1.63  1.63     0 -0.207     0 0.257
#'
#' trait_to$nv %>%
#'     filter(trait == 1, rep == 1) %>%
#'     # filter(time == max(time))
#'     # filter(time < 1000) %>%
#'     ggplot(aes(time, N)) +
#'     geom_line(aes(color = spp)) +
#'     scale_color_viridis_d(guide = FALSE)
#'
#'
#'
#' zz <- perturb(trait_to, 100e3L, 0, d = 0.27)
#'
#' # all.equal(zz$start$N, zz$end$N)
#' all.equal(zz$start$V, zz$end$V)
#'
#'
#' # LEFT OFF --> even if you reduce d to 0.27 after the fact, the traits remain
#' # stable (Ns don't remain stable with even small change in d, as expected)
#'
#'
#'
#'
#'
#'
#'
#'
#' # This function combines id labels for all combinations that have the exact
#' # same outcomes. Only used for "1 sub".
#' # It also processes all id columns to be parsed.
#' process_ids <- function(z) {
#'
#'     parseable <- function(.id) {
#'         .id = paste(.id)
#'         IDD1 = str_sub(.id, 1, 1) %>%
#'             {case_when(. == "+" ~ paste0("{}", .),
#'                        . == "-" ~ paste0(., "''"),
#'                        . == "0" ~ "0",
#'                        . == "X" ~ "0^{'+'}",
#'                        TRUE ~ NA_character_)}
#'         IDD2 = str_sub(.id, 2, 2) %>%
#'             {case_when(. == "-" ~ paste0(., "''"),
#'                        . == "X" ~ "0^{'+'}",
#'                        TRUE ~ .)}
#'         IDD3 = str_sub(.id, 3, 3) %>%
#'             {case_when(. == "+" ~ paste0(., "{}"),
#'                        . == "-" ~ paste0(., "''"),
#'                        . == "0" ~ "0",
#'                        . == "X" ~ "0^{'+'}",
#'                        TRUE ~ NA_character_)}
#'         if(any(is.na(IDD1)) || any(is.na(IDD3))) {
#'             print(.id[is.na(IDD1) | is.na(IDD3)])
#'             stop("Non-programmed option in `parseable`")
#'         }
#'         paste0(IDD1, IDD2, IDD3) %>%
#'             str_replace_all("''0", "0") %>%
#'             str_replace_all("\\{\\}0", "0") %>%
#'             str_replace_all("\\}0", "\\} * 0") %>%
#'             str_replace_all("00", "0 * 0") %>%
#'             str_replace_all("00", "0 * 0")
#'     }
#'
#'
#'     if (z$p_groups[1] != "'1 sub'") return(mutate(z, id = parseable(id)))
#'
#'     one_chr <- function(.x) {
#'         .xx <- .x %>%
#'             str_split("") %>%
#'             {lapply(1:length(.[[1]]), function(i) map_chr(., ~ .x[[i]]))} %>%
#'             map_chr(function(x) {
#'                 xx <- unique(x)
#'                 if (length(xx) == 1) return(xx)
#'                 if (all(xx %in% c("0", "+"))) return("X")
#'                 return(NA_character_)
#'             })
#'         if (any(is.na(.xx))) stop("Non-programmed input to `one_chr`")
#'         return(paste(.xx, collapse = ""))
#'     }
#'
#'     z %>%
#'         mutate(color_group = group_spp(V1, V2, V3) + 1L) %>%
#'         mutate_at(vars(id), paste) %>%
#'         arrange(p_groups, color_group) %>%
#'         group_by(p_groups, color_group) %>%
#'         summarize(V1 = mean(V1),
#'                   V2 = mean(V2),
#'                   V3 = mean(V3),
#'                   id = one_chr(id)) %>%
#'         ungroup() %>%
#'         mutate(id = parseable(id))
#'
#' }
#'
#' p_group_levels <- c("{} >= '2 sub'", "'1 sub'",
#'                     "'all super'",
#'                     "'2 super, 1 neutral'", "'2 neutral, 1 super'",
#'                     "'all neutral'")
#'
#'
#' trait_outcomes_p_df <- eta_sim_df %>%
#'     mutate_at(vars(starts_with("eta")),
#'               ~ factor(sign(.x), levels = -1:1, labels = c("-", "0", "+"))) %>%
#'     mutate(id = interaction(eta1, eta2, eta3, sep = "")) %>%
#'     group_by(id) %>%
#'     filter(unq_spp_filter(V1, V2, V3, .prec = 0.1)) %>%
#'     ungroup() %>%
#'     mutate(n_p = str_count(id, "\\+"),
#'            n_n = str_count(id, "\\-"),
#'            n_z = str_count(id, "0"),
#'            p_groups = case_when(n_n >= 2 ~ p_group_levels[1],
#'                                 n_n == 1  ~ p_group_levels[2],
#'                                 n_p == 3 ~ p_group_levels[3],
#'                                 n_p == 2 & n_z == 1 ~ p_group_levels[4],
#'                                 n_z == 2 & n_p == 1 ~ p_group_levels[5],
#'                                 n_z == 3 ~ p_group_levels[6],
#'                                 TRUE ~ NA_character_)) %>%
#'     select(-starts_with("n", ignore.case = FALSE)) %>%
#'     # filter(is.na(p_groups))  ## <-- check that this yields no rows
#'     mutate(p_groups = factor(p_groups, levels = p_group_levels)) %>%
#'     arrange(p_groups) %>%
#'     split(.$p_groups) %>%
#'     map(~ mutate(.x, color_group = id %>% droplevels() %>% as.integer())) %>%
#'     map_dfr(process_ids) %>%
#'     mutate_at(vars(color_group, id), factor) %>%
#'     arrange(desc(V3))
#'
#'
#'
#' .pal <- function(.a = 1) viridisLite::plasma(7, end = 0.85,
#'                                              alpha = .a)[c(6,4,1,3,5,7,2)]
#'
#'
#' trait_outcomes_p <- trait_outcomes_p_df %>%
#'     ggplot(aes(V1, V2, size = V3, color = color_group)) +
#'     geom_point(aes(fill = color_group), shape = 21) +
#'     geom_text(data = trait_outcomes_p_df %>%
#'                   filter(!grepl("^\\'all", p_groups)) %>%
#'                   group_by(p_groups, color_group, id) %>%
#'                   summarize_at(vars(V1, V2), median) %>%
#'                   ungroup() %>%
#'                   mutate(V1 = case_when(id == "{}+0+{}" ~ V1 - 0.3,
#'                                         id == "0 * 0+{}" ~ V1 - 0.3,
#'                                         id == "{}+0 * 0" ~ V1 - 0.1,
#'                                         id == "0-''-''" ~ V1 - 0.4,
#'                                         id == "-''-''+{}" ~ V1 - 0.1,
#'                                         id == "{}+-''-''" ~ V1 - 0.3,
#'                                         id == "-''+-''" ~ V1 - 0.3,
#'                                         id == "-0-''" ~ V1 - 0.4,
#'                                         TRUE ~ V1),
#'                          V2 = case_when(id == "{}+0+{}" ~ V2 + 0.05,
#'                                         id == "0 * 0+{}" ~ V2 + 0.05,
#'                                         id == "{}+0 * 0" ~ V2 + 0.05,
#'                                         id == "0+0" ~ V2 - 0.05,
#'                                         id == "0-''-''" ~ V2 - 0.6,
#'                                         id == "-''-0" ~ V2 - 0.4,
#'                                         id == "-''-''+{}" ~ V2 - 0.3,
#'                                         id == "-''-''-''" ~ V2 - 0.1,
#'                                         id == "{}+-''-''" ~ V2 - 0.6,
#'                                         id == "-0-''" ~ V2 + 0.05,
#'                                         TRUE ~ V2)) %>%
#'                   identity(),
#'               aes(label = id), size = 9 / 2.835, hjust = 0, vjust = 0, parse = TRUE,
#'               fontface = "bold", lineheight = 0.75, nudge_x = 0.2, nudge_y = 0.2) +
#'     scale_size_continuous("Trait 3", range = c(2, 6), breaks = 0.5 * 1:3) +
#'     scale_x_continuous("Trait 1", breaks = 0:2) +
#'     scale_y_continuous("Trait 2", breaks = 0:2) +
#'     coord_equal(xlim = c(0, 2.5), ylim = c(0, 2.5)) +
#'     facet_wrap(~ p_groups, nrow = 3, labeller = label_parsed) +
#'     theme(strip.text = element_text(size = 10),
#'           panel.border = element_rect(size = 0.5, fill = NA)) +
#'     scale_color_manual(values = .pal(), guide = FALSE) +
#'     scale_fill_manual(values = .pal(0.25), guide = FALSE) +
#'     NULL
#'
#'
#' if (.RESAVE_PLOTS) save_plot(trait_outcomes_p, 6, 5, "1-")
#'
#' trait_outcomes_p






# =============================================================================*
# =============================================================================*

# S1. Coexistence ----

# =============================================================================*
# =============================================================================*



coex_cap <- "Number of surviving species for
             (A) all permutations of two traits being conflicting (``$-$'')
             or non-conflicting (``$+$''),
             (B) varying values of $d_2$ when $d_1$ is kept positive
             (i.e., trait 1 is kept non-conflicting), and
             (C) surviving species through time with $d_2 = -10^{-2}$ and
             $d_2 = -10^{-4}$."






# ------------------------------*
# __ all combinations ----
# ------------------------------*

# All combinations of - or + d values

one_d_combo <- function(sign1, sign2) {

    .d <- c(0.05, 0.1) * c(sign1, sign2)

    Z <- quant_gen(q = 2, eta = etas[[2]], d = .d, max_t = 50e3L, n_reps = 24,
                   save_every = 0L, n = 100, N0 = rep(1, 100),
                   start_t = 0, sigma_V0 = 2, n_threads = .N_THREADS,
                   show_progress = FALSE)

    Z$call[["d"]] <- eval(.d)

    jacs <- jacobians(Z)

    return(list(NV = Z$nv %>%
                    mutate(trait = paste0("V", trait)) %>%
                    spread(trait, value) %>%
                    # filter(unq_spp_filter(V1, V2, V3)) %>%
                    mutate(d1 = .d[1], d2 = .d[2]),
                J = jacs))


}




if (.REDO_SIMS) {
    # Takes ~10 sec w/ q = 2 and 3 threads
    set.seed(751678518)
    all_d_sim <- crossing(sign1 = c(-1,1), sign2 = c(-1,1)) %>%
        pmap(one_d_combo)
    saveRDS(all_d_sim, rds("all_d_sim-q2"))
    d_sim_df <- map_dfr(all_d_sim, ~ .x[["NV"]])
} else {
    all_d_sim <- readRDS(rds("all_d_sim-q2"))
    d_sim_df <- map_dfr(all_d_sim, ~ .x[["NV"]])
}


#' #'
#' #' Code below shows that they're all stable.
#' #' It also shows that we get complex eigenvalues for positive d1 and d2
#' #' in about 17% of reps.
#' #'
#' map_dfr(1:length(all_d_sim),
#'         function(i) {
#'             eigs <- map_dbl(all_d_sim[[i]][["J"]],
#'                             function(.x){
#'                                 if (any(is.na(.x))) return(NA)
#'                                 max(Re(eigen(.x)$values))
#'                             })
#'             all_d_sim[[i]]$NV %>%
#'                 .[1,] %>%
#'                 select(d1, d2) %>%
#'                 mutate(e_min = min(eigs),
#'                        e_max = max(eigs))
#'         }) %>%
#'     print_big_nums()
#'
#'
#' map_dfr(1:length(all_d_sim),
#'         function(i) {
#'             eigs <- map(all_d_sim[[i]][["J"]],
#'                         function(.x){
#'                             if (any(is.na(.x))) return(NA)
#'                             eigen(.x)$values
#'                         })
#'             eigs <- do.call(c, eigs)
#'             tibble(d1 = all_d_sim[[i]]$NV[["d1"]][[1]],
#'                    d2 = all_d_sim[[i]]$NV[["d2"]][[1]],
#'                    eigs = eigs)
#'         }) %>%
#'     filter(Im(eigs) != 0)
#'
#'
#' map_dfr(1:length(all_d_sim),
#'         function(i) {
#'             eigs <- map(all_d_sim[[i]][["J"]],
#'                         function(.x){
#'                             if (any(is.na(.x))) return(NA)
#'                             eigen(.x)$values
#'                         })
#'             cmplx <- map_lgl(eigs, ~ any(Im(.x) != 0))
#'             tibble(d1 = all_d_sim[[i]]$NV[["d1"]][[1]],
#'                    d2 = all_d_sim[[i]]$NV[["d2"]][[1]],
#'                    cmplx = mean(cmplx))
#'         }) %>%
#'     identity()





coexist_all_d_p <- d_sim_df %>%
    group_by(d1, d2, rep) %>%
    summarize(N = n()) %>%
    ungroup() %>%
    mutate_at(vars(starts_with("d")),
              ~ factor(sign(.x), levels = c(-1,1), labels = c("-", "+"))) %>%
    mutate(id = interaction(d1, d2, sep = "{}"),
           id = factor(id, levels = levels(id),
                       labels = sprintf("%s{}", levels(id))),
           n_p = str_count(id, "\\+"),
           idd = interaction(n_p, id, sep = "", lex.order = TRUE, drop = TRUE),
           idd = factor(idd, levels = levels(idd),
                        labels = gsub("0|1|2|3", "", levels(idd)))) %>%
    ggplot(aes(idd, N)) +
    geom_jitter(width = 0.25, alpha = 0.5, height = 0,
                color = "dodgerblue4", fill = "dodgerblue", shape = 21) +
    scale_y_continuous("Number of surviving species", limits = c(0, 30),
                       breaks = c(0, 15, 30)) +
    scale_x_discrete(expression("Sign of" ~ d[1] ~ "and" ~ d[2]),
                     labels = rlang::parse_exprs) +
    scale_color_viridis_d(guide = FALSE, end = 0.8, option = "B") +
    theme(strip.text = element_text(size = 8),
          axis.text.x = element_text(size = 11),
          plot.margin = margin(0,0,0,t=12)) +
    NULL




# ------------------------------*
# __ slowly vary 1 ----
# ------------------------------*

# Keep one d value positive, make the last one go slowly negative

one_d_combo_vary_one <- function(.d2, .max_t) {


    .d <- c(0.05, .d2)

    Z <- quant_gen(q = length(.d), eta = etas[[2]], d = .d, max_t = .max_t,
                   n_reps = 24, save_every = 0L, n = 100, N0 = rep(1, 100),
                   start_t = 0, sigma_V0 = 2, n_threads = .N_THREADS,
                   show_progress = FALSE)

    Z$call[["q"]] <- eval(length(.d))
    Z$call[["d"]] <- eval(.d)

    jacs <- jacobians(Z)

    return(list(NV = Z$nv %>%
                    spread(trait, value) %>%
                    group_by(rep) %>%
                    summarize(n_spp = n()) %>%
                    ungroup() %>%
                    mutate(d2 = .d2) %>%
                    identity(),
                J = jacs))

}


if (.REDO_SIMS) {
    # Takes ~2 min w/ q = 2, 3 threads
    set.seed(807582316)
    one_d_sim <- tibble(.d2 = c(-10^(c(-2, -3, -4)), 0, 10^(c(-4, -3, -2))),
                        .max_t = 200e3L) %>%
        mutate(.max_t = ifelse(.d2 == -1e-4, 2e6L, .max_t)) %>%
        pmap(one_d_combo_vary_one)
    saveRDS(one_d_sim, rds("one_d_sim-q2"))
    one_d_sim_df <- map_dfr(one_d_sim, ~ .x[["NV"]])
} else {
    one_d_sim <- readRDS(rds("one_d_sim-q2"))
    one_d_sim_df <- map_dfr(one_d_sim, ~ .x[["NV"]])
}


#' #'
#' #' Code below shows that they're mostly stable, except
#' #' in some reps when d2 is 0 or 1e-4.
#' #'
#' map_dfr(1:length(one_d_sim),
#'         function(i) {
#'             eigs <- map_dbl(one_d_sim[[i]][["J"]],
#'                             function(.x){
#'                                 if (any(is.na(.x))) return(NA)
#'                                 max(Re(eigen(.x)$values))
#'                             })
#'             one_d_sim[[i]]$NV %>%
#'                 .[1,] %>%
#'                 select(d2) %>%
#'                 mutate(e_min = min(eigs),
#'                        e_max = max(eigs))
#'         }) %>%
#'     print_big_nums()
#'
#'
#' #'
#' #' If we only look at evolution (∂ V / ∂ V), then it's totally stable
#' #' with no complex eigenvalues:
#' #'
#' map_dfr(1:length(one_d_sim),
#'         function(i) {
#'             eigs <- map_dbl(one_d_sim[[i]][["J"]],
#'                             function(.x){
#'                                 nm <- nrow(.x) / (1 + 1 / 2)
#'                                 .x <- .x[1:nm, 1:nm]
#'                                 if (any(is.na(.x))) return(NA)
#'                                 max(eigen(.x)$values)
#'                             })
#'             one_d_sim[[i]]$NV %>%
#'                 .[1,] %>%
#'                 select(d2) %>%
#'                 mutate(e_min = min(eigs),
#'                        e_max = max(eigs))
#'         }) %>%
#'     print_big_nums()
#'
#'
#' #'
#' #' And when looking at only ecology (∂ N / ∂ N), we get some unstable
#' #' points at d2 = 0 and 1e-4, but no complex eigenvalues.
#' #'
#' map_dfr(1:length(one_d_sim),
#'         function(i) {
#'             eigs <- map_dbl(one_d_sim[[i]][["J"]],
#'                             function(.x){
#'                                 nqn <- nrow(.x)
#'                                 nm <- nqn / (1 + 1 / 2)
#'                                 .x <- .x[(nm+1):nqn, (nm+1):nqn]
#'                                 if (any(is.na(.x))) return(NA)
#'                                 max(eigen(.x)$values)
#'                             })
#'             one_d_sim[[i]]$NV %>%
#'                 .[1,] %>%
#'                 select(d2) %>%
#'                 mutate(e_min = min(eigs),
#'                        e_max = max(eigs))
#'         }) %>%
#'     print_big_nums()
#'
#'
#' #'
#' #' Every positive d2 produced complex eigenvalues:
#' #'
#' map_dfr(1:length(one_d_sim),
#'         function(i) {
#'             eigs <- map(one_d_sim[[i]][["J"]],
#'                         function(.x){
#'                             if (any(is.na(.x))) return(NA)
#'                             eigen(.x)$values
#'                         })
#'             cmplx <- map_lgl(eigs, ~ any(Im(.x) != 0))
#'             tibble(d2 = one_d_sim[[i]]$NV[["d2"]][[1]],
#'                    cmplx = mean(cmplx))
#'         })
#'
#'



# start complex-eigen simulations ----

#'
#' This part simulates a scenario that results in complex eigenvalues
#' at equilibrium.
#'
.d2 = 1e-4
.max_t = 500e3L

.d <- c(0.05, .d2)

#'
#' By setting the seed beforehand, we can look at reps at long time scales
#' to look at their equilibria, then do short-time-scale simulations to
#' plot the transient dynamics.
#' Trying to do both at the same time results in too-massive output.
#'

set.seed(634789)
Z_long <- quant_gen(q = length(.d), eta = etas[[2]], d = .d,
                    max_t = .max_t,
                    save_every = 0,
                    n_reps = 3, n = 100, N0 = rep(1, 100),
                    start_t = 0, sigma_V0 = 2, n_threads = .N_THREADS,
                    show_progress = TRUE)
set.seed(634789)
Z_short <- quant_gen(q = length(.d), eta = etas[[2]], d = .d,
                     max_t = 5000,
                     save_every = 1L,
                     n_reps = 3, n = 100, N0 = rep(1, 100),
                     start_t = 0, sigma_V0 = 2, n_threads = .N_THREADS,
                     show_progress = TRUE)


jacs <- jacobians(Z_long)

map(jacs, ~ which(Im(eigen(.x, only.values = TRUE)[["values"]]) != 0))


Z_short$nv %>%
    filter(rep == 2) %>%
    filter(time < 1000) %>%
    ggplot(aes(time, N)) +
    geom_line(aes(color = spp)) +
    scale_color_viridis_d(guide = FALSE)



# end complex-eigen simulations ----





lab_fun <- function(.x) {
    labs <- paste(.x)
    labs[.x != 0] <- sprintf("%s10^{%i}",
                             sign(.x[.x != 0]) %>% paste() %>% str_remove("1"),
                             abs(.x[.x != 0]) %>% log10())
    return(labs)
}

coexist_one_d_p <- one_d_sim_df %>%
    mutate(d2 = factor(d2, levels = d2 %>% sort() %>% unique(),
                       labels = d2 %>% sort() %>% unique() %>% lab_fun())) %>%
    ggplot(aes(d2, n_spp)) +
    geom_jitter(width = 0.2, alpha = 0.5, height = 0,
                color = "dodgerblue4", fill = "dodgerblue", shape = 21) +
    scale_color_viridis_d(guide = FALSE, end = 0.95, option = "D") +
    scale_x_discrete(expression("Value of" ~ d[2]), labels = rlang::parse_exprs) +
    scale_y_continuous("Number of surviving species", breaks = c(0, 15, 30),
                       limits = c(0, 30)) +
    theme(plot.margin = margin(0,0,0,t=12)) +
    NULL






# ------------------------------*
# __ Time series of d2 = -0.01 and d2 = -1e-4 ----
# ------------------------------*

# Same as above, but returns a time series
one_d_combo_vary_one_TS <- function(.d2, .max_t) {

    # .d2 = -1e-4; .max_t = 2e6
    # rm(.d2, .max_t, C, Z, .d)

    stopifnot(is.numeric(.d2) && length(.d2) == 1)

    .max_t <- as.integer(.max_t)

    .d <- c(0.05, .d2)

    Z <- quant_gen(q = 2, eta = etas[[2]], d = .d, max_t = .max_t, n_reps = 24,
                   save_every = .max_t %/% 100L, n = 100, N0 = rep(1, 100),
                   start_t = 0, sigma_V0 = 2, n_threads = .N_THREADS,
                   show_progress = TRUE)

    Z$call[["d"]] <- eval(.d)

    jacs <- jacobians(Z)

    return(list(NV = Z$nv %>%
                    filter(trait == 1) %>%
                    group_by(rep, time) %>%
                    summarize(N = n()) %>%
                    ungroup() %>%
                    mutate(d2 = .d2) %>%
                    identity(),
                J = jacs))

}



if (.REDO_SIMS) {
    # Takes 41 sec
    set.seed(1821003492)
    one_d_TS_sim <- tibble(.d = c(-1e-2, -1e-4), .max_t = c(10e3L, 2e6L)) %>%
        pmap(one_d_combo_vary_one_TS)
    saveRDS(one_d_TS_sim, rds("one_d_TS_sim"))
    one_d_TS_sim_df <- map_dfr(one_d_TS_sim, ~ .x[["NV"]]) %>%
        mutate(d2 = factor(d2, levels = c(-0.01, -1e-4)))
} else {
    one_d_TS_sim <- readRDS(rds("one_d_TS_sim"))
    one_d_TS_sim_df <- map_dfr(one_d_TS_sim, ~ .x[["NV"]]) %>%
        mutate(d2 = factor(d2, levels = c(-0.01, -1e-4)))
}


#' #'
#' #' Code below shows that they're all stable with no complex eigenvalues:
#' #'
#' map_dfr(1:length(one_d_TS_sim),
#'         function(i) {
#'             eigs <- map_dbl(one_d_TS_sim[[i]][["J"]],
#'                             function(.x){
#'                                 if (any(is.na(.x))) return(NA)
#'                                 max(eigen(.x)$values)
#'                             })
#'             one_d_TS_sim[[i]]$NV %>%
#'                 .[1,] %>%
#'                 select(d2) %>%
#'                 mutate(e_min = min(eigs),
#'                        e_max = max(eigs))
#'         }) %>%
#'     print_big_nums()




one_d_TS_p <- one_d_TS_sim_df %>%
    filter(time < 1.15e6) %>%
    mutate(time = time / 1000L) %>%
    ggplot(aes(time, N)) +
    geom_line(aes(group = rep), alpha = 0.5, size = 0.5, color = "dodgerblue") +
    geom_text(data = tibble(d2 = one_d_TS_sim_df$d2 %>% unique() %>% sort(),
                            time = c(5, 600),
                            N = 100,
                            lab = factor(d2, levels = levels(d2),
                                         labels = sprintf("d[2] == %s",
                                                          c("-10^{-2}",
                                                            "-10^{-4}")))),
              aes(label = lab), parse = TRUE, hjust = 0.5, vjust = 1) +
    facet_wrap(~ d2, nrow = 1, label = label_parsed, scales = "free_x") +
    scale_y_continuous("Number of species", breaks = c(0, 50, 100)) +
    scale_x_continuous(expression("Time (" %*% "1,000" * ")"),
                       breaks = scales::breaks_extended(n = 3)) +
    theme(strip.text = element_blank(),
          plot.margin = margin(0,0,0,t=12)) +
    NULL








# ------------------------------*
# __ combine them ----
# ------------------------------*



coexist_p <- ggarrange(coexist_all_d_p + theme(axis.title.y = element_blank()),
                       coexist_one_d_p + theme(axis.title.y = element_blank()),
                       one_d_TS_p + theme(axis.title.y = element_blank()),
                       heights = c(1, 1, 0.75),
                       left = "Number of surviving species",
                       ncol = 1, labels = LETTERS[1:3],
                       label.args = list(gp = gpar(fontsize = 14,
                                                   fontface =  "plain"),
                                         vjust = 1, hjust = -3))


if (.RESAVE_PLOTS) save_plot(coexist_p, 4, 5, "S1-")




# ===========================================================================*
# ===========================================================================*

# S2. Invasibility ----

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
#' For these simulations, q = 2, eta = -0.2
#'



one_d_invasion <- function(.x) {

    .d <- .x[[".d"]][[1]]
    .n <- .x[[".n"]][[1]]
    .inv_N0 <- .x[[".inv_N0"]][[1]]

    .eta <- -etas[[2]]

    V0 <- stable_points(.eta) %>%
        as.matrix() %>%
        split(1:nrow(.)) %>%
        rep(ceiling(.n / length(.))) %>%
        .[1:.n]
    N0 <- pop_sizes(.n, .eta, .d)

    inv_df <- tibble(.v1 = stable_points(.eta)[["V1"]][1],
                     .v2 = stable_points(.eta)[["V2"]][1] +
                         sort(c(0, -0.1 * 1:8, 0.1 * 1:8))) %>%
        pmap_dfr(function(.v1, .v2) {
            .inv <- quant_gen(q = 2, eta = .eta, d = .d, n = .n + 1,
                              save_every = 0L, max_t = 50e3L, n_reps = 1,
                              N0 = c(N0, N0[1] * .inv_N0),
                              V0 = c(V0, list(cbind(.v1, .v2))),
                              start_t = 0, sigma_V0 = 0,
                              show_progress = FALSE) %>%
                .[["nv"]] %>%
                filter(trait == 1, spp == .n+1) %>%
                nrow() %>%
                `>=`(1) %>%
                as.integer()
            tibble(dist = .v2 - stable_points(.eta)[["V2"]][1], invaded = .inv)
        }) %>%
        mutate(d = .d, res_n = .n, inv_N0 = .inv_N0)

    return(inv_df)

}


# Takes ~ 30 sec w/ 3 threads
invade_df <- crossing(.d = 0.01 * c(-1, 0, 1),
                      .n = c(2, 10),
                      .inv_N0 = c(0.01, 0.1, 1)) %>%
    split(1:nrow(.)) %>%
    mclapply(FUN = one_d_invasion, mc.cores = .N_THREADS)
invade_df <- bind_rows(invade_df)



invasion_caption <- "Successful invasion of an equilibrium community based on
                     the invader's starting distance (in trait space)
                     from the stable point.
                     Columns of sub-panels separate whether evolution was
                     conflicting, neutral, or non-conflicting.
                     Rows of sub-panels separate the invaders' starting
                     abundances ($N_{eq}$) in relation to the residents'
                     abundances ($N_{res}$); all residents had the same
                     abundance.
                     The point color indicates the number of species present
                     in the resident community."


invasion_p <- invade_df %>%
    mutate(invaded = factor(invaded, levels = 0:1, labels = c("no", "yes")),
           d = factor(d, levels = 0.01 * c(-1, 0, 1),
                      labels = c("'conflicting'", "'neutral'",
                                 "'non-conflicting'")),
           inv_N0 = factor(inv_N0, levels = c(0.01, 0.1, 1),
                           labels = sprintf("N[inv] == %s", c("frac(N[res],100)",
                                                              "frac(N[res],10)",
                                                              "N[res]"))),
           res_n = factor(res_n, levels = c(2, 10),
                          labels = sprintf("%i species", c(2, 10)))) %>%
    ggplot(aes(invaded, dist)) +
    geom_hline(yintercept = 0, linetype = 2, size = 0.5, color = "gray70") +
    geom_jitter(aes(color = res_n, fill = res_n),
                height = 0, width = 0.25, shape = 21) +
    facet_grid(inv_N0 ~ d, label = label_parsed) +
    xlab("Successful invasion") +
    ylab("Distance from stable point") +
    scale_color_viridis_d("starting community size:", begin = 0.3,
                          end = 0.7, option = "A") +
    scale_fill_viridis_d("starting community size:", begin = 0.3,
                         end = 0.7, option = "A", alpha = 0.5) +
    theme(panel.spacing = unit(2, "lines"),
          strip.text.y = element_text(angle = 0, margin = margin(0,0,0,l=6)),
          legend.position = "top",
          legend.title = element_text(size = 10)) +
    coord_flip() +
    NULL


if (.RESAVE_PLOTS) save_plot(invasion_p, 5, 4.5, "S2-")








# ===========================================================================*
# ===========================================================================*

# S3. "Filling in" of trait space ----

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


