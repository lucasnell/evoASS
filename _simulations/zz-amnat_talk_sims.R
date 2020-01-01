

suppressPackageStartupMessages({
    library(tidyverse)
    library(sauron)
    library(grid)
})

source(".Rprofile")


save_plot <- function(plot_obj, .width, .height, .prefix = "", .suffix = "") {
    fn <- gsub("_p$", "", paste(substitute(plot_obj)))
    if (.prefix != "") fn <- paste0(.prefix, "_", fn)
    if (.suffix != "") fn <- paste0(fn, "_", .suffix)
    fn <- sprintf("_simulations/zz-talk_data/%s.png", fn)
    cat(fn, "\n")
    ggsave(fn, plot_obj, width = .width, height = .height, dpi = 300, bg = "black")
    invisible(NULL)
}
color_scale <- scale_colour_viridis_d(guide = FALSE, option = "B", begin = 0.4)

unq_spp_filter <- function(..., .prec = 0.001) {
    V <- list(...)
    V <- do.call(cbind, V)
    V <- split(V, row(V))
    return(as.logical(sauron:::unq_spp_cpp(V, precision = .prec)))
}

# Making points look consistent for showing alternative states:
state_geom <- function(.dd, ...) {
    list(geom_point(data = .dd, shape = 19, alpha = 0.5, ...),
         geom_point(data = .dd, shape = 1, ...))
}


# # Looking at palettes:
# library(scales)
# show_col(viridis_pal(option = "B", begin = 0.3)(6))








# ====================================================================================*
# ====================================================================================*

# Time series ----

# ====================================================================================*
# ====================================================================================*


# ==============*
# * d = 0 (neutral) ----
# ==============*

set.seed(3)
trait_ts <- quant_gen(q = 2, eta = 0.6, d = 0, max_t = 20e3L, n_reps = 1,
                      save_every = 1L, n = 100, N0 = rep(1e-3, 100),
                      perturb_sd = 2)

# Changes in trait values through time
V_ts_p <- trait_ts$nv %>%
    filter(time < 75) %>%
    rename(Trait = trait) %>%
    ggplot(aes(time, value)) +
    geom_hline(yintercept = 0, linetype = 1, color = "gray30", size = 1) +
    geom_vline(xintercept = 0, linetype = 1, color = "gray30", size = 1) +
    geom_line(aes(color = spp), alpha = 0.5) +
    facet_wrap(~ Trait, ncol = 1, labeller = function(x) label_both(x, sep = " ")) +
    xlab("Generation") +
    ylab("Trait value") +
    color_scale +
    theme_black() +
    theme(panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          plot.margin = margin(0,0,0,0),
          strip.text = element_blank()) +
    coord_cartesian(ylim = c(0, 4.25))




# How they move through trait space (x's are starting points)
V_space_p <- trait_ts$nv %>%
    filter(time < 100) %>%
    mutate(trait = paste0("V", trait)) %>%
    spread(trait, value) %>%
    arrange(time) %>%
    ggplot(aes(V1, V2)) +
    geom_hline(yintercept = 0, linetype = 1, color = "gray30", size = 1) +
    geom_vline(xintercept = 0, linetype = 1, color = "gray30", size = 1) +
    geom_point(data = filter(trait_ts$nv, time == max(time)) %>%
                   mutate(trait = paste0("V", trait)) %>%
                   spread(trait, value) %>%
                   filter(unq_spp_filter(V1, V2)),
               color = "gray60", shape = 19, size = 4) +
    # stable_points(eta = eval(trait_ts$call[["eta"]]), return_geom = TRUE,
    #               color = "gray40", shape = 19, size = 4) +
    # unstable_points(eta = eval(trait_ts$call[["eta"]]), return_geom = TRUE,
    #                 color = "gray40", shape = 1, size = 4) +
    geom_path(aes(color = spp), alpha = 0.5) +
    geom_point(data = filter(trait_ts$nv, time == min(time)) %>%
                   mutate(trait = paste0("V", trait)) %>%
                   spread(trait, value),
               aes(color = spp), alpha = 0.5, shape = 4, size = 1) +
    xlab("Trait 1") +
    ylab("Trait 2") +
    theme_black() +
    theme(panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          plot.margin = margin(0,0,0,0)) +
    color_scale +
    coord_equal(ylim = c(0, 4.5), xlim = c(0, 4.5))


# save_plot(V_ts_p, .width = 3, .height = 6)
# save_plot(V_space_p, .width = 4, .height = 4)



# ====================================================================================*
# ====================================================================================*

# Winner values ----

# ====================================================================================*
# ====================================================================================*

# ==============*
# * 2 traits ----
# ==============*


# Simulations from this scenario doesn't coincide with calculated stable points,
# because of constraining all V >= 0
pts_q2_pe <- quant_gen(q = 2, eta = 0.6, d = 0, max_t = 2e3L, n_reps = 12,
                       save_every = 0L, n = 100, N0 = rep(1e-3, 100),
                       perturb_sd = 2, n_cores = 4, show_progress = FALSE) %>%
    .[["nv"]] %>%
    mutate(trait = factor(paste0("V", paste(trait)))) %>%
    spread("trait", "value") %>%
    filter(unq_spp_filter(V1, V2))


pts_q2_p <- ggplot(NULL, aes(V1, V2)) +
    geom_hline(yintercept = 0, linetype = 1, color = "gray30", size = 1) +
    geom_vline(xintercept = 0, linetype = 1, color = "gray30", size = 1) +
    state_geom(.dd = stable_points(-0.6), color = "firebrick3", size = 8) +
    state_geom(.dd = pts_q2_pe, color = "dodgerblue", size = 8) +
    xlab("Trait 1") +
    ylab("Trait 2") +
    theme_black() +
    theme(panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.title = element_text(size = 20, margin = margin(0,0,0,b=-4)),
          plot.margin = margin(0,0,0,0)) +
    scale_x_continuous(breaks = 1:2) +
    scale_y_continuous(breaks = 1:2) +
    coord_equal(xlim = c(0, 2.5), ylim = c(0, 2.5))




# save_plot(pts_q2_p, .width = 4.5, .height = 3.75)



# ==============*
# * 3 traits ----
# ==============*



set.seed(7890524)
pts_q3 <- map_dfr(c(-0.1, 0.1),
                  ~ quant_gen(q = 3, eta = .x, d = 0, max_t = 200e3L, n_reps = 12,
                              save_every = 0L, n = 100, N0 = rep(1e-3, 100),
                              perturb_sd = 2, n_cores = 4) %>%
                      .[["nv"]] %>%
                      mutate(eta = .x)) %>%
    mutate(trait = factor(paste0("V", paste(trait)))) %>%
    spread("trait", "value") %>%
    mutate(eta = factor(eta))


pts_q3_p <- ggplot(NULL, aes(V1, V2, size = V3)) +
    geom_hline(yintercept = 0, linetype = 1, color = "gray30", size = 1) +
    geom_vline(xintercept = 0, linetype = 1, color = "gray30", size = 1) +
    state_geom(pts_q3 %>%
                   group_by(eta) %>%
                   filter(unq_spp_filter(V1, V2, V3, .prec = 0.1)) %>%
                   arrange(V3), mapping = aes(color = eta)) +
    theme_black() +
    theme(panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.text = element_blank(),
          plot.margin = margin(0,0,0,0)) +
    scale_size_continuous(NULL, range = c(2, 10), breaks = 0.5 * 1:3) +
    guides(size = guide_legend(override.aes = list(color = "white"))) +
    scale_color_manual(values = c("firebrick3", "dodgerblue"),
                       drop = FALSE, guide = FALSE) +
    coord_equal(xlim = c(0, 2.5), ylim = c(0, 2.5)) +
    NULL



# save_plot(pts_q3_p, .width = 4.25, .height = 3.75)




# --------------*
# ... varying eta ----
# --------------*


vary_eta <- function(e1, e2, e3, ...) {
    C <- matrix(1, 3, 3)
    C[lower.tri(C)] <- c(e1, e2, e3)
    C[upper.tri(C)] <- C[lower.tri(C)]
    args <- list(q = 3, eta = C,
                 d = 0, max_t = 200e3L, n_reps = 12,
                 save_every = 0, n = 100, N0 = rep(1e-3, 100),
                 perturb_sd = 2, n_cores = 4)
    new_args <- list(...)
    for (n in names(new_args)) args[[n]] <- new_args[[n]]
    qg <- do.call(quant_gen, args)
    qg_nv <- qg %>%
        .[["nv"]] %>%
        mutate(eta1 = e1, eta2 = e2, eta3 = e3)
    return(qg_nv)
}

set.seed(1105998337)
pts_var_eta_q3 <- crossing(e1 = c(-1, 1) * 0.1,
                           e2 = c(-1, 1) * 0.1,
                           e3 = c(-1, 1) * 0.1) %>%
    pmap_dfr(vary_eta) %>%
    mutate(id = paste(paste(spp), paste(rep), sep = "_"))


pts_var_eta_q3 %>%
    mutate(trait = factor(paste0("V", paste(trait)))) %>%
    spread("trait", "value") %>%
    arrange(V3) %>%
    group_by(eta1, eta2, eta3) %>%
    summarize(n_states = sum(unq_spp_filter(V1, V2, V3))) %>%
    as.data.frame() %>%
    print(row.names = FALSE)


pts_var_eta_q3 %>%
    mutate(trait = factor(paste0("V", paste(trait)))) %>%
    spread("trait", "value") %>%
    arrange(V3) %>%
    group_by(eta1, eta2, eta3) %>%
    filter(unq_spp_filter(V1, V2, V3)) %>%
    ungroup() %>%
    ggplot(aes(V1, V2, size = V3)) +
    geom_hline(yintercept = 0, linetype = 1, color = "gray30", size = 1) +
    geom_vline(xintercept = 0, linetype = 1, color = "gray30", size = 1) +
    geom_point(alpha = 0.5, shape = 19, color = "dodgerblue") +
    geom_point(shape = 1, color = "dodgerblue") +
    facet_wrap(~ eta1 + eta2 + eta3, nrow = 2) +
    theme_black() +
    theme(panel.border = element_blank(),
          # axis.ticks = element_blank(),
          # axis.title = element_blank(),
          # axis.text = element_blank(),
          plot.margin = margin(0,0,0,0)) +
    scale_size_continuous("Trait 3", breaks = c(0.5, 1, 1.5), range = c(2, 8)) +
    # coord_equal(xlim = c(0, 2.5), ylim = c(0, 2.5)) +
    coord_equal() +
    NULL











# ====================================================================================*
# ====================================================================================*

# Varying d ----

# ====================================================================================*
# ====================================================================================*



#'
#' It appears that when `d` varies among traits, the trait with the negative
#' `d` (conflicting evolution) eventually gets maximized while the others go to ~0.
#'


one_d <- function(.d, ...) {
    args <- list(q = 2, eta = 0.6, d = .d, max_t = 200e3L, n_reps = 12,
                 save_every = 0L, n = 100, N0 = rep(1e-3, 100),
                 perturb_sd = 2, n_cores = 4, show_progress = FALSE)
    new_args <- list(...)
    for (n in names(new_args)) args[[n]] <- new_args[[n]]
    dd <- do.call(quant_gen, args)
    dd$nv %>%
        mutate(trait = factor(paste0("V", paste(trait))),
               d = paste(sign(eval(dd$call[["d"]])), collapse = "_")) %>%
        spread("trait", "value")
}


var_d <- map_dfr(list(c(0, 1e-2),
                      c(0, -1e-2),
                      c(1e-2, -1e-2)),
                 ~ one_d(.d = .x))



# npp <- one_d(.d = c(-1e-2, 1e-2, 1e-2), q = 3, max_t = 2e6L)
#
# nppts <- quant_gen(q = 3, eta = 0.6, d = c(-1e-2, 1e-2, 1e-2), max_t = 2e3L, n_reps = 12,
#                    save_every = 1L, n = 100, N0 = rep(1e-3, 100),
#                    perturb_sd = 2, n_cores = 4, show_progress = FALSE)
#
# nppts %>%
#     .[["nv"]] %>%
#     # filter(time < 100, time > 75) %>%
#     filter(time > (max(time) - 10)) %>%
#     mutate(id = factor(paste(rep, spp, sep = "_"))) %>%
#     ggplot(aes(time, value)) +
#     geom_line(aes(color = id)) +
#     facet_wrap(~ trait) +
#     color_scale +
#     theme_black() +
#     coord_cartesian(ylim = c(0, 0.2))
#
# nppts %>%
#     .[["nv"]] %>%
#     filter(time > (max(time) - 50)) %>%
#     mutate(trait = paste0("V", paste(trait))) %>%
#     spread("trait", "value") %>%
#     filter(V3 > 1.8) %>%
#     select(-V3) %>%
#     gather("trait", "value", V1:V2) %>%
#     mutate(id = factor(paste(rep, spp, sep = "_"))) %>%
#     filter(id == sample(levels(id), 1)) %>%
#     ggplot(aes(time, value)) +
#     geom_line(aes(color = id)) +
#     facet_wrap(~ trait, nrow = 2) +
#     color_scale +
#     theme_black() +
#     coord_cartesian(ylim = c(0, 0.2))
#
#
#
#
# npp %>%
#     filter(unq_spp_filter(V1, V2)) %>%
#     ggplot(aes(V1, V2)) +
#     geom_point(color = "white", size = 3) +
#     theme_black() +
#     xlim(0, 2.5) +
#     ylim(0, 2.5) +
#     coord_equal()


var_d %>%
    filter(d == "0_-1") %>%
    filter(unq_spp_filter(V1, V2)) %>%
    ggplot(aes(V1, V2)) +
    geom_point(color = "white", size = 3) +
    theme_black()


var_d %>%
    group_by(d, rep) %>%
    summarize(n_spp = n()) %>%
    ungroup() %>%
    ggplot(aes(d, n_spp)) +
    geom_point()



#'
#' If `d` for 2 traits is `+-`, then if >1 species start above the main diagonal,
#' one should
#'
#'
test_exclusion <- function(.i, .exclude, .n = 10, .seed = NULL, ...) {
    cat(.i, "")
    if (!is.null(.seed)) set.seed(.seed)
    V0 <- map(1:.n,
              function(x) {
                  a <- abs(rnorm(1, 0, 2))
                  b <- abs(rnorm(1, 0, 2))
                  c <- min(c(a, b)) * runif(1)
                  return(cbind(V1 = a, V2 = b, V3 = c))
              })
    if (.exclude) {
        inds <- sample.int(length(V0), length(V0) / 2)
        for (i in inds) V0[[i]] <- rev(V0[[i]])
    }
    args <- list(q = 3, eta = 0.6, d = c(1e-2, 1e-2, -1e-2), max_t = 200e3L, n_reps = 1,
                 save_every = 100L, n = .n, N0 = rep(1e-3, .n), V0 = V0,
                 perturb_sd = 0, show_progress = FALSE)
    new_args <- list(...)
    for (n in names(new_args)) args[[n]] <- new_args[[n]]
    qg_nv <- do.call(quant_gen, args) %>%
        .[["nv"]] %>%
        mutate(exclude = .exclude,
               rep = .i) %>%
        select(exclude, rep, everything())
    return(qg_nv)
}


set.seed(790824507)
seeds <- sample.int(2^31-1, 20)
excl_combos <- crossing(.i = 1:10, .exclude = c(TRUE, FALSE)) %>%
    mutate(.seed = seeds)
excl_df <- excl_combos %>%
    pmap_dfr(test_exclusion) %>%
    mutate(rep = factor(rep), exclude = factor(exclude))

excl_df %>%
    filter(time == max(time)) %>%
    group_by(exclude, rep) %>%
    summarize(n_spp = length(unique(spp)))

excl_df %>%
    filter(time == max(time)) %>%
    mutate(trait = paste0("V", trait)) %>%
    spread(trait, value) %>%
    group_by(exclude, rep) %>%
    summarize(n_states = sum(unq_spp_filter(V1, V2, V3)))


excl_df %>%
    filter(trait == 1, time < 14e3) %>%
    filter(rep == 1) %>%
    mutate(id = factor(paste(exclude, rep, spp, sep = "_"))) %>%
    ggplot(aes(time, N)) +
    geom_hline(yintercept = 0, linetype = 1, color = "gray30", size = 1) +
    geom_vline(xintercept = 0, linetype = 1, color = "gray30", size = 1) +
    # geom_line(aes(color = rep, group = id)) +
    geom_line(aes(color = exclude, group = id)) +
    facet_wrap(~ exclude, nrow = 2, scales = "free_y") +
    theme_black() +
    theme(strip.text = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          plot.margin = margin(0,0,0,0)) +
    scale_color_manual(values = c("#0F9E52", "#7E2EB3"), guide = FALSE)



# Shorter simulations for rep 1 with finer temporal resolution

excl_short_df <- excl_combos %>%
    filter(.i == 1) %>%
    pmap_dfr(test_exclusion, max_t = 2e3L, save_every = 1L) %>%
    mutate(rep = factor(rep), exclude = factor(exclude))

excl_short_df %>%
    mutate(trait = paste0("V", trait)) %>%
    spread(trait, value) %>%
    arrange(time) %>%
    mutate(id = factor(paste(exclude, rep, spp, sep = "_"))) %>%
    ggplot(aes(V1, V2)) +
    geom_hline(yintercept = 0, linetype = 1, color = "gray30", size = 1) +
    geom_vline(xintercept = 0, linetype = 1, color = "gray30", size = 1) +
    geom_point(data = filter(excl_short_df, time == max(time)) %>%
                   mutate(trait = paste0("V", trait)) %>%
                   spread(trait, value) %>%
                   filter(unq_spp_filter(V1, V2)),
               color = "gray60", shape = 19, size = 4) +
    geom_path(aes(color = exclude, group = id), alpha = 0.5) +
    geom_point(data = filter(excl_short_df, time == min(time)) %>%
                   mutate(trait = paste0("V", trait)) %>%
                   spread(trait, value),
               aes(color = exclude), alpha = 0.5, shape = 4, size = 1) +
    xlab("Trait 1") +
    ylab("Trait 2") +
    theme_black() +
    facet_wrap(~ exclude) +
    theme(panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          plot.margin = margin(0,0,0,0)) +
    scale_color_manual(values = c("#0F9E52", "#7E2EB3"), guide = FALSE) +
    coord_equal()


