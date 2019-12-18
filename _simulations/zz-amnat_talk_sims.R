

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
# Save a PNG file with nice arrows
save_plot_arrows <- function(plot_obj, .width, .height, .prefix = "", .suffix = "") {

    fn <- gsub("_p$", "", paste(substitute(plot_obj)))
    if (.prefix != "") fn <- paste0(.prefix, "_", fn)
    if (.suffix != "") fn <- paste0(fn, "_", .suffix)
    fn <- sprintf("_simulations/zz-talk_data/%s.png", fn)
    cat(fn, "\n")

    png(fn, width = .width, height = .height, units = "in", res = 300, bg = "black")

    grid.newpage()
    grid.draw(plot_obj)
    grid.force()
    # change shape of arrows
    grid.gedit("polyline", gp=gpar(linejoin = "mitre", lineend = "butt"))
    # change the shape in legend also
    grid.gedit("layout", gp=gpar(linejoin = "mitre", lineend = "butt"))

    dev.off()

    invisible(NULL)
}
color_scale <- scale_colour_viridis_d(guide = FALSE, option = "B", begin = 0.4)


# =======================*
# Simulations for time series ----
# =======================*

# N plays out over long-term, but traits change quickly, so running 2 versions of these
# to prevent making giant data frames:

# long_ts <- list(neg_d = NA, pos_d = NA)
# # Takes ~2 min
# set.seed(1)
# long_ts$neg_d <- quant_gen(q = 2, eta = -0.6, d = -1e-04, max_t = 5e6L, n_reps = 1,
#                      save_every = 1e4L, n = 100, N0 = rep(1e-3, 100),
#                      perturb_sd = 2)
# # Takes ~15 min
# set.seed(2)
# long_ts$pos_d <- quant_gen(q = 2, eta = -0.6, d = 1e-04, max_t = 10e6L, n_reps = 1,
#                            save_every = 1e4L, n = 100, N0 = rep(1e-3, 100),
#                            perturb_sd = 2)
# saveRDS(long_ts, "_simulations/zz-talk_data/long_ts.rds")

long_ts <- readRDS("_simulations/zz-talk_data/long_ts.rds")

short_ts <- list(neg_d = NA, pos_d = NA)
set.seed(1)
short_ts$neg_d <- quant_gen(q = 2, eta = -0.6, d = -1e-04, max_t = 2e3L, n_reps = 1,
                            save_every = 1L, n = 100, N0 = rep(1e-3, 100),
                            perturb_sd = 2)
set.seed(2)
short_ts$pos_d <- quant_gen(q = 2, eta = -0.6, d = 1e-04, max_t = 2e3L, n_reps = 1,
                            save_every = 1L, n = 100, N0 = rep(1e-3, 100),
                            perturb_sd = 2)




# ==============*
# Time series, d < 0 (conflicting) ----
# ==============*


# Species that won:
winner <- long_ts$neg_d$nv %>%
    filter(time == max(time)) %>%
    .[["spp"]] %>%
    paste() %>% as.integer() %>%
    unique()

# This shows that over long period, one species dominates
N_ts_p <- long_ts$neg_d$nv %>%
    filter(trait == 1, time < 1.4e6, N >= 1) %>%
    ggplot(aes(time, N)) +
    geom_line(aes(color = spp), alpha = 0.5, size = 0.5) +
    color_scale +
    theme_black() +
    scale_x_continuous("Generation",
                       breaks = seq(0, 1e6, 0.5e6),
                       labels = c(0, "500,000", "1,000,000")) +
    ylab("Abundance")

# Zoomed to look at other species
N_ts_zoom_p <- N_ts_p + coord_cartesian(ylim = c(0, 2000))

# save_plot(N_ts_p, .width = 6, .height = 4, .prefix = "-d")
# save_plot(N_ts_zoom_p, .width = 6, .height = 4, .prefix = "-d")


# Changes in trait values through time
V_ts_p <- short_ts$neg_d$nv %>%
    filter(time < 100) %>%
    rename(Trait = trait) %>%
    ggplot(aes(time, value)) +
    geom_line(aes(color = spp), alpha = 0.5) +
    facet_wrap(~ Trait, ncol = 1, labeller = function(x) label_both(x, sep = " ")) +
    xlab("Generation") +
    ylab("Trait value") +
    color_scale +
    theme_black()

# How they move through trait space (x's are starting points)
V_space_p <- short_ts$neg_d$nv %>%
    filter(time < 100) %>%
    mutate(trait = paste0("V", trait)) %>%
    spread(trait, value) %>%
    arrange(time) %>%
    ggplot(aes(V1, V2)) +
    stable_points(eta = eval(short_ts$neg_d$call[["eta"]]), return_geom = TRUE,
                  color = "gray40", shape = 19, size = 4) +
    # unstable_points(eta = eval(short_ts$neg_d$call[["eta"]]), return_geom = TRUE,
    #                 color = "gray40", shape = 1, size = 4) +
    geom_path(aes(color = spp), alpha = 0.5) +
    geom_point(data = filter(short_ts$neg_d$nv, time == min(time)) %>%
                   mutate(trait = paste0("V", trait)) %>%
                   spread(trait, value),
               aes(color = spp), alpha = 0.5, shape = 4, size = 1) +
    scale_x_continuous("Trait 1", limits = 4.533577 * c(-1,1)) +
    scale_y_continuous("Trait 2", limits = 4.533577 * c(-1,1)) +
    theme_black() +
    color_scale +
    coord_equal()


# And with the winner highlighted in white
V_space_winner_p <- V_space_p +
    geom_path(data = filter(short_ts$neg_d$nv, time < 100, spp == winner) %>%
                  mutate(trait = paste0("V", trait)) %>%
                  spread(trait, value) %>%
                  arrange(time),
              color = "white", size = 1) +
    geom_point(data = filter(short_ts$neg_d$nv, time == min(time), spp == winner) %>%
                   mutate(trait = paste0("V", trait)) %>%
                   spread(trait, value) %>%
                   arrange(time),
               color = "white", shape = 4, size = 2)


# save_plot(V_ts_p, .width = 4, .height = 6, .prefix = "-d")
# save_plot(V_space_p, .width = 5, .height = 5, .prefix = "-d")
# save_plot(V_space_winner_p, .width = 5, .height = 5, .prefix = "-d")







# ==============*
# Time series, d > 0 (non-conflicting) ----
# ==============*

# # Number of species that won:
# n_winners <- long_ts$pos_d$nv %>%
#     filter(time == max(time)) %>%
#     .[["spp"]] %>%
#     paste() %>% as.integer() %>%
#     unique() %>%
#     length()

# This shows that over long period, one species dominates
N_ts_p <- long_ts$pos_d$nv %>%
    filter(trait == 1) %>%
    ggplot(aes(time, N)) +
    geom_line(aes(color = spp), alpha = 0.5, size = 0.5) +
    color_scale +
    theme_black() +
    scale_x_continuous("Generation") +
    ylab("Abundance")

# save_plot(N_ts_p, .width = 6, .height = 4, .prefix = "+d")


# Changes in trait values through time
V_ts_p <- short_ts$pos_d$nv %>%
    filter(time < 100) %>%
    rename(Trait = trait) %>%
    ggplot(aes(time, value)) +
    geom_line(aes(color = spp), alpha = 0.5) +
    facet_wrap(~ Trait, ncol = 1, labeller = function(x) label_both(x, sep = " ")) +
    xlab("Generation") +
    ylab("Trait value") +
    color_scale +
    theme_black()

# How they move through trait space (x's are starting points)
V_space_p <- short_ts$pos_d$nv %>%
    filter(time < 100) %>%
    mutate(trait = paste0("V", trait)) %>%
    spread(trait, value) %>%
    arrange(time) %>%
    ggplot(aes(V1, V2)) +
    stable_points(eta = eval(short_ts$pos_d$call[["eta"]]), return_geom = TRUE,
                  color = "gray40", shape = 19, size = 4) +
    # unstable_points(eta = eval(short_ts$pos_d$call[["eta"]]), return_geom = TRUE,
    #                 color = "gray40", shape = 1, size = 4) +
    geom_path(aes(color = spp), alpha = 0.5) +
    geom_point(data = filter(short_ts$pos_d$nv, time == min(time)) %>%
                   mutate(trait = paste0("V", trait)) %>%
                   spread(trait, value),
               aes(color = spp), alpha = 0.5, shape = 4, size = 1) +
    scale_x_continuous("Trait 1", limits = 4.533577 * c(-1,1)) +
    scale_y_continuous("Trait 2", limits = 4.533577 * c(-1,1)) +
    theme_black() +
    color_scale +
    coord_equal()


# save_plot(V_ts_p, .width = 4, .height = 6, .prefix = "+d")
# save_plot(V_space_p, .width = 5, .height = 5, .prefix = "+d")




# ==============*
# Winner values, 2 traits ----
# ==============*

two_trait_pts_p <- tibble(eta = 0.6 * -1:1) %>%
    mutate(points = map(eta, stable_points, line_n = 10e3)) %>%
    unnest(cols = points) %>%
    mutate(eta = factor(eta)) %>%
    ggplot() +
    geom_point(aes(V1, V2, color = eta), size = 4) +
    scale_color_manual(expression(bold(eta)),
                       values = c("dodgerblue3", "goldenrod1", "gray70"),
                       guide = guide_legend(override.aes = list(alpha = 1, size = 4))) +
    xlab("Trait 1") +
    ylab("Trait 2") +
    theme_black() +
    theme(legend.key = element_blank(),
          legend.title = element_text(size = 20, margin = margin(0,0,0,b=-4)),
          plot.margin = margin(0,0,0,0)) +
    coord_equal()

# save_plot(two_trait_pts_p, .width = 4.5, .height = 3.75)






# ==============*
# Trajectories, 2 traits ----
# ==============*


# Omega at stable points
O_eq <- function(eta, f = 0.1, a0 = 0.5, r0 = 0.5) {
    rho <- ifelse(eta < 0, f * (1 + eta), f * (1 - eta))
    (rho / a0) * exp((r0 / rho) - 1)
}
# Change in x,y coordinates given angle and length
delta_xy <- function(len, angle) {
    cbind(len * cos(angle), len * sin(angle))
}


# Trajectories for one eta and Omega proportion, using analytical solutions
eta_Omega_traj <- function(eta, O_prop) {

    # eta = -0.6; O_prop = 0
    # rm(eta, O_prop, .max_t, Omega, pts, radius, xy_df, one_run, sim_df)

    Omega <- O_eq(eta) * O_prop

    pts <- stable_points(eta)
    radius <- sqrt(pts[[1]][1]^2 + pts[[1]][2]^2)
    C <- matrix(eta, 2, 2)
    diag(C) <- 1

    xy_mat <- crossing(dists = seq(0, radius * 1.25, length.out = 5)[-1],
                       angles = 2 * pi * seq(0, 1, length.out = 16 + 1)[-1]) %>%
        pmap_dfr(.f = ~ tibble(x = .x * cos(.y), y = .x * sin(.y))) %>%
        as.matrix()

    f <- formals(quant_gen)[["f"]]
    a0 <- formals(quant_gen)[["a0"]]
    r0 <- formals(quant_gen)[["r0"]]
    n <- 1
    add_var <- eval(formals(quant_gen)[["add_var"]])

    one_run <- function(i) {
        out <- tibble(eta = eta,
                      O_prop = O_prop,
                      id = sprintf("%.2g_%.2g_%i", eta, O_prop, i),
                      V = list(NA),
                      angle = NA_real_,
                      magnitude = NA_real_)
        delta_V <- add_var * sauron:::sel_str_cpp(V = list(xy_mat[i,,drop=F]),
                                                  N = Omega, f, a0, C, r0, d = 0)
        out$V[[1]] <- cbind(V1 = xy_mat[i,][["x"]], V2 = xy_mat[i,][["y"]])
        out$angle <- atan2(delta_V[[2]], delta_V[[1]])
        out$magnitude <- sqrt(sum(delta_V^2))
        return(out)
    }

    sim_df <- map_dfr(1:nrow(xy_mat), one_run)

    return(sim_df)

}

fitness <- function(V1, V2, eta, O_prop, f = 0.1, a0 = 0.5, r0 = 0.5) {
    Omega <- O_eq(eta) * O_prop
    C <- matrix(eta, 2, 2)
    diag(C) <- 1
    Vi <- cbind(V1, V2)
    exp(r0 - f * Vi %*% C %*% t(Vi) - a0 * exp(- Vi %*% t(Vi)) * Omega)
}


fitness_landscape <- function(eta, O_prop, max_mag = 4) {
    hm_df <- crossing(O_prop = O_prop,
                      eta = eta,
                      V1 = seq(-max_mag, max_mag, length.out = 100),
                      V2 = seq(-max_mag, max_mag, length.out = 100)) %>%
        mutate(fit = pmap_dbl(list(V1, V2, eta, O_prop),
                              ~ fitness(..1, ..2, eta = ..3, O_prop = ..4)))

    return(hm_df)
}



alt_states <- bind_rows(map_dfr(as.numeric(levels(traj_df$eta)),
                                ~ stable_points(eta = .x) %>%
                                    mutate(eta = .x, type = "stable")),
                        map_dfr(as.numeric(levels(traj_df$eta)),
                                ~ unstable_points(eta = .x) %>%
                                    mutate(eta = .x, type = "unstable")))

fl_df <- crossing(eta = 0.6 * -1:1,
                  O_prop = c(0, c(10^(c(-3, -1, 0))))) %>%
    pmap_dfr(fitness_landscape) %>%
    mutate(O_prop = factor(O_prop, levels = rev(sort(unique(O_prop)))))


traj_df <- crossing(eta = 0.6 * -1:1,
                    O_prop = c(0, c(10^(c(-3, -1, 0))))) %>%
    pmap_dfr(eta_Omega_traj) %>%
    mutate_at(vars(eta, id), factor) %>%
    mutate(O_prop = factor(O_prop, levels = rev(sort(unique(O_prop)))))


traj_df2 <- traj_df %>%
    mutate(V = pmap(list(V, angle, magnitude), function(.V, .a, .m) {
        d_V <- delta_xy(len = 0.3, angle = .a)
        Vt1 <- .V + d_V
        M <- as_tibble(cbind(time = 0:1, rbind(.V, Vt1)))
        return(M)
    })) %>%
    unnest(cols = V)



# for (eta_ in as.numeric(levels(traj_df$eta))) {
#
#     if (eta_ == 0) {
#         alt_p <- geom_path(data = filter(alt_states, eta == eta_), aes(V1, V2),
#                            size = 1, color = "white", linetype = 1)
#         size_scale <- scale_size("magnitude",
#                                  range = c(0.1, 2), trans = "log10", breaks = 10^c(-1,0),
#                                  labels = map(c(-1,0), ~bquote(10^{.(.x)})))
#     } else {
#         alt_p <- geom_point(data = filter(alt_states, eta == eta_),
#                             aes(V1, V2, shape = type), size = 4, color = "white")
#         size_scale <- scale_size("magnitude",
#                                  range = c(0.1, 2), trans = "log10", breaks = 10^c(-1,1,3),
#                                  labels = map(c(-1,1,3), ~bquote(10^{.(.x)})))
#     }
#
#     fit_land_trait_traj_p <- fl_df %>%
#         filter(eta == eta_) %>%
#         ggplot() +
#         geom_raster(aes(V1, V2, fill = fit), interpolate = FALSE) +
#         geom_path(data = filter(traj_df2, eta == eta_),
#                   aes(V1, V2, group = id, size = magnitude), color = "gray50", alpha = 1,
#                   arrow = arrow(length = unit(2, "points"), type = "closed", angle = 20)) +
#         alt_p +
#         facet_wrap(~ O_prop, nrow = 2) +
#         scale_fill_gradient2("fitness", high = "firebrick", low = "dodgerblue3",
#                              mid = "black", midpoint = 1) +
#         scale_shape_manual(values = c(19, 1), guide = FALSE) +
#         size_scale +
#         xlab("Trait 1") +
#         ylab("Trait 2") +
#         theme_black() +
#         theme(plot.margin = margin(0,0,0,0)) +
#         coord_equal()
#
#     # Filename suffix based on eta_
#     suf_ <- paste0(ifelse(eta_ > 0, "+", ifelse(eta_ < 0, "-", "0")), "eta")
#     save_plot_arrows(fit_land_trait_traj_p, .width = 6, .height = 6, .suffix = suf_)
#
# }; rm(eta_)





# ==============*
# Coexistence ----
# ==============*




eta_ <- -0.6
d_ <- 0

invasion <- function(eta_, d_, p_, max_t_ = 250e3) {
    quant_gen(q = 2, eta = eta_, d = d_, max_t = max_t_, n_reps = 1,
              save_every = 1, n = 2, perturb_sd = 0,
              N0 = c(O_eq(eta_) * p_, O_eq(eta_)),
              V0 = list(as.matrix(stable_points(eta_)[1,]),
                        as.matrix(stable_points(eta_)[1,])),
              show_progress = FALSE) %>%
        .[["nv"]] %>%
        filter(trait == 1) %>%
        mutate(spp = factor(as.integer(spp), levels = 1:2,
                            labels = c("invader", "resident")),
               prop = p_, eta = eta_, d = d_) %>%
        select(eta, d, prop, time, spp, N) %>%
        add_row(eta = eta_, d = d_, prop = p_, time = 0, spp = c("invader", "resident"),
                N = c(O_eq(eta_) * p_, O_eq(eta_))) %>%
        arrange(time, spp)
}


# Takes ~30 sec
inv_df <- crossing(eta_ = -0.6,
                   d_ = 1e-4 * -1:1,
                   p_ = c(0.5, 1, 1.5)) %>%
    pmap_dfr(invasion) %>%
    mutate_at(vars(d, prop), factor) %>%
    select(-eta)

inv_df2 <- inv_df %>%
    # filter(spp == "invader") %>%
    filter(time <= 100e3, time %% 100 == 0) %>%
    mutate(d = factor(as.numeric(paste(d)), levels = 1e-4 * -1:1,
                      labels = c("conflicting", "neutral", "non-conflicting")))

invasion_p <- inv_df2 %>%
    ggplot(aes(time, N)) +
    # geom_line(aes(color = prop), size = 1) +
    geom_line(aes(color = spp), size = 0.75, alpha = 1) +
    geom_point(data = inv_df2 %>% filter(time == 0), aes(color = spp), size = 2) +
    facet_grid(prop ~ d) +
    scale_colour_viridis_d(guide = FALSE, option = "A", begin = 0.4) +
    scale_linetype_manual(values = c(2, 3)) +
    scale_x_continuous("Generation", breaks = c(0, 40e3, 80e3),
                       labels = c("0", "40e3", "80e3")) +
    scale_y_continuous("Abundance", limits = c(0, NA),
                       breaks = c(0, 5e3, 10e3),
                       labels = c("0", "5e3", "10e3")) +
    theme_black() +
    theme(strip.text.x = element_text(size = 12),
          strip.text.y = element_blank())


save_plot(invasion_p, .width = 5, .height = 5)



# ### Three traits
#
# q <- 3
# i <- 1
# for (e in -1:1) {
#     for (d in -1:1) {
#         args <- quant_gen_args(e, d, q)
#         args$n_reps <- 8
#         cat(sprintf("\n\neta: %.2g || d: %.2g (%i / 9) \n", args$eta, args$d, i))
#         # Using already-generated seed:
#         seed <- filter(qg_seeds, n_traits == q, eta_par == e, d_par == d) %>%
#             .[["seed"]]
#         set.seed(seed)
#         qg <- do.call(quant_gen, args)
#         fn <- sprintf("_simulations/zz-talk_data/quant_gen/q%i/eta=%s__d=%s.rds", q, e, d)
#         saveRDS(qg, fn, compress = FALSE)
#         i <- i + 1
#     }
# }





qg <- quant_gen(q = 2, eta = -0.6, d = 1e-04, max_t = 10e3L, n_reps = 4,
                save_every = 1L, n = 2, perturb_sd = 0, n_cores = 4,
                V0 = list(as.matrix(stable_points(-0.6)[1,]),
                          as.matrix(stable_points(-0.6)[1,])))


qg <- quant_gen(q = 2, eta = -0.6, d = 0, max_t = 2e6L, n_reps = 4,
                save_every = 1e3L, n = 1, perturb_sd = 0, n_cores = 4,
                V0 = list(as.matrix(stable_points(-0.6)[1,])))


eta <- -0.6
d <- 1e-4

N <- function(eta, d, f = 0.1, a0 = 0.5, r0 = 0.5) {
    rho <- f * (1 + eta)
    (rho / a0) * exp((r0 / rho) - 1)
}


N(-0.6, 0)



Omega_hat <- function(eta, f = 0.1, a0 = 0.5, r0 = 0.5) {
    rho <- ifelse(eta < 0, f * (1 + eta), f * (1 - eta))
    (rho / a0) * exp((r0 / rho) - 1)
}
# St


Omega_hat(-0.6)

N <- 3950.901
Vi <- as.matrix(stable_points(-0.6)[1,])
N + N * exp(-d * Vi[1]^2 * Vi[2]^2)

qg %>%
    .[["nv"]] %>%
    filter(time == max(time)) %>%
    as.data.frame()


qg %>%
    .[["nv"]] %>%
    filter(time < 250) %>%
    mutate(id = interaction(rep, spp)) %>%
    ggplot(aes(time, value)) +
    geom_line(aes(group = id)) +
    facet_wrap(~ trait + spp)


qg %>%
    .[["nv"]] %>%
    filter(time < 100, rep == 1) %>%
    spread()
ggplot(aes(time, value)) +
    geom_path(aes(group = spp)) +
    facet_wrap(~ trait + spp)


qg %>%
    .[["nv"]] %>%
    filter(time == max(time) & N > 7500) %>%
    as.data.frame()


qg %>%
    .[["nv"]] %>%
    filter(time == min(time))


qg %>%
    .[["nv"]] %>%
    filter(trait == 1) %>%
    mutate(id = interaction(rep, spp)) %>%
    ggplot(aes(time, N)) +
    geom_line(aes(group = id, color = spp), alpha = 0.5)

