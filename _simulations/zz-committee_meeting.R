
suppressPackageStartupMessages({
    library(tidyverse)
    library(sauron)
    library(grid)
    library(cowplot)
})

source(".Rprofile")



save_plot <- function(plot_obj, .width, .height, .prefix = NULL, .suffix = NULL) {
    fn <- gsub("_p$", "", paste(substitute(plot_obj)))
    if (!is.null(.prefix)) fn <- paste0(.prefix, fn)
    if (!is.null(.suffix)) fn <- paste0(fn, .suffix)
    fn <- sprintf("_simulations/zz-comm_plots/%s.pdf", fn)
    cat(fn, "\n")
    cairo_pdf(fn, width = .width, height = .height)
    plot(plot_obj)
    dev.off()
    invisible(NULL)
}

unq_spp_filter <- function(..., .prec = 0.001) {
    V <- list(...)
    V <- do.call(cbind, V)
    V <- split(V, row(V))
    return(as.logical(sauron:::unq_spp_cpp(V, precision = .prec)))
}


.q <- 3



# =======================================================================================`
# =======================================================================================`

#               Trait outcomes ----

# =======================================================================================`
# =======================================================================================`

set.seed(1652339501)
etas <- runif(.q * ((.q - 1) / 2), 0.1, 0.4)



one_eta_combo <- function(sign1, sign2, sign3) {

    eta <- matrix(0, .q, .q)
    eta[lower.tri(eta)] <- abs(etas) * c(sign1, sign2, sign3)
    # sample(c(-1, 1), .q, TRUE)
    eta <- eta + t(eta)
    diag(eta) <- 1

    trait_ts <- quant_gen(q = .q, eta = eta, d = 0, max_t = 20e3L, n_reps = 24,
                          save_every = 0L, n = 100, N0 = rep(1, 100),
                          start_t = 0, perturb_sd = 2, n_cores = 3,
                          show_progress = FALSE)

    trait_ts$nv %>%
        mutate(trait = paste0("V", trait)) %>%
        spread(trait, value) %>%
        filter(unq_spp_filter(V1, V2, V3)) %>%
        select(starts_with("V")) %>%
        mutate(eta1 = eta[2,1], eta2 = eta[3,1], eta3 = eta[3,2])
}



set.seed(145746935)
eta_sim_df <- crossing(sign1 = c(-1,1), sign2 = c(-1,1), sign3 = c(-1,1)) %>%
    pmap_dfr(one_eta_combo)

eta_sim_df






trait_outcomes_p <- eta_sim_df %>%
    mutate_at(vars(starts_with("eta")),
              ~ factor(sign(.x), levels = c(-1,1), labels = c("-", "+"))) %>%
    mutate(id = interaction(eta1, eta2, eta3, sep = "{}"),
           id = factor(id, levels = levels(id),
                       labels = sprintf("%s{}", levels(id)))) %>%
    ggplot(aes(V1, V2, size = V3)) +
    geom_point(shape = 1, color = "dodgerblue3") +
    geom_point(shape = 19, color = "dodgerblue3", alpha = 0.5) +
    scale_size_continuous("Trait 3", range = c(2, 6), breaks = 0.5 * 1:3) +
    scale_x_continuous("Trait 1", breaks = 0:2) +
    scale_y_continuous("Trait 2", breaks = 0:2) +
    coord_equal(xlim = c(0, 2.25), ylim = c(0, 2.25)) +
    facet_wrap(~ id, nrow = 3, labeller = label_parsed) +
    theme(strip.text = element_text(size = 8),
          panel.border = element_rect(size = 0.5, fill = NA)) +
    NULL


save_plot(trait_outcomes_p, 6.5, 6, "1-")





# =======================================================================================`
# =======================================================================================`

#               Coexistence ----

# =======================================================================================`
# =======================================================================================`


set.seed(532052696)
ds <- exp(runif(.q, -6, -2))

eta <- matrix(0, .q, .q)
eta[lower.tri(eta)] <- abs(etas)
eta <- eta + t(eta)
diag(eta) <- 1



one_d_combo <- function(sign1, sign2, sign3) {

    .d <- ds * c(sign1, sign2, sign3)

    Z <- quant_gen(q = .q, eta = eta, d = .d, max_t = 20e3L, n_reps = 24,
              save_every = 0L, n = 100, N0 = rep(1, 100),
              start_t = 0, perturb_sd = 2, n_cores = 3,
              show_progress = FALSE)

    Z$nv %>%
        mutate(trait = paste0("V", trait)) %>%
        spread(trait, value) %>%
        # filter(unq_spp_filter(V1, V2, V3)) %>%
        mutate(d1 = .d[1], d2 = .d[2], d3 = .d[3])

}


set.seed(751678518)
d_sim_df <- crossing(sign1 = c(-1,1), sign2 = c(-1,1), sign3 = c(-1,1)) %>%
    pmap_dfr(one_d_combo)

coexistence1_p <- d_sim_df %>%
    group_by(d1, d2, d3, rep) %>%
    summarize(N = n()) %>%
    # group_by(d1, d2, d3) %>%
    # summarize(N = list(unique(N))) %>%
    ungroup() %>%
    # unnest(cols = N) %>%
    mutate_at(vars(starts_with("d")),
              ~ factor(sign(.x), levels = c(-1,1), labels = c("-", "+"))) %>%
    mutate(id = interaction(d1, d2, d3, sep = "{}"),
           id = factor(id, levels = levels(id),
                       labels = sprintf("%s{}", levels(id)))) %>%
    ggplot(aes(id, N)) +
    geom_jitter(aes(color = id), width = 0.25, alpha = 0.5) +
    theme(strip.text = element_text(size = 8)) +
    scale_y_continuous("Number of surviving species") +
    scale_x_discrete(expression("Sign of" ~ d[1] * "," ~ d[2] * "," ~ d[3]),
                     labels = rlang::parse_exprs) +
    scale_color_viridis_d(guide = FALSE, end = 0.8, option = "B") +
    theme(panel.grid.major.x = element_line(color = "gray80", size = 0.5),
          axis.text.x = element_text(size = 11)) +
    NULL

save_plot(coexistence1_p, 6.5, 3, "2-")




# Keep 2 positive, make the last one go slowly negative
one_d_combo2 <- function(.d3, .max_t) {

    .d <- c(0.1, 0.1, .d3)

    Z <- quant_gen(q = .q, eta = eta, d = .d, max_t = .max_t, n_reps = 24,
                   save_every = 0L, n = 100, N0 = rep(1, 100),
                   start_t = 0, perturb_sd = 2, n_cores = 3,
                   show_progress = FALSE)

    Z$nv %>%
        spread(trait, value) %>%
        group_by(rep) %>%
        summarize(n_spp = n()) %>%
        ungroup() %>%
        mutate(d3 = .d3) %>%
        identity()

}

# Takes ~5 min
set.seed(807582316)
d_sim_df2 <- tibble(.d3 = c(-10^(c(-2, -3, -4)), 10^(c(-4, -3, -2))),
                    .max_t = c(rep(200e3L, 2), 2e6L, rep(200e3L, 3))) %>%
    pmap_dfr(one_d_combo2)


coexistence2_p <- d_sim_df2 %>%
    mutate(d3 = factor(d3, levels = d3 %>% sort() %>% unique(),
                       labels = d3 %>% sort() %>% unique() %>%
                           {sprintf("%s10^{%i}", sign(.) %>% paste() %>% str_remove("1"),
                                    abs(.) %>% log10())})) %>%
    ggplot(aes(d3, n_spp)) +
    geom_jitter(aes(color = rep), height = 0, width = 0.2, alpha = 0.75) +
    scale_color_viridis_d(guide = FALSE, end = 0.95, option = "D") +
    scale_x_discrete(expression("Value of" ~ d[3]), labels = rlang::parse_exprs) +
    scale_y_continuous("Number of surviving species") +
    NULL


save_plot(coexistence2_p, 6.5, 3, "3-")






# =======================================================================================`
# =======================================================================================`

#               Invasibility ----

# =======================================================================================`
# =======================================================================================`



invade_test <- function(r, .V1, .V2, .dsigns, .V3 = 0) {
    Z <- quant_gen(q = .q, eta = eta, d = ds*.dsigns, max_t = 20e3L, n_reps = 1,
                   save_every = 100L, n = 2, N0 = c(10.9196300066, 10.9196300066 * 0.1),
                   V0 = list(cbind(2, 0, 0), cbind(.V1, .V2, .V3)),
                   start_t = 0, perturb_sd = 0,
                   show_progress = FALSE) %>%
        .[["nv"]] %>%
        mutate(trait = paste0("V", trait)) %>%
        spread(trait, value) %>%
        mutate(spp = factor(spp, levels = 1:2, labels = c("native", "invader")),
               rep = paste(.V1, .V2, .V3, sep = "__"),
               dsigns = paste(.dsigns, collapse = "__"))
    return(Z)
}


# Takes ~2 min
invade_df <- crossing(.V1 = seq(0, 4, 0.1),
                      .V2 = seq(0, 4, 0.1),
                      .dsigns = list(c(-1,1,1), c(1,1,1))) %>%
    mutate(r = 1:n()) %>%
    pmap_dfr(invade_test)

invasion_p_df <- invade_df %>%
    mutate(dsigns = factor(dsigns, levels = c("-1__1__1", "1__1__1"),
                           labels = c("-{}+{}+{}", "+{}+{}+{}"))) %>%
    group_by(dsigns, rep) %>%
    summarize(invaded = factor(sum(time == max(time) & spp == "invader"),
                               levels = 0:1, labels = c("no", "yes")),
              V1 = str_split(rep, "__")[[1]][[1]] %>% as.numeric(),
              V2 = str_split(rep, "__")[[1]][[2]] %>% as.numeric()) %>%
    ungroup() %>%
    select(-rep)


invasion_p <- invasion_p_df %>%
    ggplot(aes(V1, V2)) +
    geom_tile(aes(fill = invaded)) +
    geom_point(data = tibble(V1 = c(0, 2, 2), V2 = c(2, 0, 0),
                             dsigns = c("+{}+{}+{}", "+{}+{}+{}", "-{}+{}+{}")),
               size = 4, color = "black", shape = 4) +
    facet_wrap(~ dsigns, labeller = label_parsed, nrow = 1) +
    xlab("Invader trait 1") +
    ylab("Invader trait 2") +
    scale_fill_viridis_d("successful", begin = 0.3, end = 0.9, option = "C") +
    coord_equal()


save_plot(invasion_p, 6.5, 3, "4-")




# =======================================================================================`
# =======================================================================================`

#               Conditional coexistence ----

# =======================================================================================`
# =======================================================================================`



# d_sim_df %>% distinct(d1, d2, d3)

d_sim_df %>%
    filter(d1 < 0, d2 > 0, d3 > 0, rep == 1) %>%
    select(-starts_with("d"), -rep) %>%
    .[["N"]]


set.seed(1245262468)
V0_coexist <- lapply(1:100,
                     function(i) {
                         V0 <- rbind(numeric(3))
                         V0[,2] <- abs(rnorm(1, 0, 2))
                         V0[,1] <- runif(1) * V0[,2]
                         return(V0)
                     })
V0_exclude <- lapply(1:100,
                     function(i) {
                         V0 <- rbind(abs(rnorm(3, 0, 2)))
                         V0[,3] <- 0
                         return(V0)
                     })


cond_coexist_test <- function(.V0, .lab) {
    Z <- quant_gen(q = .q, eta = eta, d = ds*c(-1,1,1), max_t = 20e3L, n_reps = 1,
                   save_every = 10L, n = 100, N0 = rep(1, 100),
                   V0 = .V0,
                   start_t = 0, perturb_sd = 0,
                   show_progress = FALSE) %>%
        .[["nv"]] %>%
        mutate(trait = paste0("V", trait)) %>%
        spread(trait, value) %>%
        select(-rep) %>%
        mutate(V0 = .lab)
    return(Z)
}




cond_coexist_df <- tibble(.V0 = list(V0_coexist, V0_exclude),
                          .lab = c("coexistence", "exclusion")) %>%
    pmap_dfr(cond_coexist_test)





cond_coexist_p1 <- cond_coexist_df %>%
    ggplot(aes(time, N, color = spp)) +
    geom_line(alpha = 0.5) +
    facet_wrap(~ V0, scales = "free") +
    scale_color_viridis_d(begin = 0.1, end = 0.9, option = "C", guide = FALSE) +
    ylab("Abundance") +
    xlab("Time")


cond_coexist_p2 <- cond_coexist_df %>%
    arrange(time, spp) %>%
    ggplot(aes(V1, V2, color = spp)) +
    # geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray70") +
    geom_polygon(data = tibble(x = c(-1, 5, -1), y = c(-1, 5, 5)), aes(x,y),
                 color = NA, fill = "gray80") +
    geom_path(alpha = 0.5) +
    geom_point(data = cond_coexist_df %>% filter(time == min(time)),
               shape = 1) +
    geom_point(data = cond_coexist_df %>% filter(time == max(time)) %>%
                   group_by(V0) %>%
                   filter(unq_spp_filter(V1, V2)) %>%
                   ungroup(),
               shape = 4, size = 4, color = "black") +
    scale_color_viridis_d(begin = 0.1, end = 0.9, option = "C", guide = FALSE) +
    facet_wrap(~ V0) +
    coord_equal(xlim = c(0, 4.15), ylim = c(0, 4.15)) +
    ylab("Trait 2") +
    xlab("Trait 1")



cond_coexist_p <- plot_grid(cond_coexist_p1, cond_coexist_p2, ncol = 1)


save_plot(cond_coexist_p, 6.5, 6, "5-")
