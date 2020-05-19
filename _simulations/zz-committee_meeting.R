
suppressPackageStartupMessages({
    library(tidyverse)
    library(sauron)
    library(grid)
    library(cowplot)
})

.REDO_SIMS <- FALSE
.N_THREADS <- 2

source(".Rprofile")

# Where to save rds files
rds <- function(.x) {
    sprintf("~/GitHub/Wisconsin/sauron/_simulations/zz-comm_rds/%s.rds",
           gsub("\\.rds$", "", gsub("^\\/", "", .x)))
}



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

group_spp <- function(..., .prec = 0.001) {
    V <- list(...)
    V <- do.call(cbind, V)
    V <- split(V, row(V))
    return(sauron:::group_spp_cpp(V, precision = .prec))
}




.q <- 3



# =======================================================================================`
# =======================================================================================`

#               I. Trait outcomes ----

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
                          start_t = 0, perturb_sd = 2, n_threads = .N_THREADS,
                          show_progress = FALSE)

    trait_ts$nv %>%
        mutate(trait = paste0("V", trait)) %>%
        spread(trait, value) %>%
        filter(unq_spp_filter(V1, V2, V3)) %>%
        select(starts_with("V")) %>%
        mutate(eta1 = eta[2,1], eta2 = eta[3,1], eta3 = eta[3,2])
}


if (.REDO_SIMS) {
    # Takes ~2.2 min
    set.seed(145746935)
    eta_sim_df <- crossing(sign1 = -1:1, sign2 = -1:1, sign3 = -1:1) %>%
        pmap_dfr(one_eta_combo)
    saveRDS(eta_sim_df, rds("eta_sim"))
} else {
    eta_sim_df <- readRDS(rds("eta_sim"))
}


# eta_sim_df




p_group_levels <- c("all super", "all sub", "all neutral", "1 of each",
                    "2 super, 1 sub", "2 super, 1 neutral", "2 sub, 1 super", "2 sub, 1 neutral",
                    "2 neutral, 1 super", "2 neutral, 1 sub")

# Because the "1 of each" category has 6 levels but only 3 of which have unique
# trait end points, I'm creating this function to deal with this:
one_each_fun <- function(.x) {

    if (.x$p_groups[1] != "1 of each") return(mutate_at(.x, vars(id), paste))

    .x %>%
        mutate(color_group = group_spp(V1, V2, V3) + 1L) %>%
        group_by(p_groups, color_group) %>%
        summarize(V1 = V1[1],
                  V2 = V2[1],
                  V3 = V3[1],
                  id = paste(id, collapse = "\n")) %>%
        ungroup()

}

trait_outcomes_p_df <- eta_sim_df %>%
    mutate_at(vars(starts_with("eta")),
              ~ factor(sign(.x), levels = -1:1, labels = c("−", "0", "+"))) %>%
    mutate(id = interaction(eta1, eta2, eta3, sep = "")) %>%
    group_by(id) %>%
    filter(unq_spp_filter(V1, V2, V3, .prec = 0.1)) %>%
    ungroup() %>%
    mutate(n_p = str_count(id, "\\+"),
           n_n = str_count(id, "\\−"),
           n_z = str_count(id, "0"),
           p_groups = case_when(n_p == 3 ~ p_group_levels[1],
                                n_n == 3 ~ p_group_levels[2],
                                n_z == 3 ~ p_group_levels[3],
                                n_p == 1 & n_z == 1 & n_n == 1 ~ p_group_levels[4],
                                n_p == 2 & n_n == 1 ~ p_group_levels[5],
                                n_p == 2 & n_z == 1 ~ p_group_levels[6],
                                n_n == 2 & n_p == 1 ~ p_group_levels[7],
                                n_n == 2 & n_z == 1 ~ p_group_levels[8],
                                n_z == 2 & n_p == 1 ~ p_group_levels[9],
                                n_z == 2 & n_n == 1 ~ p_group_levels[10],
                                TRUE ~ NA_character_)) %>%
    select(-starts_with("n", ignore.case = FALSE)) %>%
    # filter(is.na(p_groups))  ## <-- check that this yields no rows
    mutate(p_groups = factor(p_groups, levels = p_group_levels)) %>%
    split(.$p_groups) %>%
    map(~ mutate(.x, color_group = id %>% droplevels() %>% as.integer())) %>%
    map_dfr(one_each_fun) %>%
    mutate_at(vars(color_group, id), factor) %>%
    arrange(desc(V3))





trait_outcomes_p <- trait_outcomes_p_df %>%
    ggplot(aes(V1, V2, size = V3, color = color_group)) +
    geom_point(aes(fill = color_group), shape = 21) +
    geom_text(data = trait_outcomes_p_df %>%
                  filter(!grepl("^all", p_groups)) %>%
                  group_by(p_groups, color_group, id) %>%
                  summarize_at(vars(V1, V2), median) %>%
                  ungroup() %>%
                  mutate(V1 = ifelse(id == "00+", V1 - 0.2, V1),
                         V2 = ifelse(id == "0+0", V2 - 0.2, V2)),
              aes(label = id), size = 9 / 2.835, hjust = 0, vjust = 0,
              fontface = "bold", lineheight = 0.75, nudge_x = 0.2, nudge_y = 0.2) +
    scale_size_continuous("Trait 3", range = c(2, 6), breaks = 0.5 * 1:3) +
    scale_x_continuous("Trait 1", breaks = 0:2) +
    scale_y_continuous("Trait 2", breaks = 0:2) +
    coord_equal(xlim = c(0, 2.5), ylim = c(0, 2.5)) +
    facet_wrap(~ p_groups, nrow = 3) +
    theme(strip.text = element_text(size = 8),
          panel.border = element_rect(size = 0.5, fill = NA)) +
    scale_color_manual(values = viridisLite::plasma(6, end = 0.85)[c(1,3,6,2,4,5)],
                       guide = FALSE) +
    scale_fill_manual(values = viridisLite::plasma(6, end = 0.85, alpha = 0.25)[c(1,3,6,2,4,5)],
                       guide = FALSE) +
    NULL


# save_plot(trait_outcomes_p, 6.5, 5.0, "1-")





# =======================================================================================`
# =======================================================================================`

#               II. Coexistence ----

# =======================================================================================`
# =======================================================================================`


set.seed(532052696)
ds <- exp(runif(.q, -6, -2))

eta <- matrix(0, .q, .q)
eta[lower.tri(eta)] <- abs(etas)
eta <- eta + t(eta)
diag(eta) <- 1



# ------------------------------*
# __ all combinations ----
# ------------------------------*

# All combinations of - or + d values

one_d_combo <- function(sign1, sign2, sign3) {

    .d <- ds * c(sign1, sign2, sign3)

    Z <- quant_gen(q = .q, eta = eta, d = .d, max_t = 20e3L, n_reps = 24,
              save_every = 0L, n = 100, N0 = rep(1, 100),
              start_t = 0, perturb_sd = 2, n_threads = .N_THREADS,
              show_progress = FALSE)

    Z$nv %>%
        mutate(trait = paste0("V", trait)) %>%
        spread(trait, value) %>%
        # filter(unq_spp_filter(V1, V2, V3)) %>%
        mutate(d1 = .d[1], d2 = .d[2], d3 = .d[3])

}

# Takes ~11 sec
set.seed(751678518)
d_sim_df <- crossing(sign1 = c(-1,1), sign2 = c(-1,1), sign3 = c(-1,1)) %>%
    pmap_dfr(one_d_combo)


lab_df <- tibble(idd = c(1, 2, 5, 8) - 0.25, idd_e = c(1, 4, 7, 8) + 0.25,
                 N = c(rep(10, 3), 50),
                 lab = sprintf("%i", 0:3))

coexist_all_d_p <- d_sim_df %>%
    group_by(d1, d2, d3, rep) %>%
    summarize(N = n()) %>%
    ungroup() %>%
    mutate_at(vars(starts_with("d")),
              ~ factor(sign(.x), levels = c(-1,1), labels = c("-", "+"))) %>%
    mutate(id = interaction(d1, d2, d3, sep = "{}"),
           id = factor(id, levels = levels(id),
                       labels = sprintf("%s{}", levels(id))),
           n_p = str_count(id, "\\+"),
           idd = interaction(n_p, id, sep = "", lex.order = TRUE, drop = TRUE),
           idd = factor(idd, levels = levels(idd),
                        labels = gsub("0|1|2|3", "", levels(idd)))) %>%
    ggplot(aes(idd, N)) +
    geom_jitter(aes(color = id), width = 0.25, alpha = 0.5, height = 0) +
    geom_segment(data = lab_df, aes(xend = idd_e, yend = N)) +
    geom_text(data = lab_df, aes(x = (idd + idd_e) / 2, y = N + 2, label = lab),
              hjust = 0.5, vjust = 0) +
    # facet_wrap(~ n_p, scales = "free_x") +
    scale_y_continuous("Number of surviving species", limits = c(0, 55)) +
    scale_x_discrete(expression("Sign of" ~ d[1] * "," ~ d[2] * "," ~ d[3]),
                     labels = rlang::parse_exprs) +
    scale_color_viridis_d(guide = FALSE, end = 0.8, option = "B") +
    theme(strip.text = element_text(size = 8),
          axis.text.x = element_text(size = 11)) +
    NULL




# ------------------------------*
# __ slowly vary 1 ----
# ------------------------------*

# Keep 2 d values positive, make the last one go slowly negative

one_d_combo_vary_one <- function(.d3, .max_t) {

    .d <- c(0.1, 0.1, .d3)

    Z <- quant_gen(q = .q, eta = eta, d = .d, max_t = .max_t, n_reps = 24,
                   save_every = 0L, n = 100, N0 = rep(1, 100),
                   start_t = 0, perturb_sd = 2, n_threads = .N_THREADS,
                   show_progress = FALSE)

    Z$nv %>%
        spread(trait, value) %>%
        group_by(rep) %>%
        summarize(n_spp = n()) %>%
        ungroup() %>%
        mutate(d3 = .d3) %>%
        identity()

}


if (.REDO_SIMS) {
    # Takes ~5 min
    set.seed(807582316)
    one_d_sim_df <- tibble(.d3 = c(-10^(c(-2, -3, -4)), 0, 10^(c(-4, -3, -2))),
                        .max_t = 200e3L) %>%
        mutate(.max_t = ifelse(.d3 == -1e-4, 2e6L, .max_t)) %>%
        pmap_dfr(one_d_combo_vary_one)
    saveRDS(one_d_sim_df, rds("one_d_sim"))
} else {
    one_d_sim_df <- readRDS(rds("one_d_sim"))
}




lab_fun <- function(.x) {
    labs <- paste(.x)
    labs[.x != 0] <- sprintf("%s10^{%i}",
                             sign(.x[.x != 0]) %>% paste() %>% str_remove("1"),
                             abs(.x[.x != 0]) %>% log10())
    return(labs)
}

coexist_one_d_p <- one_d_sim_df %>%
    mutate(d3 = factor(d3, levels = d3 %>% sort() %>% unique(),
                       labels = d3 %>% sort() %>% unique() %>% lab_fun())) %>%
    ggplot(aes(d3, n_spp)) +
    geom_jitter(aes(color = rep), height = 0, width = 0.2, alpha = 0.75) +
    scale_color_viridis_d(guide = FALSE, end = 0.95, option = "D") +
    scale_x_discrete(expression("Value of" ~ d[3]), labels = rlang::parse_exprs) +
    scale_y_continuous("Number of surviving species") +
    NULL




# ------------------------------*
# __ combine them ----
# ------------------------------*


coexist_p <- plot_grid(NULL, coexist_all_d_p + theme(axis.title.y = element_blank()),
                       NULL, coexist_one_d_p + theme(axis.title.y = element_blank()),
                       align = "vh", ncol = 2, rel_widths = c(0.07, 1),
                       labels = c("A", "", "B", "")) +
    draw_label("Number of surviving species", x = 0.05, y = 0.5,
               vjust = 0.5, angle = 90, hjust = 0.5, size = 12)

# save_plot(coexist_p, 6.5, 6, "2-")






# =======================================================================================`
# =======================================================================================`

#               III. Invasibility ----

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


if (.REDO_SIMS) {
    # Takes ~2 min
    invade_df <- crossing(.V1 = seq(0, 4, 0.1),
                          .V2 = seq(0, 4, 0.1),
                          .dsigns = list(c(-1,1,1), c(1,1,1))) %>%
        mutate(r = 1:n()) %>%
        pmap_dfr(invade_test)
    saveRDS(invade_df, rds("invade"))
} else {
    invade_df <- readRDS(rds("invade"))
}





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


# save_plot(invasion_p, 6.5, 3, "3-")




# =======================================================================================`
# =======================================================================================`

#               IV. Conditional coexistence ----

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



# Just takes a few seconds
cond_coexist_df <- tibble(.V0 = list(V0_coexist, V0_exclude),
                          .lab = c("coexistence", "exclusion")) %>%
    pmap_dfr(cond_coexist_test)




# Time series for coexistence and exclusion
cond_coexist_ts_p <- cond_coexist_df %>%
    ggplot(aes(time, N, color = spp)) +
    geom_line(alpha = 0.5) +
    facet_wrap(~ V0, scales = "free") +
    scale_color_viridis_d(begin = 0.1, end = 0.9, option = "C", guide = FALSE) +
    ylab("Abundance") +
    scale_x_continuous("Time", breaks = seq(0, 20e3, 5e3),
                       labels = format(seq(0, 20e3, 5e3), big.mark = ",")) +
    theme(strip.text = element_text(size = 12, margin = margin(0,0,0,b = 24)))


# Movement through trait space for coexistence and exclusion:
cond_coexist_sp_p <- cond_coexist_df %>%
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
    coord_equal(xlim = c(0, 4.25), ylim = c(0, 4.25)) +
    ylab("Trait 2") +
    xlab("Trait 1") +
    theme(strip.text = element_blank())





cond_coexist_p <- plot_grid(cond_coexist_ts_p, cond_coexist_sp_p, ncol = 1,
                            rel_heights = c(0.85, 1), labels = LETTERS[1:2],
                            label_y = c(0.83,1))

cond_coexist_p

# save_plot(cond_coexist_p, 5, 5, "4-")


