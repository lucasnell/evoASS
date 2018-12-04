
suppressPackageStartupMessages({
    library(dplyr)
    library(magrittr)
    library(tidyr)
    library(ggplot2)
    library(purrr)
    library(evoASS)
})
# Used in my system for better plotting (if statement checks for macOS):
if ((!is.null(Sys.info()[["sysname"]]) && Sys.info()[["sysname"]] == "Darwin") ||
    (!is.null() && grepl("^darwin", R.version$os))) {
    options("device" = "quartz")
    grDevices::graphics.off()
}


n <- 100
q <- 4
args <- list(
    n_reps = 100,
    V0 = lapply(1:n, function(i) matrix(1, 1, q)),
    N0 = rep(1, n),
    f = 0.1,  # cost of the trait on the growth rate
    g = 0.5,  # benefit of the trait on density dependence
    eta = 0.9,  # the non-additive effects of traits on `r`
    r0 = 0.5,
    d = -0.1,  # changes how the focal line is affected by other lines' trait values
    add_var = rep(0.1, n),
    mut_sd = 0.1,
    keep_pos = FALSE,
    start_t = 100,
    max_t = 1e5,
    min_N = 1e-4,
    save_every = 100,
    show_progress = TRUE,
    n_cores = 4)

qg <- do.call(quant_gen, args)
qg

# evoASS:::F_t_cpp(V0, N0, f, g, {C <- matrix(0.9, q, q); diag(C) <- 1; C}, r0, d)[[1]]


# # With no perturbation:
# qg0 <- quant_gen(n_reps, V0, N0,
#                  f, g, eta, r0, d, add_var, delta * 0, start_t,
#                  max_t, min_N, save_every, show_progress, n_cores)
# # With no perturbation and starting traits at -1:
# qg01 <- quant_gen(n_reps, V0 = lapply(V0, `*`, e2 = -1), N0,
#                   f, g, eta, r0, d, add_var, delta * 0, start_t,
#                   max_t, min_N, save_every, show_progress, n_cores)

annot_ <- annotate(geom = "text",
         x = filter(qg$nv, time == args$max_t, trait == "1")$value %>% range() %>% mean(),
         y = filter(qg$nv, time == args$max_t, trait == "2")$value %>% range() %>% mean(),
         hjust = 0.5, vjust = 0.5,
          label = with(args, sprintf(paste(c("f =", "g =", "eta =", "r0 =",
                                  "d =", "add_var ="),
                                "%.2g", sep = " ", collapse = "\n"),
                          f, g, eta, r0, d, add_var[1])))



# qgp <-
qg %>%
    .[["nv"]] %>%
    filter(time == args$max_t) %>%
    mutate(trait = factor(paste0("T", paste(trait)))) %>%
    spread("trait", "value") %>%
    ggplot(aes(T1, T2, size = T3, color = T4)) +
    geom_point(shape = 16, alpha = 0.75) +
    scale_color_continuous(low = "gray80", high = "black") +
    NULL


qgp <- qg %>%
    .[["nv"]] %>%
    filter(time == args$max_t) %>%
    mutate(trait = factor(paste0("T", paste(trait)))) %>%
    spread("trait", "value") %>%
    distinct(T1, T2, T3, .keep_all = TRUE) %>%
    ggplot(aes(T1, T2)) +
    geom_point(shape = 1) +
    scale_color_brewer(palette = "Dark2") +
    scale_linetype_manual(values = rep(1, n), guide = FALSE) +
    # scale_x_continuous(limits = c(0, NA)) +
    # scale_y_continuous(limits = c(0, NA)) +
    NULL


qgp +
    geom_point(aes(T1, T3), color = "dodgerblue", shape = 1) +
    # geom_line(aes(T1, T2, group = rep)) +
    # geom_point(aes(T1, T3, group = rep)) +
    # geom_point(aes(T2, T3, group = rep)) +
    # annot_ +
    NULL

qg %>%
    .[["nv"]] %>%
    filter(time == args$max_t) %>%
    mutate(trait = as.integer(trait)) %>%
    distinct(rep, time, trait, value, .keep_all = TRUE) %>%
    ggplot(aes(trait, value)) +
    geom_line(aes(group = rep), alpha = 0.5) +
    geom_point(shape = 1) +
    scale_color_brewer(palette = "Dark2") +
    scale_linetype_manual(values = rep(1, n), guide = FALSE)

# qgp +
#     geom_point(data = qg0$nv %>%
#                    filter(time == max_t) %>%
#                    distinct(trait, value, .keep_all = TRUE) %>%
#                    mutate(trait = factor(paste0("T", paste(trait)))) %>%
#                    spread("trait", "value"),
#                shape = 1) +
#     geom_point(data = qg01$nv %>%
#                    filter(time == max_t) %>%
#                    distinct(trait, value, .keep_all = TRUE) %>%
#                    mutate(trait = factor(paste0("T", paste(trait)))) %>%
#                    spread("trait", "value"),
#                shape = 1)
#
# qg0$nv %>%
#     filter(time == max_t) %>%
#     mutate(trait = factor(paste0("T", paste(trait)))) %>%
#     spread("trait", "value") %>%
#     ggplot(aes(T1, T2)) +
#     geom_point() +
#     scale_color_brewer(palette = "Dark2") +
#     scale_linetype_manual(values = rep(1, n), guide = FALSE) +
#     # scale_x_discrete(breaks = seq(10, n_reps, 10)) +
#     NULL
#
#
#
# scale_ <- scale_x_continuous(bquote(.("Time (") %*% 1000 * .(")")),
#                              breaks = seq(0, max_t, max_t / 2),
#                              labels = seq(0, max_t / 1000, max_t / 2000))
#
# r = 1
# qg$nv %>%
#     # filter(rep == r, trait == 1) %>%
#     ggplot(aes(time, value)) +
#     geom_line(aes(group = interaction(spp, rep)), alpha = 0.5) +
#     facet_wrap(~ trait) +
#     scale_color_brewer(palette = "Dark2") +
#     scale_linetype_manual(values = rep(1, n), guide = FALSE) +
#     scale_
#
#
# qg$nv %>%
#     # filter(rep == r) %>%
#     distinct(rep, time, spp, .keep_all = TRUE) %>%
#     ggplot(aes(time, N)) +
#     geom_line(aes(group = rep), alpha = 0.25) +
#     facet_wrap(~ spp) +
#     scale_
#
#
qg$nv %>%
    filter(time == max_t) %>%
    ggplot(aes(value)) +
    # geom_point(position = position_jitter(height = 0)) +
    geom_freqpoly(color = "dodgerblue", size = 0.75, bins = 20) +
    geom_rug() +
    facet_wrap(~ trait, ncol = 1) +
    scale_color_brewer(palette = "Dark2") +
    scale_linetype_manual(values = rep(1, n), guide = FALSE) +
    # scale_x_discrete(breaks = seq(10, n_reps, 10)) +
    NULL


