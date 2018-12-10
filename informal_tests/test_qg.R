
suppressPackageStartupMessages({
    library(dplyr)
    library(magrittr)
    library(tidyr)
    library(ggplot2)
    library(purrr)
    library(sauron)
})
# Used in my system for better plotting (if statement checks for macOS):
if ((!is.null(Sys.info()[["sysname"]]) && Sys.info()[["sysname"]] == "Darwin") ||
    (!is.null() && grepl("^darwin", R.version$os))) {
    options("device" = "quartz")
    grDevices::graphics.off()
}


n <- 100
q <- 3
args <- list(
    n_reps = 100,
    V0 = lapply(1:n, function(i) matrix(1, 1, q)),
    N0 = rep(1, n),
    f = 0.1,  # cost of the trait on the growth rate
    g = 0.5,  # benefit of the trait on density dependence
    eta = 0.01,  # the non-additive effects of traits on `r`
    r0 = 0.5,
    d = -0.1,  # changes how the focal line is affected by other lines' trait values
    add_var = rep(0.1, n),
    mut_sd = 0.1,
    keep_pos = FALSE,
    start_t = 100,
    max_t = 1e6,
    min_N = 1e-4,
    save_every = 1000,
    show_progress = TRUE,
    n_cores = 4)

qg <- do.call(quant_gen, args)
qg


qg %>%
    .[["nv"]] %>%
    filter(time == max(time))

qg %>%
    .[["nv"]] %>%
    filter(time == max(time)) %>%
    mutate(trait = factor(paste0("T", paste(trait)))) %>%
    spread("trait", "value") %>%
    ggplot(aes(T1, T2, size = T3)) +
    geom_point(shape = 16, alpha = 0.5) +
    scale_color_manual(values = c("dodgerblue", "firebrick")) +
    NULL

# # Start with one, see if it goes anywhere:
# re_quantgen <- function(N_, V_, args_, ...) {
#     other_args <- list(...)
#     args_[names(other_args)] <- other_args
#     args_[["N0"]] <- N_
#     if (!inherits(V_, "list")) V_ <- list(V_)
#     args_[["V0"]] <- V_
#     args_[["mut_sd"]] <- 0
#     args_[["start_t"]] <- 0
#     args_[["save_every"]] <- 0
#     do.call(quant_gen, args_)$nv
# }
# # Start with all, see where they go:
# combined_quantgen <- function(qg_obj, args_, ...) {
#
#     other_args <- list(...)
#     args_[names(other_args)] <- other_args
#
#     NV_df <- qg_obj %>%
#         .[["nv"]] %>%
#         filter(time == max(time)) %>%
#         mutate(trait = paste0("V", trait)) %>%
#         spread("trait", "value") %>%
#         dplyr::select(N, starts_with("V", ignore.case = FALSE))
#
#     N_ <- pmax(NV_df$N / nrow(NV_df), 0.1)
#
#     V_ <- NV_df %>%
#         dplyr::select(-N) %>%
#         as.matrix() %>%
#         split(row(.))
#
#     args_[["N0"]] <- N_
#     args_[["V0"]] <- V_
#
#     return(do.call(quant_gen, args_))
# }
#
# qg_c <- combined_quantgen(qg, args, max_t = 1e5, save_every = 0, mut_sd = 0.001)
#
# qg_c$nv %>%
#     filter(rep == 1, trait == 1) %>%
#     ggplot(aes(time, value)) +
#     geom_line(aes(group = interaction(spp, rep)), alpha = 0.5) +
#     facet_wrap(~ trait) +
#     scale_color_brewer(palette = "Dark2") +
#     scale_linetype_manual(values = rep(1, n), guide = FALSE)
#
# qg_c %>%
#     .[["nv"]] %>%
#     mutate(trait = factor(paste0("T", paste(trait)))) %>%
#     spread("trait", "value") %>%
#     ggplot(aes(T1, T2, size = T3)) + # , color = T4)) +
#     geom_point(shape = 16, alpha = 0.5) +
#     scale_color_continuous(low = "gray80", high = "black") +
#     NULL

jacobian <- function(N, V, eps = -1e-6) {
    args_ <- c(list(N = N[[1]], V = V, eps = eps),
               args[c("f", "g", "r0", "d", "add_var")])
    args_$add_var <- args_$add_var[1:length(N)]  # it's assumed they're all the same
    args_$C <- matrix(args$eta, length(V[[1]]), length(V[[1]]))
    diag(args_$C) <- 1
    jacobian_mat <- do.call(sauron:::jacobian_cpp, args_)
    return(jacobian_mat)
}


# # Start with one, see if it goes anywhere:
# re_quantgen <- function(N_, V_, args_, ...) {
#     other_args <- list(...)
#     args_[names(other_args)] <- other_args
#     args_[["N0"]] <- N_
#     if (!inherits(V_, "list")) V_ <- list(V_)
#     args_[["V0"]] <- V_
#     args_[["mut_sd"]] <- 0
#     args_[["start_t"]] <- 0
#     args_[["save_every"]] <- 0
#     args_[["add_var"]] <- args_[["add_var"]][1:length(N_)]
#     z <- do.call(quant_gen, args_)$nv
#     return(range(z$value - rep(do.call(c, V_), args_$n_reps)))
# }
#
# # All zeros:
# qg_re <- qg %>%
#     .[["nv"]] %>%
#     filter(time == max(time)) %>%
#     group_by(rep) %>%
#     summarize(N = list(unique(N)),
#               V = map(unique(spp), ~ list(value[spp == .x]))) %>%
#     mutate(diff = map2(N, V, ~ re_quantgen(.x, .y, args, show_progress = FALSE, max_t = 1e5))) %>%
#     select(rep, diff) %>%
#     unnest()
# sum(qg_re$diff != 0)


jacobs <- qg %>%
    .[["nv"]] %>%
    filter(time == max(time)) %>%
    group_by(rep) %>%
    summarize(N = list(unique(N)),
              V = map(unique(spp), ~ list(value[spp == .x]))) %>%
    mutate(jacobs = map2(N, V, ~ jacobian(.x, .y))) %>%
    .[["jacobs"]]

eigens <- map_dbl(jacobs, ~ eigen(.x)$values[1])

set.seed(999)
qg2 <- qg %>%
    .[["nv"]] %>%
    filter(time == max(time)) %>%
    mutate(value = runif(n(), -2, 2))
jacobs2 <- qg2 %>%
    group_by(rep) %>%
    summarize(N = list(unique(N)),
              V = map(unique(spp), ~ list(value[spp == .x]))) %>%
    mutate(jacobs = map2(N, V, ~ jacobian(.x, .y))) %>%
    .[["jacobs"]]
eigens2 <- map_dbl(jacobs2, ~ eigen(.x)$values[1])
qg2 <- qg2 %>%
    mutate(trait = factor(paste0("T", paste(trait)))) %>%
    spread("trait", "value") %>%
    mutate(stable = map_dbl(as.integer(rep), ~ eigens2[.x]))

# # map_dbl(jacobs, ~ eigen(.x)$values[2])
# # hist(map_dbl(jacobs, ~ eigen(.x)$values[3]))
# eigen_df <- map_dfr(jacobs, ~ as_data_frame(rbind(eigen(.x)$values)))
# sum(eigen_df[,1] > 0)
# sum(eigen_df[,2] > 0)
# sum(eigen_df[,3] > 0)
#
# jacobs[[1]]
# jacobs[1:2]
# map(jacobs[1:2], ~ eigen(.x)[c("values", "vectors")])

qg %>%
    .[["nv"]] %>%
    filter(time == max(time)) %>%
    mutate(trait = factor(paste0("T", paste(trait)))) %>%
    spread("trait", "value") %>%
    # mutate(stable = map_lgl(as.integer(rep), ~ eigens[.x] < 1)) %>%
    mutate(stable = map_dbl(as.integer(rep), ~ eigens[.x])) %>%
    ggplot(aes(T1, T2, size = T3)) +
    geom_point(aes(fill = stable), shape = 21, color = "gray50") +
    # geom_point(data = qg2, aes(fill = stable), shape = 21, color = "gray50") +
    # scale_color_manual(values = c("dodgerblue", "firebrick")) +
    scale_fill_gradient2(midpoint = 1.0, low = "dodgerblue", high = "firebrick") +
    NULL

qg2 %>%
    ggplot(aes(T1, T2, size = T3)) +
    geom_point(aes(fill = stable), shape = 21, color = "gray50") +
    # scale_color_manual(values = c("dodgerblue", "firebrick")) +
    scale_fill_gradient2(midpoint = 1.0, low = "dodgerblue", high = "firebrick") +
    NULL


# sauron:::F_t_cpp(V0, N0, f, g, {C <- matrix(0.9, q, q); diag(C) <- 1; C}, r0, d)[[1]]


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



qgp <- qg %>%
    .[["nv"]] %>%
    filter(time == args$max_t) %>%
    mutate(trait = factor(paste0("T", paste(trait)))) %>%
    spread("trait", "value") %>%
    ggplot(aes(T1, T2, size = T3)) + # , color = T4)) +
    geom_point(shape = 16, alpha = 0.5) +
    scale_color_continuous(low = "gray80", high = "black") +
    NULL
qgp


# qg %>%
#     .[["nv"]] %>%
#     filter(time == args$max_t) %>%
#     mutate(trait = as.integer(trait)) %>%
#     distinct(rep, time, trait, value, .keep_all = TRUE) %>%
#     ggplot(aes(trait, value)) +
#     geom_line(aes(group = rep), alpha = 0.5) +
#     geom_point(shape = 1) +
#     scale_color_brewer(palette = "Dark2") +
#     scale_linetype_manual(values = rep(1, n), guide = FALSE)
#
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


