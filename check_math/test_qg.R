
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
get_sim <- function(eta, d) {
    return(readRDS(sprintf("results/%s-eta_%s-d.rds", eta, d)))
}



#' | eta  |  d   | results                                                              |
#' |:----:|:----:|:---------------------------------------------------------------------|
#' | neg  | neg  | two alternative states, no coexistence                               |
#' | neg  | zero | two alternative states, with coexistence                             |
#' | neg  | pos  | two alternative states, with coexistence                             |
#' | zero | neg  | neutrally stable "shell", no coexistence                             |
#' | zero | zero | neutrally stable "shell", with coexistence                           |
#' | zero | pos  | neutrally stable "shell", with coexistence                           |
#' | pos  | neg  | neutrally stable ring, no coexistence †                              |
#' | pos  | zero | neutrally stable ring, with coexistence                              |
#' | pos  | pos  | neutrally stable ring, with coexistence                              |
#'
#' † Default parameterization
#'



qg <- get_sim("n", "n")
qg %>%
    .[["nv"]] %>%
    filter(time == max(time)) %>% .[["value"]] %>% round(digits = 8) %>% unique()


V <- qg %>%
    .[["nv"]] %>%
    filter(rep == 1, time == 1) %>% {
        N <<- {distinct(., spp, N)}$N
        .
    } %>%
    mutate(trait = factor(paste0("T", paste(trait)))) %>%
    spread("trait", "value") %>%
    select(starts_with("T", FALSE)) %>%
    as.matrix()


qg %>%
    .[["nv"]] %>%
    filter(time == max(time)) %>%
    mutate(trait = factor(paste0("T", paste(trait)))) %>%
    spread("trait", "value") %>%
    ggplot(aes(T1, T2, size = T3)) +
    geom_point(shape = 16, alpha = 0.5) +
    scale_color_manual(values = c("dodgerblue", "firebrick")) +
    NULL


qg %>%
    .[["nv"]] %>%
    filter(time == max(time)) %>%
    mutate(trait = factor(paste0("T", paste(trait)))) %>%
    spread("trait", "value") %>%
    ggplot(aes(T1, T2, size = T3)) +
    geom_point(shape = 16, alpha = 0.5) +
    facet_wrap(~ rep, nrow = 10) +
    scale_color_manual(values = c("dodgerblue", "firebrick")) +
    scale_size_continuous(range = c(0.5, 2), breaks = -1:1) +
    theme(strip.background = element_blank(), strip.text = element_blank()) +
    NULL


qg %>%
    .[["nv"]] %>%
    filter(time == max(time)) %>%
    group_by(rep) %>%
    summarize(N = length(unique(spp))) %>%
    .[["N"]] %>%
    range()


qg %>%
    .[["nv"]] %>%
    filter(rep == 1, time < 1e4) %>%
    distinct(time, spp, N) %>%
    ggplot(aes(time, N)) +
    geom_line(aes(group = spp), alpha = 0.5) +
    # geom_line(aes(y = value, group = spp), alpha = 0.5) +
    # facet_wrap(~ ) +
    NULL


#'
#' #' Create list for running a set of simulations:
#' #'
#' quant_gen_pars <- function(eta_sign, d_sign, n = 100, q = 3) {
#'     # Other parameter values that I'll keep constant:
#'     args <- list(
#'         n_reps = 100,
#'         V0 = rep(list(matrix(0, 1, q)), n),
#'         N0 = rep(1, n),
#'         f = 0.1,  # cost of the trait on the growth rate
#'         g = 0.5,  # benefit of the trait on density dependence
#'         r0 = 0.5,
#'         add_var = rep(0.5, n),
#'         mut_sd = 1,
#'         keep_pos = FALSE,
#'         start_t = 0,
#'         max_t = 1e6,
#'         min_N = 1e-4,
#'         save_every = 1000,
#'         show_progress = TRUE,
#'         n_cores = 4)
#'     # the non-additive effects of traits on `r`:
#'     if (grepl("^p", eta_sign, TRUE)) {
#'         args$eta <- 0.01
#'     } else if (grepl("^z", eta_sign, TRUE)) {
#'         args$eta <- 0
#'     } else {
#'         args$eta <- -0.01
#'     }
#'     # changes how the focal line is affected by other lines' trait values:
#'     if (grepl("^p", d_sign, TRUE)) {
#'         args$d <- 1e-4
#'     } else if (grepl("^z", d_sign, TRUE)) {
#'         args$d <- 0
#'     } else {
#'         args$d <- -0.1
#'     }
#'
#'     return(args)
#' }
#'
#'
#' set.seed(1713696743)
#' for (e in c("n", "z", "p")) {
#'     for (d in c("n", "z", "p")) {
#'         cat(sprintf("\n\neta: %s || d: %s \n", e, d))
#'         args <- quant_gen_pars(e, d)
#'         qg <- do.call(quant_gen, args)
#'         saveRDS(qg, sprintf("results/%s-eta_%s-d.rds", e, d))
#'     }
#' }



jacobian <- function(N, V) {
    args_ <- c(list(N = N, V = V),
               args[c("f", "g", "d", "add_var")])
    args_$add_var <- args_$add_var[1:length(N)]  # it's assumed they're all the same
    args_$C <- matrix(args$eta, length(V[[1]]), length(V[[1]]))
    diag(args_$C) <- 1
    jacobian_mat <- do.call(sauron:::jacobian_cpp, args_)
    return(jacobian_mat)
}


jacobs <- qg %>%
    .[["nv"]] %>%
    filter(time == max(time)) %>%
    group_by(rep, spp) %>%
    summarize(N = N[[1]],
              V = list(value)) %>%
    group_by(rep) %>%
    summarize(N = list(N),
              V = list(V)) %>%
    ungroup() %>%
    mutate(jacobs = map2(N, V, ~ jacobian(.x, .y))) %>%
    .[["jacobs"]] %>%
    identity()

eigens <- map_dbl(jacobs, ~ eigen(.x)$values[1])

set.seed(999)
qg2 <- qg %>%
    .[["nv"]] %>%
    filter(time == max(time)) %>%
    mutate(value = runif(n(), -2, 2))
jacobs2 <- qg2 %>%
    group_by(rep) %>%
    group_by(rep, spp) %>%
    summarize(N = N[[1]],
              V = list(value)) %>%
    group_by(rep) %>%
    summarize(N = list(N),
              V = list(V)) %>%
    ungroup() %>%
    mutate(jacobs = map2(N, V, ~ jacobian(.x, .y))) %>%
    .[["jacobs"]] %>%
    identity()
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


