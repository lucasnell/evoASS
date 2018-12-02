

n <- 40
q <- 3

qg <- evoASS:::quantgen_cpp(n_reps = 2,
                            V0 = lapply(1:n, function(i) matrix(1, 1, q)),
                            N0 = rep(10, n),
                            f = 0.1,
                            g = 0.5,
                            nu = 0.05,
                            r0 = 0.25,
                            d = -0.001,
                            add_var = rep(1, n),
                            delta = 1,
                            start_t = 1e2,
                            max_t = 1e4,
                            min_N= 1e-4,
                            save_every = 100,
                            show_progress = FALSE,
                            n_cores = 4)

str(qg)


dplyr::as_data_frame(qg$N_V)
