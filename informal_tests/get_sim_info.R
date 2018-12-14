get_sim_info <- function(sim_i) {

    sims <- readr::read_csv("informal_tests/simulated_data.csv",
                            col_types = readr::cols(.default = readr::col_double()))

    info <- list2env(as.list(sims[sim_i,c("f", "g", "r0", "d", "eta")]))

    with(info, {
        N = as.numeric(sims[sim_i,colnames(sims)[grepl("^N", colnames(sims))]])
        V = matrix(as.numeric(sims[sim_i,colnames(sims)[grepl("^V", colnames(sims))]]),
                   length(N))
        C = matrix(eta, ncol(V), ncol(V))
        diag(C) = 1
        CCC = C + t(C)
        sigma2 = 0.01
    })

    return(as.list(info))
}
