

#' Seeds for quantitative genetics simulations.
#'
#' By storing seeds for each simulation, I can retrieve individual runs to look
#' at them in greater detail.
#'
#' @format A tibble with 18 rows and 4 variables:
#' \describe{
#'   \item{n_traits}{Number of traits: 2 or 3.}
#'   \item{eta_par}{Sign of `eta` parameter: -1, 0, or 1.}
#'   \item{d_par}{Sign of `d` parameter: -1, 0, or 1.}
#'   \item{seed}{Seed to use for this particular set of simulations.}
#' }
#'
#'
#' @section How data were generated:
#'
#' ```r
#' library(dplyr)
#' qg_seeds <- tibble(n_traits = rep(2:3, each = 9),
#'                    eta_par = rep(rep(-1:1, each = 3), 2),
#'                    d_par = rep(rep(-1:1, 3), 2),
#'                    seed = numeric(18))
#' set.seed(2028111205)
#' qg_seeds$seed[qg_seeds$n_traits==2] <- sample.int(.Machine$integer.max, 9)
#' set.seed(1713696743)
#' qg_seeds$seed[qg_seeds$n_traits==3] <- sample.int(.Machine$integer.max, 9)
#'
#' usethis::use_data(qg_seeds)
#' ```
#'
"qg_seeds"


