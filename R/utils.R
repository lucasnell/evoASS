

#' Black theme to use for presentations
#'
#' Implemented only for `adapt_dyn` objects thus far.
#'
#' @param base_size base font size
#' @param base_family base font family
#' @export
#'
theme_black = function(base_size = 12, base_family = "") {

    theme_classic(base_size = base_size, base_family = base_family) %+replace%

        theme(
            # Specify axis options
            axis.line = element_blank(),
            axis.text.x = element_text(size = base_size*0.8, color = "white",
                                       lineheight = 0.9),
            axis.text.y = element_text(size = base_size*0.8, color = "white",
                                       lineheight = 0.9),
            axis.ticks = element_line(color = "white", size  =  0.2),
            axis.title.x = element_text(size = base_size, color = "white",
                                        margin = margin(10, 0, 0, 0)),
            axis.title.y = element_text(size = base_size, color = "white", angle = 90,
                                        margin = margin(0, 10, 0, 0)),
            axis.ticks.length = unit(0.3, "lines"),
            # Specify legend options
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.margin = margin(0,0,0,0),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = NULL,
            legend.key.width = NULL,
            legend.text = element_text(size = base_size*0.8, color = "white"),
            legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0,
                                        color = "white"),
            legend.position = "right",
            legend.text.align = NULL,
            legend.title.align = NULL,
            legend.direction = "vertical",
            legend.box = NULL,
            # Specify panel options
            panel.background = element_rect(fill = "black", color  =  NA),
            panel.border = element_rect(fill = NA, color = "white"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(0.5, "lines"),
            # Specify facetting options
            strip.background = element_blank(),
            strip.text.x = element_text(size = base_size*0.8, color = "white"),
            strip.text.y = element_text(size = base_size*0.8, color = "white",
                                        angle = -90),
            # Specify plot options
            plot.background = element_rect(color = "black", fill = "black"),
            plot.title = element_text(size = base_size*1.2, color = "white"),
            plot.margin = margin(0,0,0,0)

        )

}



#' Random deviates from truncated normal distribution
#'
#' Note: truncated at zero.
#'
#' @param n Number of deviates.
#' @param mu Mean of distribution.
#' @param sigma Standard deviation of distribution.
#'
#' @return
#' @export
#'
trnorm <- function(n, mu, sigma) {

    stopifnot(inherits(mu, "numeric") && inherits(sigma, "numeric"))
    stopifnot(length(mu) %in% c(1, n))
    stopifnot(length(sigma) %in% c(1, n))

    ml <- length(mu)
    sl <- length(sigma)

    if (ml == 1 && sl == 1) {
        out <- trunc_rnorm_cpp(n, mu, sigma)
    } else if (ml > 1 && sl == 1) {
        out <- trunc_rnorm_mu_cpp(mu, sigma)
    } else if (ml == 1 && sl > 1) {
        out <- trunc_rnorm_sigma_cpp(mu, sigma)
    } else {
        out <- trunc_rnorm_mu_sigma_cpp(mu, sigma)
    }

    return(out)

}
