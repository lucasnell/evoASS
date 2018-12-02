
#' Generic method to perturb an object.
#'
#' Implemented only for `adapt_dyn` objects thus far.
#'
#' @param obj The object to perturb.
#' @param ... Additional arguments.
#' @export
#'
perturb <- function(obj, ...) UseMethod("perturb")
