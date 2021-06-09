#' Compare multiple multidimensional reliability estimates
#'
#' This function calculates multidimensional reliability coefficients such as
#' maximal reliability, correlated factors reliability, stratified alpha,
#' multidimensional parallel reliability, second-order factor reliability,
#' bifactor reliability, and related omega hierarchical and sub-dimensional
#' reliability coefficients.
#'
#' @param x observed item scores or their covariances
#' @param until The number of items up to the first sub-construct
#' @return several multidimensional reliability, omega hiearchical, and subdimensional reliability
#'
multirel <- function(x, until) {
  out <- list(maximal_reliability = maximal_reliability(x, until),
              correlated_factors = correlated_factors(x, until),
              stratified_alpha = stratified_alpha(x, until),
              multi_parallel = multi_parallel(x, until),
              second_order = second_order(x, until),
              bifactor = bifactor(x, until))
  invisible(out)
}
