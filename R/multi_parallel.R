#' Obtain multidimensional parallel reliability
#'
#' Multidimensional parallel reliability is derived from the multidimensional parallel model (Cho, 2016).
#' This is equivalent to entering the correlation instead of the covariance into the stratified alpha formula.
#'
#' @param x observed item scores or their covariances
#' @param until The number of items up to the first sub-construct
#' @param print If TRUE, the result is printed to the screen.
#' @return a multidimensional parallel reliability estimate
#' @references Cho, E. (2016). Making reliability reliable: A systematic
#' approach to reliability coefficients. Organizational Research Methods, 19(4),
#'  651–682.
multi_parallel <- function(x, until, print = FALSE) {
  r <- stats::cov2cor(get_cov(x))
  return(stratified_alpha(r, until, mp = TRUE, print = print))
}
