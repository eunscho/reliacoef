#' Obtain stratified alpha reliability estimates of two-dimensional data
#'
#' Obtain stratified alpha reliability estimates of two-dimensional data. It
#' is a multidimensional version of coefficient alpha. Items should be grouped by each
#' sub-dimension. For example, the first dimension consists of the first item
#' through the until item.
#'
#' @param x observed item scores or their covariances
#' @param until The number of items up to the first sub-construct
#' @param mp if TRUE, print multidimensional parallel reliability estimates
#' @return multidimensional_reliability a stratified alpha reliability estimate
#' @return subdimensional_reliability reliability estimates of subdimensions
#' @return omega_hierarchical general factor saturation
#' @examples stratified_alpha(Osburn_moderate, 4)
#' @examples stratified_alpha(Sijtsma2a, c(2, 4))
#' @references Rajaratnam, N., Cronbach, L. J., & Gleser, G. C. (1965).
#' Generalizability of stratified-parallel tests. Psychometrika, 30(1), 39–56.
#' https://doi.org/10.1007/BF02289746
#' @references Cho, E. (2016). Making reliability reliable: A systematic approach
#' to reliability coefficients. Organizational Research Methods, 19(4), 651-682.
#' https://doi.org/10.1177/1094428116656239
#' @author Eunseong Cho, \email{bene@kw.ac.kr}
#'
stratified_alpha <- function(x, until, mp = FALSE) {
  m <- get_cov(x)
  n <- nrow(x)
  grp_start <- c(1, until + 1)
  grp_end <- c(until, n)

  sub_alpha <- sub_var <- vector("double", length(grp_start))
  for (i in seq_along(grp_start)) {
     sub_alpha[i] <- reliacoef::alpha(m[grp_start[i]:grp_end[i],
                                       grp_start[i]:grp_end[i]],
                                      print = FALSE)
     sub_var[i] <- sum(m[grp_start[i]:grp_end[i],
                           grp_start[i]:grp_end[i]])
  }


  sum <- vector("double", length = 1)
  for (i in seq_along(grp_start)) {
     sum <- sum + sub_var[i] * (1 - sub_alpha[i])
  }
  multi_rel <- 1 - sum / sum(m)

  # to obtain hierarchical omega
  btw_cov <- m
  for (i in seq_along(grp_start)) { # to conveniently compute the average correlation using na.rm
    btw_cov[grp_start[i]:grp_end[i], grp_start[i]:grp_end[i]] <- NA
  }
  avg_btw_cov <- mean(btw_cov, na.rm = TRUE)
  omega_h <- n^2 * avg_btw_cov / sum(m)

  # output
  out <- list(multidimensional_reliability = multi_rel,
              omega_hiearchical = omega_h,
              subdimensional_reliability = sub_alpha
              )
  if (mp) {
    cat("multidimensional parallel reliability                    ", multi_rel, "\n")
    cat("omega_hierarchical (from multidimensional parallel  )    ", omega_h, "\n")
    cat("Sub-dimensional reliability (standardized alpha)         ", sub_alpha, "\n")
  } else {
    cat("stratified alpha (multidimensional reliability)          ", multi_rel, "\n")
    cat("omega_hierarchical (from multidimensional tau-equivalent)", omega_h, "\n")
    cat("Sub-dimensional reliability (coefficient alpha)          ", sub_alpha, "\n")
  }
  invisible(out)
}
