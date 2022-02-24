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
#' @param print If TRUE, the result is printed to the screen.
#' @return multidimensional_reliability a stratified alpha reliability estimate
#' @return subdimensional_reliability reliability estimates of subdimensions
#' @return omega hierarchical general factor saturation
#' @examples stratified_alpha(Osburn_moderate, 4)
#' @references Rajaratnam, N., Cronbach, L. J., & Gleser, G. C. (1965).
#' Generalizability of stratified-parallel tests. Psychometrika, 30(1), 39–56.
#' https://doi.org/10.1007/BF02289746
#' @references Cho, E. (2016). Making reliability reliable: A systematic approach
#' to reliability coefficients. Organizational Research Methods, 19(4), 651-682.
#' https://doi.org/10.1177/1094428116656239
#' @author Eunseong Cho, \email{bene@kw.ac.kr}
#'
stratified_alpha <- function(x, until, mp = FALSE, print = TRUE) {
  m <- get_cov(x)
  grp_start <- c(1, until + 1)
  grp_end <- c(until, nrow(m))
  item <- round(nrow(m)/length(grp_start))
  
  # to obtain hierarchical omega
  btw_cov <- m
  for (i in seq_along(grp_start)) { # to conveniently compute the average correlation using na.rm
    btw_cov[grp_start[i]:grp_end[i], grp_start[i]:grp_end[i]] <- NA
  }
  avg_btw_cov <- mean(btw_cov, na.rm = TRUE)
  omegah <- nrow(m)^2 * avg_btw_cov / sum(m)
  
  sub_rel <- sub_omegah <- sub_var <- vector("double", length(grp_start))
  for (i in seq_along(grp_start)) {
    sub_mat <- m[grp_start[i]:grp_end[i], grp_start[i]:grp_end[i]]
    sub_rel[i] <- reliacoef::alpha(sub_mat, print = FALSE)
    sub_var[i] <- sum(sub_mat)
    sub_omegah[i] <- item^2 * avg_btw_cov / sub_var[i]
  }

  sum <- vector("double", length = 1)
  for (i in seq_along(grp_start)) {
     sum <- sum + sub_var[i] * (1 - sub_rel[i])
  }
  rel <- 1 - sum / sum(m)

  # output
  out <- list(rel = rel,
              omegah = omegah,
              sub_rel = sub_rel,
              sub_omegah = sub_omegah)
  
  if (print) {
    if (mp) {
    cat("multidimensional parallel reliability                    ", rel, "\n")
    } else {
    cat("stratified alpha (multidimensional reliability)          ", rel, "\n")
    }
    cat("omegahierarchical (from multidimensional tau-equivalent) ", omegah, "\n")
    cat("Sub-dimensional reliability (coefficient alpha)          ", sub_rel, "\n")
    cat("Sub-dimensional omega hiearchical                        ", sub_omegah, "\n")
   }

  invisible(out)
}
