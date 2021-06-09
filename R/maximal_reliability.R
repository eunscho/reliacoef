#' Obtain maximal reliability estimates of two-dimensional data
#'
#' Obtain maximal reliability estimates of two-dimensional data. Items should be
#' grouped by each sub-dimension. For example, the first dimension consists of
#' the first item through the until item.
#'
#' @param x observed item scores or their covariances
#' @param until The number of items up to the first sub-construct
#' @return a maximal reliability estimate
#' @examples maximal_reliability(Osburn_moderate, 4)
#' @references Li, H., Rosenthal, R., & Rubin, D. B. (1996). Reliability of
#' measurement in psychology: From Spearman-Brown to maximal reliability.
#' Psychological Methods, 1(1), 98-107. https://doi.org/10.1037/1082-989X.1.1.98
#' @author Eunseong Cho, \email{bene@kw.ac.kr}
#'
maximal_reliability <- function(x, until) {
  m <- get_cov(x)
  n <- nrow(m)
  grp_start <- c(1, until + 1)
  grp_end <- c(until, n)
  k <- length(grp_start)
  cor <- stats::cov2cor(m)
  diag(cor) <- NA # to conveniently compute the average correlation using na.rm
  r <- vector("double", length = k) # within subdimension average correlation

  for (i in 1:k) {
    r[i] <- mean(cor[grp_start[i]:grp_end[i], grp_start[i]:grp_end[i]], na.rm = TRUE)
  }

  for (i in 1:k) { # to conveniently compute the average correlation using na.rm
    cor[grp_start[i]:grp_end[i], grp_start[i]:grp_end[i]] <- NA
  }
  rho <- mean(cor, na.rm = TRUE) / sqrt(cumprod(r)[length(r)]) # notation of Li et al. (1996)
  delta <- k / (1 + (k - 1) * rho)
  sum <- vector("double", length = 1)
  for (i in 1:k) {
    sum <- sum + (grp_end[i] - grp_start[i] + 1) * r[i] / (1 - r[i])
  }

  out <- sum / (delta + sum)
  cat("maximal reliability (multidimensional)                   ", out, "\n")
  invisible(out)
}
