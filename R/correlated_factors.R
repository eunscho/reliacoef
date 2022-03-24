#' Obtain correlated factors reliability estimates of two-dimensional data
#'
#' Obtain correlated factors reliability estimates of two-dimensional data. It is a kind
#' of multidimensional CFA reliability, and it is the reliability derived from
#' the correlated factors model. Items should be grouped by each sub-dimension. For
#' example, the first dimension consists of the first item through the until item.
#'
#' @param x observed item scores or their covariances
#' @param nobs number of observations (i.e., sample size)
#' @param until The number of items up to the first sub-construct
#' @param print If TRUE, the result is printed to the screen.
#' @return a correlated factors reliability estimate
#' @examples correlated_factors(Osburn_moderate, 4)
#' @examples correlated_factors(Sijtsma2a, c(2, 4))
#' @references Cho, E. (2016). Making reliability reliable: A systematic
#' approach to reliability coefficients. Organizational Research Methods, 19(4),
#' 651-682. https://doi.org/10.1177/1094428116656239
#' @author Eunseong Cho, \email{bene@kw.ac.kr}
#'
correlated_factors <- function(x, nobs = NULL, until, print = TRUE) {
  stopifnot(requireNamespace("lavaan"))
  if(nrow(x) == ncol(x)) {
    m <- x
    nobs <- ifelse(is.null(nobs), 200, nobs)
  } else {
    m <- cov(x)
    nobs <- nrow(x)
  }
  colnames(m) <- rownames(m) <- character(length = nrow(m))
  grp_start <- c(1, until + 1)
  grp_end <- c(until, nrow(m))
  model_str <- vector("character")

  for (i in 1:nrow(m)) {
    rownames(m)[i] <- paste0("V", i)
  }

  for (i in 1:nrow(m)) {
    if (any(i == grp_start)) {
      which <- which(i == grp_start)
      for (j in i:grp_end[which]) {
        if (j == i) {
          model_str <- paste0(model_str, "\n G", which, "=~ NA*V", j)
        } else {
          model_str <- paste0(model_str, " + b", j, "*V", j)
        }
      } # end of for (j in i:grp_end[which])
    } # end of if(any(i == grp_start))
  } # end of for (i in 1:nrow(m))

  for (i in seq_along(grp_start)) {
    model_str <- paste0(model_str, "\n G", i, " ~~ 1 * G", i)
  }

  for (i in 1:nrow(m)) { # to prevent negative errors
    model_str <- paste0(model_str, "\n V", i, " ~~ e", i, "*V", i, "\n e", i, "> 0")
  }

  lav_out <- lavaan::cfa(model_str, sample.cov = m, sample.nobs = 500)
  if (lavaan::inspect(lav_out, what = "converged")) {
    fit <- lavaan::inspect(lav_out, what = "fit")
    est <- lavaan::inspect(lav_out, what = "est")
    implied <- lavaan::inspect(lav_out, what = "implied")[[1]]
    rel <- 1 - sum(est$theta)/sum(implied)
    # to obtain subdimensional reliability
    sub_rel <- vector("double", length = length(grp_start))
    for (i in seq_along(grp_start)) { # to conveniently compute the average correlation using na.rm
      sub_uniq <- sum(est$theta[grp_start[i]:grp_end[i], grp_start[i]:grp_end[i]])
      sub_implied <- sum(implied[grp_start[i]:grp_end[i], grp_start[i]:grp_end[i]])
      sub_rel[i] <- 1 - sub_uniq / sub_implied
    }
    # output
    out <- list(rel = rel,
                sub_rel = sub_rel,
                fit = fit,
                est = est)
  } else { # NA for non-convergent solutions
    out <- NA
  }

  if (print) {
    cat("correlated factors reliability (multidimensional CFA)    ", rel, "\n")
    cat("Sub-dimensional reliability (congeneric reliability)     ", sub_rel, "\n")
  }
  invisible(out)
}
