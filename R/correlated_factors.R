#' Obtain correlated factors reliability estimates of two-dimensional data
#'
#' Obtain correlated factors reliability estimates of two-dimensional data. It is a kind
#' of multidimensional CFA reliability, and it is the reliability derived from
#' the correlated factors model. Items should be grouped by each sub-dimension. For
#' example, the first dimension consists of the first item through the until item.
#'
#' @param x observed item scores or their covariances
#' @param until The number of items up to the first sub-construct
#' @param nonneg_loading if TRUE, constraint loadings to nonnegative values
#' @param print If TRUE, the result is printed to the screen.
#' @return a correlated factors reliability estimate
#' @examples correlated_factors(Osburn_moderate, 4)
#' @examples correlated_factors(Sijtsma2a, c(2, 4))
#' @references Cho, E. (2016). Making reliability reliable: A systematic
#' approach to reliability coefficients. Organizational Research Methods, 19(4),
#' 651-682. https://doi.org/10.1177/1094428116656239
#' @author Eunseong Cho, \email{bene@kw.ac.kr}
#'
correlated_factors <- function(x, until, nonneg_loading = FALSE, print = TRUE) {
  stopifnot(requireNamespace("lavaan"))
  m <- get_cov(x)
  n <- nrow(m)
  colnames(m) <- rownames(m) <- character(length = n)
  grp_start <- c(1, until + 1)
  grp_end <- c(until, n)
  model_str <- vector("character")

  for (i in 1:n) {
    rownames(m)[i] <- paste0("V", i)
  }

  for (i in 1:n) {
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
  } # end of for (i in 1:n)

  for (i in seq_along(grp_start)) {
    model_str <- paste0(model_str, "\n G", i, " ~~ 1 * G", i)
  }

  for (i in 1:n) { # to prevent negative errors
    model_str <- paste0(model_str, "\n V", i, " ~~ e", i, "*V", i, "\n e", i, "> 0.0000001")
    if (!any(i == grp_start) & nonneg_loading) { # to prevent negative loadings
      model_str <- paste0(model_str, "\n b", i, " > .0000001")
    }
  }

  fit <- lavaan::cfa(model_str, sample.cov = m, sample.nobs = 500)
  if (lavaan::inspect(fit, what = "converged")) {
    theta <- lavaan::inspect(fit, what = "est")$theta
    implied <- lavaan::inspect(fit, what = "implied")[[1]]
    multi_rel <- 1 - sum(theta)/sum(implied)
    # fit indices & estimates
    fit_indices <- lavaan::inspect(fit, what = "fit")
    est <- lavaan::inspect(fit, what = "est")
    # output
    out <- list(multidimensional_reliability = multi_rel,
                fit_indices = fit_indices,
                estimates = est)
  } else { # NA for non-convergent solutions
    out <- NA
  }

  if (print) {
    cat("correlated factors reliability (multidimensional CFA)    ", multi_rel, "\n")
  }
  invisible(out)
}
