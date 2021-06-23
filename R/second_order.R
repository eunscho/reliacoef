#' Obtain second-order factor reliability estimates of two-dimensional data
#'
#' Obtain second-order factor reliability estimates of two-dimensional data. It
#' is a kind of multidimensional CFA reliability, and it is the reliability
#' derived from the second-order factor model. Items should be grouped by each
#' sub-dimension. For example, the first dimension consists of the first item
#' through the until item.
#'
#' @param x observed item scores or their covariances
#' @param until The number of items up to the first sub-construct
#' @param only_first_grp TRUE if only_first_grp. In this case, this function does not compute a
#' reliability estimate.
#' @param nonneg_loading if TRUE, constraint loadings to nonnegative values
#' @param print If TRUE, the result is printed to the screen.
#' @return a second-order factor reliability estimate
#' @examples second_order(Osburn_moderate, 4)
#' @examples second_order(Sijtsma2a, c(2, 4))
#' @references Cho, E. (2016). Making reliability reliable: A systematic
#' approach to reliability coefficients. Organizational Research Methods, 19(4),
#' 651-682. https://doi.org/10.1177/1094428116656239
#' @author Eunseong Cho, \email{bene@kw.ac.kr}
#'
second_order <- function(x, until, only_first_grp = FALSE,
                         nonneg_loading = FALSE, print = TRUE) {
  stopifnot(requireNamespace("lavaan"))
  if (only_first_grp == TRUE) {
    out <- NA
  } else {
    m <- get_cov(x)
    n <- nrow(m)
    grp_start <- c(1, until + 1)
    grp_end <- c(until, n)
    colnames(m) <- rownames(m) <- character(length = n)
    model_str <- vector("character")

    for (i in 1:n) {
      rownames(m)[i] <- paste0("V", i)
      if (any(i == grp_start)) {
        which <- which(i == grp_start)
        for (j in i:grp_end[which]) {
          if (j == i) {
            model_str <- paste0(model_str, "\n FoF", which, "=~ NA*V", j)
          } else {
            model_str <- paste0(model_str, " + b", j, "*V", j)
          }
        } # end of for (j in i:grp_end[which])
      } # end of if(any(i == grp_start))
    } #end of for (i in 1:n)
    colnames(m) <- rownames(m)

    for (i in seq_along(grp_start)) {
      if (i == 1) {
        model_str <- paste0(model_str, "\n SoF =~ NA * FoF1")
      } else {
        model_str <- paste0(model_str, " + a", i, " * FoF", i)
      }
    }

    for (i in 2:length(grp_start)) {
      model_str <- paste0(model_str, "\n a", i, " > .0000001")
      model_str <- paste0(model_str, "\n a", i, " < .9999999")
      model_str <- paste0(model_str, "\n FoF", i, " ~~ 1 * FoF", i)
    }

    model_str <- paste0(model_str, "\n SoF ~~ 1 * SoF")

    for (i in 1:n) { # to prevent negative errors
      model_str <- paste0(model_str, "\n V", i, " ~~ e", i, "*V", i, "\n e", i, "> 0.0000001")
      if (!any(i == grp_start) & nonneg_loading) { # to prevent negative loadings
        model_str <- paste0(model_str, "\n b", i, " > .0000001")
      }
    }

    fit <- lavaan::cfa(model_str, sample.cov = m, sample.nobs = 500)
    if (lavaan::inspect(fit, what = "converged")) {
      # to obtain multidimensional reliability
      theta <- lavaan::inspect(fit, what = "est")$theta
      implied <- lavaan::inspect(fit, what = "implied")[[1]]
      multi_rel <- 1 - sum(theta)/sum(implied)

      # to obtain subdimensional reliability
      subdim_rel <- vector("double", length = length(grp_start))
      for (i in seq_along(grp_start)) { # to conveniently compute the average correlation using na.rm
        sum_uniq <- sum(theta[grp_start[i]:grp_end[i], grp_start[i]:grp_end[i]])
        sum_implied <- sum(implied[grp_start[i]:grp_end[i], grp_start[i]:grp_end[i]])
        subdim_rel[i] <- 1 - sum_uniq / sum_implied
      }

      # to obtain hierarchical omega
      lambda <- lavaan::inspect(fit, what = "est")$lambda[, 1:2]
      beta <- lavaan::inspect(fit, what = "est")$beta[1:2, 3]
      omega_h <- sum(lambda %*% beta) ^ 2 / sum(implied)

      # fit indices & estimates
      fit_indices <- lavaan::inspect(fit, what = "fit")
      est <- lavaan::inspect(fit, what = "est")

      # output
      out <- list(multidimensional_reliability = multi_rel,
                  subdimensional_reliability = subdim_rel,
                  hierarchica_omega = omega_h,
                  fit_indices = fit_indices,
                  estimates = est)

    } else { # NA for non-convergent solutions
      out <- NA
    }
  }

  if (print) {
    cat("second-order factor reliability (multidimensional CFA)   ", multi_rel, "\n")
    cat("omega_hierarchical (from second-order factor model)      ", omega_h, "\n")
    cat("Sub-dimensional reliability (congeneric reliability)     ", subdim_rel, "\n")
  }
  invisible(out)
}
