#' Obtain bifactor reliability estimates of two-dimensional data
#'
#' Obtain bifactor reliability estimates of two-dimensional data. It is a kind
#' of multidimensional CFA reliability, and it is the reliability derived from
#' the bifactor model. Items should be grouped by each sub-dimension. For
#' example, the first dimension consists of the first item through the until item.
#'
#' @param x observed item scores or their covariances
#' @param until The number of items up to the first sub-construct
#' @param only_first_grp whether it is only_first_grp and Cui's simulation, which has only one group factor
#' @param nonneg_loading if TRUE, constraint loadings to nonnegative values
#' @return a bifactor reliability estimate
#' @examples bifactor(Osburn_moderate, 4)
#' @references Cho, E. (2016). Making reliability reliable: A systematic
#' approach to reliability coefficients. Organizational Research Methods, 19(4),
#' 651-682. https://doi.org/10.1177/1094428116656239
#' @author Eunseong Cho, \email{bene@kw.ac.kr}
#'
bifactor <- function(x, until, only_first_grp = FALSE, nonneg_loading = FALSE) {
  stopifnot(requireNamespace("lavaan"))
  m <- get_cov(x)
  n <- nrow(m)
  grp_start <- c(1, until + 1)
  grp_end <- c(until, n)
  colnames(m) <- rownames(m) <- character(length = n)

  for (i in 1:n) {
    rownames(m)[i] <- paste0("V", i)
    # general factor
    if (i == 1) {
      model_str <- paste("F =~ NA*V1")
    } else {
      model_str <- paste0(model_str, "+ a", i, "*V", i)
    }
  }
  if (only_first_grp) { # when only the first group exists (Tang and Cui, 2012)
    for (i in 1:until[1]) {
      if (i == 1) {
        model_str <- paste0(model_str, "\n G1 =~ NA*V1")
      } else {
        model_str <- paste0(model_str, " + b", i, "*V", i)
      }
    }
    model_str <- paste0(model_str, "\n F ~~ 1*F \n G1 ~~ 1*G1")
    model_str <- paste0(model_str, "\n F ~~ 0*G1")
  } else {# group factors in common cases
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
    model_str <- paste0(model_str, "\n F ~~ 1*F")
    for (i in seq_along(grp_start)) {
      model_str <- paste0(model_str, "\n G", i, " ~~ 1*G", i)
      model_str <- paste0(model_str, "\n F ~~ 0 * G", i)
      if (i < length(grp_start)) {
        for (j in (i + 1):length(grp_start)) {
          model_str <- paste0(model_str, "\n G", i, "~~ 0 * G", j)
        } # end of for (j in (i + 1):length(grp_start))
      } # end of if (i < length(grp_start))
    } # end of for (i in seq_along(grp_start))
  } # end of else


  if (nonneg_loading) { # prevent negative loadings
    for (i in 2:n) {
      model_str <- paste0(model_str, "\n a", i, " > .0000001")
      if (!only_first_grp & !any(i == grp_start)) {
        model_str <- paste0(model_str, "\n b", i, " > .0000001")
      }
    }
  }

  for (i in 1:n) { # to prevent negative errors
    model_str <- paste0(model_str, "\n V", i, " ~~ e", i, "*V", i, "\n e", i, "> 0.0000001")
  }

  fit <- lavaan::cfa(model_str, sample.cov = m, sample.nobs = 500)

  if (lavaan::inspect(fit, what  = "converged")) {
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
    sum_gen_loading <- sum(lavaan::inspect(fit, what = "est")$lambda[, 1])
    omega_h <- sum_gen_loading^2 / sum(implied)

    # fit indices & estimates
    fit_indices <- lavaan::inspect(fit, what = "fit")
    est <- lavaan::inspect(fit, what = "est")

    # output
    out <- list(multidimensional_reliability = multi_rel,
                omega_hierarchical = omega_h,
                subdimensional_reliability = subdim_rel,
                fit_indices = fit_indices,
                estimates = est)
  } else {# NA for non-convergent solutions
    out <- NA
  }
  cat("bifactor reliability (multidimensional CFA)              ", multi_rel, "\n")
  cat("omega_hierarchical (from bifactor model)                 ", omega_h, "\n")
  cat("Sub-dimensional reliability (congeneric reliability)     ", subdim_rel, "\n")
  invisible(out)
}
