#' Unidimensional confirmatory factor analysis
#'
#' @param cov observed covariances
#' @param what e.g., "est", "std", "fit"
#' @param sample_size number of sample observations
#' @param nonneg_loading if TRUE, constraint loadings to nonnegative values
#' @param nonneg_error if TRUE, constraint loadings to positive values
#' @param taueq if TRUE, a tau-equivalent model is estimated
#' @param parallel if TRUE, a parallel model is estimated
#' @examples uni_cfa(Graham1)
#' @import lavaan
#' @export uni_cfa
#' @return parameter estimates of unidimensional cfa model
uni_cfa <- function(cov, what = "est", sample_size = 500, nonneg_loading = FALSE,
                    nonneg_error = TRUE, taueq = FALSE, parallel = FALSE) {
  stopifnot(requireNamespace("lavaan"))
  k <- nrow(cov)
  rownames(cov) <- character(length = k)
  for (i in 1:k) {
    rownames(cov)[i] <- paste0("V", i)
    if (i == 1) {
      model_str <- paste("F =~ NA*V1")
    } else if (taueq | parallel) { # tau-equivalent or parallel
      model_str <- paste0(model_str, " + equal('F=~V1')*V", i)
    } else {# congeneric
      model_str <- paste0(model_str, " + l", i, "*V", i)
      }
  }
  colnames(cov) <- rownames(cov)
  model_str <- paste0(model_str, " \n F ~~ 1*F", collapse = "\n")
  if (parallel) {
    for (i in 1:k) { # all errors are constained to be equal
      model_str <- paste0(model_str, "\n V", i, " ~~ e*V", i)
    }
  } else if (!taueq) { # congeneric
    for (i in 1:k) { # to prevent negative errors
      if (nonneg_error) {
        model_str <- paste0(model_str, "\n V", i, " ~~ e", i, "*V", i, "\n e", i,
                            "> 0")
      }
      if (i > 1 & nonneg_loading) { # prevent negative loadings
        model_str <- paste0(model_str, "\n l", i, "> .0")
      }
    }
  }
  fit <- lavaan::cfa(model_str, sample.cov = cov, sample.nobs = sample_size)
  if (lavaan::inspect(fit, what = "converged")) {
    out <- lavaan::inspect(fit, what = what)
  } else {
    out <- NA
  }
  return(out)
}
