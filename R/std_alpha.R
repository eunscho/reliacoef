#' Obtain standardized alpha
#'
#' @param x a dataframe or a matrix (unidimensional)
#' @return standardized alpha reliability estimate
#' @export std_alpha
#' @examples std_alpha(Graham1)
#'
std_alpha <- function(x) {
  m <- get_cov(x, cor = TRUE)
  n <- nrow(m)/(nrow(m) - 1)
  off <- m
  diag(off) <- 0
  out <- n * sum(off)/sum(m)
  cat("standardized alpha                                                ", out,
      "\n")
  invisible(out)
}
