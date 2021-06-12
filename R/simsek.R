#' Obtain Şimşek-Noyan's theta (principal component analysis reliability)
#'
#' Şimşek-Noyan's (2013) theta is the multidimensional  principal component analysis (PCA)
#' reliability. It is a multidimensional generalization of Kaiser-Caffrey's alpha.
#' @param x a dataframe or a matrix (multidimensional)
#' @param dim the number of dimensions
#' @return Şimşek-Noyan's theta
#' @export Şimşek
#' @examples Şimşek(Osburn_moderate)
#' @references Armor, D. J. (1974). Theta reliability and factor scaling.
#' In H. L. Costner (Ed.), Sociological methodology (pp. 17-50). Jossey-Bass.
#' @references Kaiser, H. F., & Caffrey, J. (1965). Alpha factor analysis.
#' Psychometrika, 30(1), 1-14.
#' @references  Şimşek, G. G., & Noyan, F. (2013). McDonald’s ωt, Cronbach’s α,
#' and Generalized θ for Composite Reliability of Common Factors Structures.
#' Communications in Statistics - Simulation and Computation, 42(9), 2008–2025.
#' https://doi.org/10.1080/03610918.2012.689062
#'
simsek <- function(x, dim) {
  m <- get_cov(x)
  k <- nrow(m)
  sum_eigen <- sum(eigen(stats::cov2cor(m))$values[1:dim])
  out <- k / (k - dim) * (1 - dim / sum_eigen)
  cat("Simsek-Noyan's theta (principal component analysis reliability) ", out,
      "\n")
  invisible(out)
}
