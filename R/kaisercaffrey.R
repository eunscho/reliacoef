#' Obtain Kaiser-Caffrey's alpha (principal component analysis reliability)
#'
#' Kaiser-Caffrey's (1965) alpha is the principal component analysis (PCA)
#' reliability. They presented this formula in the context of factor analysis,
#' but Bentler (1968) showed that it was in fact PCA reliability. Armor (1974),
#' citing Bentler (1968), referred to this formula as theta, and some studies
#' refer to it as Armor's theta. Kaiser and Caffrey (1965) labeled this formula
#' alpha, and people may have mistaken it for coefficient alpha. See Vehkalahti
#' (2000) and Cho(in press) for further explanation of this formula.
#' @param x a dataframe or a matrix (unidimensional)
#' @return Kaiser-Caffrey's alpha
#' @export kaisercaffrey
#' @examples kaisercaffrey(Graham1)
#' @references Armor, D. J. (1974). Theta reliability and factor scaling.
#' In H. L. Costner (Ed.), Sociological methodology (pp. 17-50). Jossey-Bass.
#' @references Bentler, P. M. (1968). Alpha-maximized factor analysis (alphamax)
#' : Its relation to alpha and canonical factor analysis. Psychometrika, 33(3),
#' 335-345.
#' @references Cho, E. (in press). Neither Cronbach's alpha nor McDonald's
#' omega: A comment on Sijtsma and Pfadt. Psychometrika.
#' @references Kaiser, H. F., & Caffrey, J. (1965). Alpha factor analysis.
#' Psychometrika, 30(1), 1-14.
#' @references Vehkalahti, K. (2000). Reliability of measurement scales:
#' Tarkkonen's general method supersedes Cronbach's alpha. University of
#' Helsinki.
#'
kaisercaffrey <- function(x) {
  matrix <- get_cov(x)
  k <- nrow(matrix)
  first_eigen <- eigen(stats::cov2cor(matrix))$values[1]
  out <- k / (k - 1) * (1 - 1 / first_eigen)
  cat("Kaiser-Caffrey's alpha (principal component analysis reliability) ", out,
      "\n")
  invisible(out)
}
