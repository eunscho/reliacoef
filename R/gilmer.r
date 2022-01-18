#' Obtain the Gilmer-Feldt reliability coefficient
#'
#' It is a unidimensional reliability coefficient
#' based on a congeneric model. The congeneric model is a model that allows
#' the length, discrimination, or importance of items to be different,
#' and is the least restrictive model among the models derived from
#' the classical test theory. The Gilmer-Feldt coefficient has the advantage
#' of being less computational than congeneric reliability (Joreskog 1971)
#' which uses confirmatory factor analysis.However, the Gilmer-Feldt coefficient
#'  derives a value very close to congeneric reliability (Cho in press).
#'  Feldt and Charter (2003) offers a user-friendly review of the Gilmer-Feldt
#'  coefficient.
#'
#' @param x A dataframe or a matrix (unidimensional)
#' @return The Gilmer-Fedlt coefficient
#' @export gilmer
#' @examples gilmer(Graham1)
#' @references Cho, E. (in press). Neither Cronbach's alpha nor McDonald's omega: A comment on Sijtsma and Pfadt. Psychometrika.
#' @references Feldt, L. S., & Charter, R. A. (2003). Estimation of internal consistency reliability when test parts vary in effective length. Measurement and Evaluation in Counseling and Development, 36(1), 23-27
#' @references Gilmer, J. S., & Feldt, L. S. (1983). Reliability estimation for a test with parts of unknown lengths. Psychometrika, 48(1), 99–111.
#' @references Jöreskog, K. G. (1971). Statistical analysis of sets of congeneric tests. Psychometrika, 36(2), 109–133.
gilmer <- function(x) {
  m <- get_cov(x)
  k <- nrow(m)
  total <- sum(m)
  D <- numeric(k)
  nondiag <- rowSums(m) - diag(m)
  key_row <- m[order(nondiag)[k], ] # A, B, C, D in Feldt and Charter
  for (i in 1:k) {
      if (nondiag[i] == max(nondiag)) {
        D[i] <- 1
      } else {
        D[i] <- (nondiag[i] - key_row[i]) / (max(nondiag) - key_row[i])
      }
  }
  W <- cumsum(D^2)[k]
  Q <- sum(D)^2
  out <- (Q / (Q - W)) * (sum(nondiag) / total)
  cat("Gilmer-Feldt's reliability coefficient                            ", out,
      "\n")
  invisible(out)
}
