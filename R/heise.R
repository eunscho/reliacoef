#' Obtain Heise-Bohrnstedt's Omega (another factor analysis reliability)
#'
#' Heise-Bohrnstedt's (1970) Omega is an factor analysis (FA) reliability.
#' This formula is different from the FA reliability we use today and yields a
#' larger value (Cho, in press). McDonald (1999) referred to all FA reliability
#' as omega, and capitalized omega was used to distinguish it from McDonald's
#' omega.
#' @param x a dataframe or a matrix (unidimensional)
#' @return Heise-Bohrnstedt's Omega
#' @import psych
#' @export heise
#' @examples heise(Graham1)
#' @references Cho, E. (in press). Neither Cronbach's alpha nor McDonald's
#' omega: A comment on Sijtsma and Pfadt. Psychometrika.
#' @references Heise, D. R., & Bohrnstedt, G. W. (1970). Validity, invalidity,
#' and reliability. Sociological Methodology, 2, 104-129.
#' @references McDonald, R. P. (1999). Test theory: A unified treatment.
#' Lawrence Erlbaum.
#'
heise <- function(x) {
  stopifnot(requireNamespace("psych"))
  m <- get_cov(x)
  k <- nrow(m)
  sum <- cumsum(diag(m) * psych::fa(stats::cov2cor(m), nfactors = k)$communality)[k]
  out <- 1 - (sum(diag(m)) - sum) / sum(m)
  cat("Heise-Borhnstedt's Omega                                          ", out,
      "\n")
  invisible(out)
}
