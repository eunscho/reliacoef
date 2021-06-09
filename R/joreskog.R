#' Obtain Jöreskog's congeneric reliability (Unidimensional CFA reliability)
#'
#' Congeneric reliability is a reliability coefficient derived from
#' unidimensional confirmatory factor analysis (CFA).
#'
#'  Features: Congeneric reliability is a unidimensional reliability
#' coefficient based on a unidimensional confirmatory factor analysis (CFA)
#' model.
#'
#' Name: Congeneric reliability is called by a variety of names, general users
#' usually call it composite reliability, and reliability researchers often call
#' it omega. One of the reasons for this confusion is that studies that first
#' proposed this coefficient (Jöreskog 1971) did not give this formula a name
#'  (Cho 2016). Jöreskog (1971) proposed a matrix-form formula, and the
#'  commonly known non-matrix formula appears in Werts et al. (1974).
#'
#' Frequency of use: Congeneric reliability is the second most commonly used
#' reliability coefficient after coefficient alpha (Cho 2016)
#'
#' Accuracy: Congeneric reliability is the most accurate reliability coefficient
#' along with the Feldt-Gilmer coefficient (Cho in press)
#'
#' Computation: This function uses maximum likelihood as estimation,
#' unstandardized covariance matrix as input, and lavaan package as software.
#'
#' @param x a dataframe or a matrix (unidimensional)
#' @param nonneg_loading if TRUE, constraint loadings to nonnegative values
#' @return congeneric reliability coefficient
#' @export joreskog
#' @examples joreskog(Graham1)
#' @references Cho, E. (2016). Making reliability reliable: A systematic
#' approach to reliability coefficients. Organizational Research Methods, 19(4),
#' 651-682.
#' @references Cho, E. (in press). Neither Cronbach's alpha nor McDonald's
#' omega: A comment on Sijtsma and Pfadt. Psychometrika.
#' @references Jöreskog, K. G. (1971). Statistical analysis of sets of
#' congeneric tests. Psychometrika, 36(2), 109-133.
#' @references Werts, C. E., Linn, R. L., & Jöreskog, K. G. (1974). Intraclass
#' reliability estimates: Testing structural assumptions. Educational and
#' Psychological Measurement, 34, 25-33.
#' @seealso [gilmer()] for the Gilmer-Feldt coefficient
#' @seealso [feldt()] for classical congeneric reliability coefficient
#' @seealso [psych::omega()] for a related function of the package psych
#' @seealso [MBESS::ci.reliability()] for a related function of the package MBESS
#' @seealso [Lambda4::omega.tot()] for a related function of the package Lambda4
#' @import matrixcalc
#'
joreskog <- function(x, nonneg_loading = FALSE) {
    stopifnot(requireNamespace("matrixcalc"))
    cov <- get_cov(x)
    if (!matrixcalc::is.positive.definite(cov)) {
      out <- NA
    } else {
      est <- uni_cfa(cov, nonneg_loading = nonneg_loading)
      if (any(is.na(est))) {
        out <- NA
      } else {
        sum_lambda <- sum(est$lambda)
        sum_theta <- sum(est$theta)
        out <- sum_lambda^2/(sum_lambda^2 + sum_theta)
      }
    }
    cat("Joreskog's congeneric (unidimensional CFA) coefficient            ", out,
        "\n")
    invisible(out)
}
