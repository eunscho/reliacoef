#' Obtain Hancock's H (CFA version of maximal reliability)
#'
#' It is the confirmatory factor analysis (CFA) version of maximal
#' reliability. This coefficient takes the standardized factor loading as the
#' reliability of each item, and finds the weight that maximizes the
#' reliability. Hence, Hancock's H shows a different result than the reliability
#' estimator using conventional unit weights.
#' @param x a dataframe or a cov (unidimensional)
#' @param nonneg_loading if TRUE, constraint loadings to nonnegative values
#' @return Hancock's H
#' @export hancock
#' @import maxtrixcalc
#' @examples hancock(Graham1)
#' @references Cho, E. (in press). Neither Cronbach's alpha nor McDonald's
#' omega: A comment on Sijtsma and Pfadt. Psychometrika.
#' @references Hancock, G., & Mueller, R. O. (2001). Rethinking construct
#' reliability within latent variable systems. In R. Cudeck, S. du Toit, & D.
#' Sörbom (Eds.), Structural equation modeling: Present and future-A festschrift
#' in honor of Karl Jöreskog (pp. 195-216). Scientific Software International.
#' @references Li, H., Rosenthal, R., & Rubin, D. B. (1996). Reliability of
#' measurement in psychology: From Spearman-Brown to maximal reliability.
#' Psychological Methods, 1(1), 98-107.
#' @references McNeish, D. (2017). Thanks coefficient alpha, we’ll take it from
#' here. Psychological Methods, 23(3), 412-433.
hancock <- function(x, nonneg_loading = FALSE) {
    stopifnot(requireNamespace("matrixcalc"))
    cov <- get_cov(x)
    if (!matrixcalc::is.positive.definite(cov)) {
        out <- NA
    } else {
      est <- uni_cfa(cov, what = "std", nonneg_loading = nonneg_loading)
      if (any(is.na(est))) {
        out <- NA
      } else {
        prop_lambda <- est$lambda^2 / (1 - est$lambda^2)
        out <- 1 / (1 + 1 / sum(prop_lambda))
      }
    }
    cat("Hancock's H (unidimensional CFA version of maximal reliability)   ", out,
        "\n")
    invisible(out)
}
