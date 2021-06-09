#' Obtain Ten Berge and Zegers (1978) mu4
#'
#' Obtain Ten Berge and Zegers (1978) mu4.
#'
#' The original formula and the formula of psych's tenberge() are different.
#' There is a high possibility that the original formula is incorrect and psych's
#' version is correct. According to Equation (4) of the original article, mu
#' should increase monotonically (e.g., mu4>=mu3), but if the original formula is
#' followed, it may decrease in some cases. The formula of the original paper is
#' 2h, but changing it to 2^h solves this problem. This function follow the latter
#' interpretation.
#' @author Eunseong Cho, \email{bene@kw.ac.kr}
#' @param x a dataframe or a matrix (unidimensional)
#' @export mu4
#' @references Ten Berge, J. M. F., & Zegers, F. E. (1978). A series of lower
#' bounds to the reliability of a test. Psychometrika, 43(4), 575-579.
#' @examples mu4(Graham1)

mu4 <- function(x) {
    m <- get_cov(x)
    n <- nrow(m)/(nrow(m) - 1)
    off <- m
    diag(off) <- 0
    numerator <- sum(off) + sqrt(sum(off^2) +
                                     sqrt(sum(off^4) +
                                              sqrt(sum(off^8) +
                                                       sqrt(n * sum(off^16)))))
    out <- numerator/sum(m)
    cat("Ten Berge and Zegers' mu4                                         ", out,
        "\n")
    invisible(out)
}
