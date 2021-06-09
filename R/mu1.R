#' Obtain Ten Berge and Zegers' (1978) mu1
#'
#' Obtain Ten Berge and Zegers' (1978) mu1. mu1 equals Guttman's (1945) lambda2.
#' @author Eunseong Cho, \email{bene@kw.ac.kr}
#' @param x a dataframe or a matrix (unidimensional)
#' @references Guttman, L. (1945). A basis for analyzing test-retest reliability.
#'  Psychometrika, 10(4), 255-282.
#' @references Ten Berge, J. M. F., & Zegers, F. E. (1978). A series of lower
#' bounds to the reliability of a test. Psychometrika, 43(4), 575-579.
mu1 <- function(x) {
    m <- get_cov(x)
    n <- nrow(m)/(nrow(m) - 1)
    off <- m
    diag(off) <- 0
    numerator <- sum(off) + sqrt(n * sum(off^2))
    out <- numerator/sum(m)
    cat("Guttman's lambda2 (mu1)                                           ", out,
        "\n")
    invisible(out)
}
