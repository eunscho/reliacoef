#' Obtain Ten Berge and Zegers' (1978) mu2
#'
#' Obtain Ten Berge and Zegers' (1978) mu2.
#' @author Eunseong Cho, \email{bene@kw.ac.kr}
#' @param x a dataframe or a matrix (unidimensional)
#' @export mu2
#' @references Ten Berge, J. M. F., & Zegers, F. E. (1978). A series of lower
#' bounds to the reliability of a test. Psychometrika, 43(4), 575-579.
#' @examples mu2(Graham1)
mu2 <- function(x) {
    m <- get_cov(x)
    n <- nrow(m)/(nrow(m) - 1)
    off <- m
    diag(off) <- 0
    numerator <- sum(off) + sqrt(sum(off^2) + sqrt(n * sum(off^4)))
    out <- numerator/sum(m)
    cat("Ten Berge and Zegers' mu2                                         ", out,
        "\n")
    invisible(out)
}
