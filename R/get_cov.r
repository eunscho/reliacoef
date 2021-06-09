#' Obtain the covariance matrix
#'
#' If the input data is a square matrix, it is converted into a matrix,
#' otherwise the covariance matrix is obtained.
#'
#' @param x A dataframe or a matrix
#' @param cor if TRUE, return correlation matrix. if FALSE, return covariance matrix
#' @return The covariance or correlation matrix
get_cov <- function(x, cor = FALSE) {
    if (nrow(x) == ncol(x)) {
        if (cor) {
            out <- stats::cov2cor(as.matrix(x))
        } else {
            out <- as.matrix(x)
        }
    } else {
        if (cor) {
            out <- stats::cor(x, use = "pairwise.complete.obs")
        } else {
            out <- stats::cov(x, use = "pairwise.complete.obs")
        }
    }
    return(out)
}
