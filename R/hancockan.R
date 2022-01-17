#' Obtain Ten Berge and Zegers (2004) mu3
#'
#' Obtain Ten Berge and Zegers (2004) mu3.
#'
#' Hancock and An (2020) published a reliability estimator they called the
#' "closed-form omega", which is designed to approximate a unidimensional
#' confirmatory factor analysis reliability estimator (i.e., joreskog).
#'
#' @author Eunseong Cho, \email{bene@kw.ac.kr}
#' @param x a dataframe or a matrix (unidimensional)
#' @export hancockan
#' @import utils
#' @references Hancock, G. R., & An, J. (2020). A Closed-Form Alternative for
#' Estimating ω Reliability under Unidimensionality. Measurement:
#' Interdisciplinary Research and Perspectives, 18(1), 1-14.
#' https://doi.org/10.1080/15366367.2019.1656049
#' @examples hancockan(Graham1)

hancockan <- function(x) {
  m <- get_cov(x)
  n.item <- ncol(m)
  n.lamb = (n.item - 1) * (n.item - 2)/2
  #~~~ place holders
  lamb.sq.num = matrix(0,n.item, n.lamb)
  lamb.sq.denom = matrix(0,n.item, n.lamb)
  #~~~ calculate loadings based on the covariance of all the other items
  for (i in 1:n.item) {
    A = utils::combn(setdiff((1:n.item), i), 2)
    for (b in seq(1, (n.lamb * 2 - 1), by = 2)) {
      lamb.sq.num[i, (b + 1)/2] = m[i, A[b]] * m[i, A[b + 1]]
      lamb.sq.denom[i, (b + 1)/2] = m[A[b], A[b + 1]]
    }
  }
  lamb = sqrt(rowSums(lamb.sq.num)/rowSums(lamb.sq.denom))
  tot.var = sum(m)
  omega.cal = (sum(lamb)^2)/tot.var

  return(omega.cal)
}
