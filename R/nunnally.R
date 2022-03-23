#' Obtain bottom-up approach multidimensional reliability estimates
#'
#' Among several approaches to estimating multidimensional reliability, these
#' estimators use a bottom-up approach. That is, the test score is divided into 
#' sub-dimensional or sub-test scores. Multidimensional reliability is 
#' obtained by estimating the reliability of each subtest score and combining 
#' them. Different estimates can be obtained depending on how each subtest 
#' reliability is estimated. These estimators use the general formula first 
#' proposed by Jum Nunnally.
#'
#' @param x observed item scores or their covariances
#' @param until The number of items up to the first sub-construct
#' @param method There are three options. joreskog, mu, kaisercaffrey. mu uses mu4.
#' @param print If TRUE, the result is printed to the screen.
#' @return a bifactor reliability estimate
#' @examples nunnally(Osburn_moderate, 4, method = "mu")
#' @examples nunnally(Sijtsma2a, c(2, 4), method = "kaisercaffrey")
#' @references Nunnally, J. C., & Bernstein, I. H. (1994). Psychometric theory 
#' (3rd ed). McGraw-Hill.
#' @author Eunseong Cho, \email{bene@kw.ac.kr}
#'
nunnally <- function(x, until, method = "joreskog", print = TRUE) {
  m <- get_cov(x)
  n <- ncol(x)
  dim <- length(until) + 1
  grp_start <- c(1, until + 1)
  grp_end <- c(until, n)
  subdim_rel <- subprod <- vector("double", dim)
  for (i in 1:dim) {
    subvar <- m[grp_start[i]:grp_end[i], grp_start[i]:grp_end[i]]
    if(method == "joreskog") {
      subdim_rel[i] <- joreskog(subvar, print = F)
    } else if (method == "mu") {
      subdim_rel[i] <- mu4(subvar, print = F)
    } else if (method == "kaisercaffrey") {
      subdim_rel[i] <- kaisercaffrey(subvar, print = F)
    }
    subprod[i] <- sum(subvar) * (1 - subdim_rel[i])
  }
  multi_rel <- 1 - sum(subprod) / sum(m)
  out <- list(multi_rel = multi_rel, subdim_rel = subdim_rel)
  if (print) {
    cat("Multidimensional reliability using the bottom-up approach", multi_rel, "\n")
    cat("Sub-dimensional reliability                              ", subdim_rel, "\n")
  }
  invisible(out)
}