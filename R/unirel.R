#' Obtain various unidimensional reliablity coefficients
#'
#' R packages that provide unidimensional reliability coefficients include
#' psych. These are usually separate functions, so comparing
#' multiple reliability estimates can take a long time. This function shows the
#'  reliability coefficients newly added in this package together with the
#'  reliability coefficients provided by psych. This function can
#'  provide Table 2 of Cho(in press) immediately.
#'
#'@author Eunseong Cho, \email{bene@kw.ac.kr}
#'@param x a dataframe or a matrix (unidimensional)
#'@param psych.include Whether to include reliability coefficients
#'(GLB.algebraic, GLB.fa) provided by the package psych
#'@export unirel
#' @references Cho, E. (in press). Neither Cronbach's alpha nor McDonald's
#' omega: A comment on Sijtsma and Pfadt. Psychometrika.
#' @examples unirel(Graham1)
#' @seealso [psych::GLB.algebraic()] for a related function of the package psych
#' @seealso [psych::GLB.fa()] for a related function of the package psych

unirel <- function(x, psych.include = TRUE) {
    ##########################################################################
    # Obtaining reliability estimates
    ###########################################################################
    r <- get_cov(x, cor = TRUE) # correlation
    alpha <- reliacoef::alpha(x)
    std_alpha <- std_alpha(x)
    lambda2 <- mu1(x)
    mu2 <- mu2(x)
    mu3 <- mu3(x)
    mu4 <- mu4(x)
    feldt <- feldt(x)
    gilmer <- gilmer(x)
    joreskog <- joreskog(x)
    hancock <- hancock(x)
    heise <- heise(x)
    kaisercaffrey <- kaisercaffrey(x)
    if (psych.include) {
        stopifnot(requireNamespace("psych"))
        #omega_psych <- psych::omega(r)$omega.tot
        glb.algebraic <- psych::glb.algebraic(x)$glb[1]
        glb.fa <- psych::glb.fa(r)$glb[1]
    } else {
        glb.algebraic <- glb.fa <- NULL
    }
    ##########################################################################
    # Printingpsych estimates
    ##########################################################################
    #cat("omega total (nfactors = 1) obtained from psych package            ", omega_psych, "\n
    cat("GLB.algebraic (greatest lower bound) obtained from psych package  ", glb.algebraic, "\n")
    cat("GLB.fa (greatest lower bound) obtained from psych package         ", glb.fa, "\n")
    ##########################################################################
    # Invisible return
    ##########################################################################
    unirel <- list(alpha = alpha, std_alpha = std_alpha, lambda2 = lambda2,
                   mu2 = mu2, mu3 = mu3, mu4 = mu4, feldt = feldt,
                   gilmer = gilmer, joreskog = joreskog, hancock = hancock,
                   heise = heise, kaisercaffrey = kaisercaffrey,
                   #omega_psych = omega_psych,
                   glb.algebraic = glb.algebraic, glb.fa = glb.fa)
    invisible(unirel)
}
