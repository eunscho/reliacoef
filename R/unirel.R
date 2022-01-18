#' Obtain various unidimensional reliablity coefficients
#'
#' R packages that provide unidimensional reliability coefficients include
#' Lambda4 and psych. These are usually separate functions, so comparing
#' multiple reliability estimates can take a long time. This function shows the
#'  reliability coefficients newly added in this package together with the
#'  reliability coefficients provided by Lambda4 and psych. This function can
#'  provide Table 2 of Cho(in press) immediately.
#'
#'@author Eunseong Cho, \email{bene@kw.ac.kr}
#'@param x a dataframe or a matrix (unidimensional)
#'@param Lambda4.include Whether to include the reliability coefficients
#'provided by the package Lambda4 (lambda5, lambda6, max_lambda, lambda4_max,
#'lambda4_75)
#'@param psych.include Whether to include reliability coefficients
#'(GLB.algebraic, GLB.fa) provided by the package psych
#'@export unirel
#' @references Cho, E. (in press). Neither Cronbach's alpha nor McDonald's
#' omega: A comment on Sijtsma and Pfadt. Psychometrika.
#' @examples unirel(Graham1)
#' @seealso [Lambda4::lambda5()] for a related function of the package Lambda4
#' @seealso [Lambda4::lambda6()] for a related function of the package Lambda4
#' @seealso [Lambda4::quant.lambda4()] for a related function of the package Lambda4
#' @seealso [psych::GLB.algebraic()] for a related function of the package psych
#' @seealso [psych::GLB.fa()] for a related function of the package psych

unirel <- function(x, Lambda4.include = TRUE, psych.include = TRUE) {
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
    if (Lambda4.include) {
        stopifnot(requireNamespace("Lambda4"))
        #omega_Lambda4 <- Lambda4::omega.tot(x)[[1]]
        lambda5 <- as.numeric(Lambda4::lambda5(x, missing = "pairwise"))
        lambda6 <- as.numeric(Lambda4::lambda6(x, missing = "pairwise"))
        lambda4_max <- Lambda4::quant.lambda4(x, missing = "pairwise",
                                              quantiles = 1)$lambda4.quantile
        lambda4_75 <- Lambda4::quant.lambda4(x, missing = "pairwise",
                                             quantiles = 0.75)$lambda4.quantile
    } else {
        lambda5 <- lambda6 <- lambda4_max <- lambda4_75 <- NULL
    }
    if (psych.include) {
        stopifnot(requireNamespace("psych"))
        #omega_psych <- psych::omega(r)$omega.tot
        glb.algebraic <- psych::glb.algebraic(x)$glb[1]
        glb.fa <- psych::glb.fa(r)$glb[1]
    } else {
        glb.algebraic <- glb.fa <- NULL
    }
    ##########################################################################
    # Printing Lambda4 & psych estimates
    ##########################################################################
    #cat("omega total (factors = 1) obtained from Lambda4 package           ", omega_Lambda4, "\n")
    #cat("omega total (nfactors = 1) obtained from psych package            ", omega_psych, "\n")
    cat("Guttman's lambda5                                                 ", lambda5, "\n")
    cat("Guttman's lambda6                                                 ", lambda6, "\n")
    cat("Maximum among all possible split-half reliability(lambda4)        ", lambda4_max, "\n")
    cat("75th percentile among all possible split-half reliability(lambda4)", lambda4_75, "\n")
    cat("GLB.algebraic (greatest lower bound) obtained from psych package  ", glb.algebraic, "\n")
    cat("GLB.fa (greatest lower bound) obtained from psych package         ", glb.fa, "\n")
    ##########################################################################
    # Invisible return
    ##########################################################################
    unirel <- list(alpha = alpha, std_alpha = std_alpha, lambda2 = lambda2,
                   mu2 = mu2, mu3 = mu3, mu4 = mu4, feldt = feldt,
                   gilmer = gilmer, joreskog = joreskog, hancock = hancock,
                   heise = heise, kaisercaffrey = kaisercaffrey,
                   #omega_Lambda4 = omega_Lambda4,
                   lambda5 = lambda5,lambda6 = lambda6,
                   lambda4_max = lambda4_max, lambda4_75 = lambda4_75,
                   #omega_psych = omega_psych,
                   glb.algebraic = glb.algebraic, glb.fa = glb.fa)
    invisible(unirel)
}
