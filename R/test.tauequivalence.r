#' Test the essential tau-equivalence of the data
#' 
#' The goodness-of-fit indices are compared with the parameter estimates of the 
#' essential tau-equivalence model and the cognate model. It can also be used for 
#' the purpose of investigating the unidimensionality of the data.
#' @export test.tauequivalence
#' @author Eunseong Cho, \email{bene@kw.ac.kr}
#' @param data a dataframe or a matrix (unidimensional)
#' @return taueq_cfi
#' @return taueq_tli
#' @return taueq_rmsea
#' @return taueq_df
#' @return taueq_pvalue
#' @return taueq_chisq
#' @return conge_cfi
#' @return conge_tli
#' @return conge_rmsea
#' @return conge_df
#' @return conge_pvalue
#' @return conge_chisq
#' @return diff_df
#' @return diff_chisq
#' @return diff_pvalue
#' @examples test.tauequivalence(Graham1)
#' @references Graham, J. M. (2006). Congeneric and (essentially) tau-equivalent 
#' estimates of score reliability what they are and how to use them. Educational 
#' and Psychological Measurement, 66(6), 930-944.
#' @references Cho, E. (2016). Making reliability reliable: A systematic 
#' approach to reliability coefficients. Organizational Research Methods, 19(4), 
#' 651-682. 
#' @references Cho, E., & Kim, S. (2015). Cronbach’s coefficient alpha: Well 
#' known but poorly understood. Organizational Research Methods, 18(2), 207-230. 
#'
test.tauequivalence <- function(data) {
        cov <- get_cov(data)
        taueq_fit  <- uni_cfa(cov, what = "fit", taueq = TRUE) 
        conge_fit <- uni_cfa(cov, what = "fit", taueq = FALSE) 
    ############################################################
    #   Converting lavaan's result to variables
    ############################################################
        taueq_cfi = taueq_fit[9]
        taueq_tli = taueq_fit[10]
        taueq_rmsea = taueq_fit[23]
        taueq_df = taueq_fit[4]
        taueq_pvalue = taueq_fit[5]
        taueq_chisq = taueq_fit[3]
        conge_cfi = conge_fit[9]
        conge_tli = conge_fit[10]
        conge_rmsea = conge_fit[23]
        conge_df = conge_fit[4]
        conge_pvalue = conge_fit[5]
        conge_chisq = conge_fit[3]
        diff_df = taueq_df - conge_df
        diff_chisq = taueq_chisq - conge_chisq
        diff_pvalue = 1- stats::pchisq(diff_chisq, diff_df)
    ############################################################
    #   print
    ############################################################
        cat("Parameter estimates of the tau-equivalent model\n")
        print(uni_cfa(cov, what = "est", taueq = TRUE))
        cat("Parameter estimates of the congeneric model\n")
        print(uni_cfa(cov, what = "est", taueq = FALSE))
        cat("                     CFI  TLI  RMSEA df chisquare pvalue\n")
        cat(paste("tau-equivalent (A) ", round(taueq_cfi, 3), 
                                         round(taueq_tli, 3), 
                                         round(taueq_rmsea, 3),
                                         round(taueq_df,3), "  ",
                                         round(taueq_chisq, 3),
                                         round(taueq_pvalue, 3), "\n"))
        cat(paste("congeneric     (B) ", round(conge_cfi, 3), 
                                         round(conge_tli, 3), 
                                         round(conge_rmsea, 3),
                                         round(conge_df,3), "  ",
                                         round(conge_chisq, 3),
                                         round(conge_pvalue, 3), "\n"))
        cat(paste("Difference (A - B)                   ", 
                                         round(diff_df, 3), "  ",
                                         round(diff_chisq, 3),
                                         round(diff_pvalue, 3)))
    ############################################################
    #   invisible return
    ############################################################
        test_taueq <- list(taueq_cfi = taueq_cfi,
                           taueq_tli = taueq_tli,
                           taueq_rmsea = taueq_rmsea,
                           taueq_df = taueq_df,
                           taueq_pvalue = taueq_pvalue,
                           taueq_chisq = taueq_chisq,
                           conge_cfi = conge_cfi,
                           conge_tli = conge_tli,
                           conge_rmsea = conge_rmsea,
                           conge_df = conge_df,
                           conge_pvalue = conge_pvalue,
                           conge_chisq = conge_chisq,
                           diff_df = diff_df,
                           diff_chisq = diff_chisq,
                           diff_pvalue = diff_pvalue
        )
        invisible(test_taueq)
}
