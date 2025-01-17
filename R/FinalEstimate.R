## * FinalEstimate (documentation)
#' @title calculate de-biased estimate
#'
#' @param Info.d Information at all decision analyses up to stage where study was stopped (should include information at final analysis if stopped at final analysis)
#' @param Info.i Information at all interim analyses up to stage where study was stopped
#' @param ck,ck.unrestricted decision boundaries for all decision analyses up to stage where study was stopped (should include final boundary if stopped at final analysis).
#' ck is possibly with restriction (when cNotBelowFixedc=TRUE) and ck.unrestricted always without.
#' @param lk lower bounds up to stage where study was stopped
#' @param uk upper bounds up to stage where study was stopped
#' @param reason.interim motivation for stopping or continuing at interim. Use to handle special cases (skipped interim because reach Imax, ...)
#' @param kMax maximum number of analyses
#' @param estimate naive estimate (e.g. using  ML or REML).
#' @param lower,upper lower and upper bounds for the MUE used to perform the root finding.
#' @param statistic naive test statistic (e.g. using  ML or REML).
#' @param method  method 1, 2 or 3
#' @param bindingFutility [logical]  whether the futility stopping rule is binding.
#' @param cNotBelowFixedc [logical] whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
#' @param continuity.correction [logical] whether to add the probability of stopping between ck and ck.uncorrected to ensure continuity of the p-value across stages.
#' When used the p-value will always be greater than this probability of stopping bettwen ck and ck.uncorrected.
#' @param tolerance [numeric] acceptable discrepancy to the objective level when evaluating the confidence intervals and median unbiased estimate.

## * FinalEstimate (code)
#' @export
FinalEstimate <- function(Info.d,  
                          Info.i,  
                          ck,
                          ck.unrestricted,
                          lk,  
                          uk,
                          reason.interim,
                          kMax, 
                          estimate, lower, upper,
                          statistic,
                          method,
                          bindingFutility,
                          cNotBelowFixedc,
                          continuity.correction,
                          tolerance,
                          FCT.p_value){ 

    f <- function(delta){
        do.call(FCT.p_value, list(Info.d=Info.d,
                                  Info.i=Info.i,
                                  ck=ck,
                                  ck.unrestricted=ck.unrestricted,
                                  lk=lk,
                                  uk=uk,
                                  reason.interim=reason.interim,
                                  kMax=kMax,
                                  statistic=statistic,
                                  delta=delta,
                                  method=method,
                                  bindingFutility=bindingFutility,
                                  cNotBelowFixedc=cNotBelowFixedc,
                                  continuity.correction=continuity.correction)
                )
    }

    se <- estimate/statistic ## equal to sqrt(1/Info.d[length(Info.d)]) except when decreasing information

    if(!is.na(lower) & !is.na(upper)){
        lowerBound <- lower
        upperBound <- upper
    }else{
        lowerBound <- estimate - 4*se
        upperBound <- estimate + 4*se
    }

    res <- try(stats::uniroot(function(x){f(x) - 0.5},
                              lower = lowerBound[1],
                              upper = upperBound[1],
                              extendInt = "yes",
                              tol = 1e-10),
               silent = TRUE)
    if(inherits(res,"try-error")){
        res <- suppressWarnings(stats::optim(fn = function(x){(f(x) - 0.5)^2},
                                             par = (lowerBound[1] + upperBound[1])/2,
                                             method = "Nelder-Mead"))
        res$iter <- unname(res$counts["function"])
        res$root <- unname(res$par)
        res$f.root <- res$value
    }

    if(abs(res$f.root)>tolerance){
        res$root <- NA
    }

    ## ** export
    out <- res$root
    attr(out,"error") <- res$f.root
    attr(out,"iter") <- res$iter
    return(out)
}
