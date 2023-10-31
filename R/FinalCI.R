## * FinalCI (documentation)
#' @title calculate confidence intervals at the end of the study
#' 
#' @param Info.d Information at all decision analyses up to stage where study was stopped
#' @param Info.i Information at all interim analyses up to stage where study was stopped
#' @param ck,ck.unrestricted decision boundaries for all decision analyses up to stage where study was stopped (should include final boundary if stopped at final analysis).
#' ck is possibly with restriction (when cNotBelowFixedc=TRUE) and ck.unrestricted always without.
#' @param lk lower bounds up to stage where study was stopped
#' @param uk upper bounds up to stage where study was stopped
#' @param reason.interim motivation for stopping or continuing at interim. Use to handle special cases (skipped interim because reach Imax, ...)
#' @param kMax maximum number of analyses
#' @param conf.level confidence level (to get a 100*(1-alpha)\% CI)
#' @param estimate naive estimate (e.g. using  ML or REML).
#' @param method  method 1, 2 or 3
#' @param bindingFutility [logical]  whether the futility stopping rule is binding.
#' @param cNotBelowFixedc [logical] whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
#' @param continuity.correction [logical] whether to add the probability of stopping between ck and ck.uncorrected to ensure continuity of the p-value across stages.
#' When used the p-value will always be greater than this probability of stopping bettwen ck and ck.uncorrected.
#' @param tolerance [numeric] acceptable discrepancy to the objective level when evaluating the confidence intervals and median unbiased estimate.
#' @param FCT.p_value [function] function used to compute the p-value.

## * FinalCI (code)
#' @export
FinalCI <- function(Info.d,  
                    Info.i,  
                    ck,  
                    ck.unrestricted,  
                    lk,  
                    uk,
                    reason.interim,
                    kMax, 
                    conf.level,
                    estimate,
                    method,
                    bindingFutility,
                    cNotBelowFixedc,
                    continuity.correction,
                    tolerance,
                    conclusion,                    
                    FCT.p_value){

    alpha <- 1 - conf.level

    ## ** objective function
    f <- function(delta){
        do.call(FCT.p_value, list(Info.d=Info.d,
                                  Info.i=Info.i,
                                  ck=ck,
                                  ck.unrestricted=ck.unrestricted,
                                  lk=lk,
                                  uk=uk,
                                  reason.interim=reason.interim,
                                  kMax=kMax,
                                  estimate=estimate,
                                  delta=delta,
                                  method=method,
                                  bindingFutility=bindingFutility,
                                  cNotBelowFixedc=cNotBelowFixedc,
                                  continuity.correction=continuity.correction)
                )
    }

    ## ** initialization
    lowerBound <- c(lbnd = estimate - 4*sqrt(1/Info.d[length(Info.d)]),
                    ubnd = estimate - 0.5*sqrt(1/Info.d[length(Info.d)]))
    upperBound <- c(lbnd = estimate + 0.5*sqrt(1/Info.d[length(Info.d)]),
                    ubnd = estimate + 4*sqrt(1/Info.d[length(Info.d)]))
    if(conclusion=="efficacy"){
        lowerBound <- pmax(0,lowerBound)
        upperBound <- pmax(0,upperBound)        
    }else if(conclusion=="futility"){
        lowerBound <- pmin(0,lowerBound)
        upperBound <- pmin(0,upperBound)        
    }

    f_lowerBound <- c(f(lowerBound[1]),f(lowerBound[2]))
    f_upperBound <- c(f(upperBound[1]),f(upperBound[2]))

    ## ** lower bound of the CI
    ## if(sign(f_lowerBound[1] - alpha/2) != sign(f_upperBound[1] - alpha/2)) ## not good enougth: f(x) is a bit stochastic and can lead to error due to same sign on bounds
    lbnd <- try(stats::uniroot(function(x){(f(x) - alpha/2)},
                               lower = lowerBound[1],
                               upper = upperBound[1],
                               tol = 1e-10), silent = TRUE)
    if(inherits(lbnd,"try-error")){
        lbnd <- suppressWarnings(stats::optim(fn = function(x){(f(x) - alpha/2)^2},
                                              par = (lowerBound[1] + upperBound[1])/2,
                                              method = "Nelder-Mead"))
        ## Advarselsbesked:
        ## I stats::optim(fn = function(x) { :
        ##   one-dimensional optimization by Nelder-Mead is unreliable:
        ## use "Brent" or optimize() directly
        ## optimize requires lower,upper values which we are unable to provide here
        lbnd$iter <- unname(lbnd$counts["function"])
        lbnd$root <- lbnd$par
        lbnd$f.root <- lbnd$value
    }

    if(abs(lbnd$f.root)>tolerance){
        lbnd$root <- NA
    }

    ## ** lower bound of the CI
    ## if(sign(1 - f_lowerBound[2] - alpha/2)!=sign(1 - f_upperBound[2] - alpha/2)) ## not good enougth: f(x) is a bit stochastic and can lead to error due to same sign on bounds
    ubnd <- try(stats::uniroot(function(x){(1 - f(x) - alpha/2)},
                           lower = lowerBound[2],
                           upper = upperBound[2],
                           tol = 1e-10), silent = TRUE)

    if(inherits(ubnd,"try-error")){
        ubnd <- suppressWarnings(stats::optim(fn = function(x){(1 - f(x) - alpha/2)^2},
                                              par = (lowerBound[2] + upperBound[2])/2,
                                              method = "Nelder-Mead"))
        ## Advarselsbesked:
        ## I stats::optim(fn = function(x) { :
        ##   one-dimensional optimization by Nelder-Mead is unreliable:
        ## use "Brent" or optimize() directly
        ## optimize requires lower,upper values which we are unable to provide here
        ubnd$iter <- unname(ubnd$counts["function"])
        ubnd$root <- ubnd$par
        ubnd$f.root <- ubnd$value
    }

    if(abs(ubnd$f.root)>tolerance){
        ubnd$root <- NA
    }

    ## ** export
    out <- c(lower = lbnd$root, upper = ubnd$root)
    attr(out,"error") <- c(lower = lbnd$f.root, upper = ubnd$f.root)
    attr(out,"iter") <- c(lower = lbnd$iter, upper = ubnd$iter)
    return(out)
}
