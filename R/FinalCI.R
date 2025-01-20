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
#' @param statistic naive test statistic (e.g. using  ML or REML).
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
                    statistic,
                    method,
                    bindingFutility,
                    cNotBelowFixedc,
                    continuity.correction,
                    tolerance,
                    conclusion,                    
                    FCT.p_value){

    alpha <- 1 - conf.level
    se <- estimate/statistic ## equal to sqrt(1/Info.d[length(Info.d)]) except when decreasing information

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
                                  statistic=statistic,
                                  delta=delta,
                                  method=method,
                                  bindingFutility=bindingFutility,
                                  cNotBelowFixedc=cNotBelowFixedc,
                                  continuity.correction=continuity.correction)
                )
    }

    ## ** initialization
    lowerBound <- c(lbnd = estimate - 4*se,
                    ubnd = estimate - 0.5*se)
    upperBound <- c(lbnd = estimate + 0.5*se,
                    ubnd = estimate + 4*se)

    ## ** lower bound of the CI
    lbnd <- try(stats::uniroot(function(x){(f(x) - alpha/2)},
                               lower = lowerBound[1],
                               upper = upperBound[1],
                               extendInt = "upX",
                               tol = 1e-10), silent = TRUE)
    
    if(!inherits(lbnd,"try-error") && !is.null(attr(lbnd$f.root,"error")) && max(attr(lbnd$f.root,"error"))>(tolerance/10)){
        ## f is stochastic due to numerical approximations so by chance is may be above tolerance
        ## replicating the evaluation should provide a more stable value
        lbnd$f.root <- mean(abs(replicate(f(lbnd$root) - alpha/2, n = 10)))
    }

    if(inherits(lbnd,"try-error") || abs(lbnd$f.root)>tolerance){
        
        if(inherits(lbnd,"try-error")){
            starter <- (lowerBound[2] + upperBound[2])/2
        }else{
            starter <- lbnd$root
        }
        ## artificially increase the objective function otherwise optim just set the solution to 0 or 1 as alpha/2 is quite small
        lbnd <- suppressWarnings(stats::optim(fn = function(x){1e3*(f(x) - alpha/2)^2},
                                              par = starter,
                                              upper = upperBound[1],                               
                                              method = "L-BFGS-B"))
        lbnd$iter <- unname(lbnd$counts["function"])
        lbnd$root <- unname(lbnd$par)
        lbnd$f.root <- mean(abs(replicate(f(lbnd$root) - alpha/2, n = 10)))
    }

    if(inherits(lbnd$f.root,"try-error") || abs(lbnd$f.root)>tolerance){
        lbnd$root <- NA
    }

    ## ** upper bound of the CI
    ubnd <- try(stats::uniroot(function(x){(1 - f(x) - alpha/2)},
                               lower = lowerBound[2],
                               upper = upperBound[2],
                               extendInt = "downX",
                               tol = 1e-10), silent = TRUE)

    if(!inherits(ubnd,"try-error") && !is.null(attr(ubnd$f.root,"error")) && max(attr(ubnd$f.root,"error"))>(tolerance/10)){
        ## f is stochastic due to numerical approximations so by chance is may be above tolerance
        ## replicating the evaluation should provide a more stable value
        ubnd$f.root <- mean(abs(replicate(1 - f(ubnd$root) - alpha/2, n=10)))
    }

    if(inherits(ubnd,"try-error") || abs(ubnd$f.root)>tolerance){ 
        if(inherits(ubnd,"try-error")){
            starter <- (lowerBound[2] + upperBound[2])/2
        }else{
            starter <- ubnd$root
        }
        ## artificially increase the objective function otherwise optim just set the solution to 0 or 1 as alpha/2 is quite small
        ubnd <- suppressWarnings(stats::optim(fn = function(x){1e3*(1 - f(x) - alpha/2)^2},
                                              par = starter,
                                              lower = lowerBound[2],
                                              method = "L-BFGS-B"))
        ## 1e3*(1 - f(ubnd$par) - alpha/2)^2
        ubnd$iter <- unname(ubnd$counts["function"])
        ubnd$root <- unname(ubnd$par)
        ubnd$f.root <- mean(abs(replicate(1 - f(ubnd$root) - alpha/2, n=10)))
    }


    if(inherits(ubnd$f.root,"try-error") || abs(ubnd$f.root)>tolerance){
        ubnd$root <- NA
    }

    ## ** export
    out <- c(lower = lbnd$root, upper = ubnd$root)
    attr(out,"error") <- c(lower = lbnd$f.root, upper = ubnd$f.root)
    attr(out,"iter") <- c(lower = lbnd$iter, upper = ubnd$iter)
    return(out)
}
