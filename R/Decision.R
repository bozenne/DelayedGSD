## * Decision (documentation)
#' @title Take Decision in a GSD
#' @description Maps the statistical results at an interim, decision or final analysis into a decision regarding whether to stop recruitment (interim) or whether to reject the null hypothesis (decision/final). Stopping boundaries are updated based on observed information and correct p-values, confidence intervals and point estimates are given.
#' 
#' @param object object of class \code{delayedGSD}, typically output from \code{\link{CalcBoundaries}}.
#' @param PositiveIsGood [logical] whether a positive effect is considered beneficial (TRUE/FALSE)
#' @param k [integer] index of the analysis.
#' @param type.k [character] type of analysis: \code{"interim"} (after continuing recruitment),
#' \code{"decision"} (after stopping recruitment for efficacy or futility),
#' or \code{"final"} (after reaching the last stage of the trial).
#' @param trace [logical] should the execution of the function be traced?
#' 
#' @details Function that maps the statistical results at interim (or final) analysis into decision to reject, or continue or stop inclusind subjects.
#' @return ff
#' @author Paul Blanche
#' 
## * Decision (example)
#' @examples
#'
#' #### Planning #####
#' theAlpha <- 0.025
#' theBeta <- 0.2
#' theDelta <- 1.5
#' theK <- 2
#' theN <- 82
#' 
#' myBound0 <- CalcBoundaries(kMax=theK,
#'                      sided=1,
#'                      alpha=theAlpha,
#'                      beta=theBeta,
#'                      InfoR.i=c(0.5,1),
#'                      rho_alpha=2,
#'                      rho_beta=2,
#'                      method=1,
#'                      delta=theDelta,
#'                      InfoR.d=0.55)
#'
#' #### Simulate data ####
#' ## generate data with all data for in case trial completes
#' set.seed(10)
#' theData <- GenData(n=theN*2,delta=theDelta*0.8,ar=5)  
#' 
#' theAR <- 10  #accrual rate (pt per month)
#' theDelay <- 0.7500001  #time in months to process data
#' tau.i <- theData$d$t3[theN + ceiling(theAR*theDelay)] #time point at which to do IA
#'
#' theObsData <- SelectData(theData$d, t = tau.i, Delta.t = theDelay)
#'
#' #### Analyse data at interim ####
#' myLMM <- analyzeData(theObsData)
#' myBound1 <- updateBoundaries(myBound0, lmm = myLMM, k = 1, type.k = "interim")
#' myInterim1 <- Decision(myBound1, k = 1, type.k = "interim")
#' myInterim1$conclusion

## * Decision (documentation)
#' @export
Decision <- function(object,
                     k,
                     type.k,
                     PositiveIsGood = TRUE, # whether positive effect is good (i.e. positive trial)
                     trace = TRUE){
  
    ## ** identify stage of the trial
    kMax <- object$kMax
    
    if(trace>0){
        if(type.k == "final"){
            cat("Decision to be taken at the final analysis  (stage ",k,") \n",sep="")
        }else if(type.k == "decision"){
            cat("Decision to be taken at the decision analysis of stage ",k," \n",sep="")
        }else if(type.k == "interim"){
            cat("Decision to be taken at the interim analysis of stage ",k," \n",sep="")
        }
    }

    
    ## ** extract information, test statistic, and boundaries
    Info.max <- object$planned$Info.max
    method <- object$method
#browser()
    outInfo <- coef(object, type = "information", planned = FALSE)
    outBound <- coef(object, type = "boundary", planned = FALSE)
    Info.i <- outInfo[,"Interim"]
    Info.d <- outInfo[,"Decision"]
    uk <- outBound[,"Ebound"]
    lk <- outBound[,"Fbound"]
    ck <- outBound[,"Cbound"]

    Z <- object$delta[k+(type.k=="decision"),"statistic"]

    if(object$alternative=="less"){
        Z <- -Z
        if(trace){message("Negative effect is good: change sign of Z")}
    }

    ## ** decision at interim
    if(type.k=="interim"){

        ## Special cases at decision dealt with in updateBoundary
        ## if(k>1 && (Info.i[k] < Info.i[k-1])){ 
        ## } else if((Info.i[k] >= Info.max) || (Info.d[k] >= Info.max)){
        if(Z>uk[k]){
            object$conclusion["interim",k] <- "stop"
            object$conclusion["reason.interim",k] <- "efficacy"
        } else if (Z<lk[k]){
            object$conclusion["interim",k] <- "stop"
            object$conclusion["reason.interim",k] <- "futility"
        } else {
            object$conclusion["interim",k] <- "continue"
            object$conclusion["reason.interim",k] <- "no boundary crossed"
        }
    }
    
    ## ** decision at decision analysis
    if(type.k=="decision"){
                                        
        if(Z > ck[k]){
            if(method == 3 && object$conclusion["reason.interim",k] == "futility"){
                ## with method 3 we cannot conclude efficacy if we stopped for futility
                object$conclusion["decision",k] <- "futility"
                object$conclusion["comment.decision",k] <- "stop for futility at interim"
            }else{
                object$conclusion["decision",k] <- "efficacy"
            }
        } else {
            object$conclusion["decision",k] <- "futility"
        }

    }

    ## ** decision at the final analysis  
    if(type.k == "final"){

        if(Z > ck[k]){
            object$conclusion["decision",k] <- "efficacy"
        } else {
            object$conclusion["decision",k] <- "futility"
        }

    }

    ## ** export  
    return(object)
}

