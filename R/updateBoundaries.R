## * updateBoundaries (documentation)
#' @title Update Boundaries of a GSD
#' @description Recompute the boundaries based according to the current information.
#'
#' @param object Object of type \code{delayedGSD}, typically output from \code{\link{CalcBoundaries}}.
#' @param Info.i [numeric vector of size k] Optional argument used to update the information at interim (only past or current information). 
#' @param Info.d [numeric vector of size k] Optional argument used to update the information at decision or final (observed or predicted information).
#' @param k [integer] Index of the analysis.
#' @param type.k [character] Type of analysis: \code{"interim"} (after continuing recruitment),
#' \code{"decision"} (after stopping recruitment for efficacy or futility),
#' or \code{"final"} (after reaching the last stage of the trial).
#' @param trace [logical] should the execution of the function be traced?
#' @param update.stage [logical] should the arguments \code{k} and \code{type.k} be used to update to stage of the trial?
#' 

## * updateBoundaries (example)
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
#'
#' #### Analyse data at the first interim ####
#' theInterimData <- SelectData(theData$d, t = tau.i)
#' 
#' myLMM <- analyzeData(theInterimData)
#' myBound1 <- update(myBound0, delta = myLMM)
#' print(myBound1)
#' print(myBound1, planned = FALSE)
#' print(myBound1, planned = "only")
#'
#' par(mfrow = c(1,2))
#' plot(myBound1, planned = "only")
#' plot(myBound1)
#' 
#' #### Analyse data at the final ####
#' theFinalData <- SelectData(theData$d, t = 1e7, Delta.t = theDelay) 
#' 
#' myLMM <- analyzeData(theFinalData)
#' myBound2 <- update(myBound1, delta = myLMM)
#' plot(myBound2)

## * updateBoundaries (code)
#' @export
updateBoundaries <- function(object, delta, Info.i, Info.d, k, type.k, update.stage, trace = FALSE){

    kMax <- object$kMax
    Info.max <- object$planned$Info.max
    uk <- object$uk
    lk <- object$lk
    ck <- object$ck
    bindingFutility <- object$bindingFutility
    method <- object$method

    ## ** update information
    if(missing(Info.i)){
        Info.i <- object$Info.i
    }else{
        if(type.k=="decision"){
            object$Info.i[k+1] <- Info.i
        }else{
            object$Info.i[k] <- Info.i
        }
    }
    if(missing(Info.d)){
        Info.d <- object$Info.d
    }else{
        if(type.k=="decision"){
            object$Info.d[k+1] <- Info.d
        }else{
            object$Info.d[k] <- Info.d
        }
    }
    if(!missing(delta)){
        if(type.k=="decision"){
            object$delta$estimate[k+1] <- delta
        }else{
            object$delta$estimate[k] <- delta
        }
    }
    
    ## ** update boundaries
   #browser()
    ## *** at interim
    if(type.k == "interim"){

        if(k>1 && (Info.i[k] < Info.i[k-1])){
            
            ## continue to the next interim when information decreased
            object$lk[k]  <- -Inf
            object$uk[k]  <- +Inf
            object$ck[k]  <- NA
            object$conclusion["interim",k] <- "continue"
            object$conclusion["reason.interim",k] <- "decreasing information"
            object$alphaSpent[k] <- object$alphaSpent[k-1]
            object$betaSpent[k] <- object$betaSpent[k-1]
            
        }else if((Info.i[k] >= Info.max) || (!is.na(Info.d[k]) && Info.d[k] >= Info.max)){
            ## stop and do decision when Imax is already reached or anticipated to be reached
            object$lk[k]  <- -Inf
            object$uk[k]  <- +Inf
            object$ck[k]  <- NA
            object$conclusion["interim",k] <- "stop"
            object$conclusion["reason.interim",k] <- "Imax reached"
            object$alphaSpent[k] <- object$alpha
            object$betaSpent[k] <- object$beta
        }else{
            ## usual evaluation
            if(method==1){
                newBounds <- updateMethod1(rho_alpha = object$planned$rho_alpha,
                                           rho_beta = object$planned$rho_beta,
                                           alpha = object$alpha, alphaSpent = object$alphaSpent,
                                           beta = object$beta, betaSpent = object$betaSpent, 
                                           Kmax = object$kMax,
                                           Info.max = object$planned$Info.max,
                                           uk = object$uk,
                                           lk = object$lk,
                                           k = k, type.k = type.k, ImaxAnticipated = FALSE,
                                           InfoR.i = object$Info.i/object$planned$Info.max,
                                           InfoR.d = object$Info.d/object$planned$Info.max,
                                           delta = object$planned$delta, ## object$delta$estimate, 
                                           alternative = object$alternative,
                                           binding = bindingFutility,
                                           Trace = trace,
                                           cMin = object$cMin,
                                           PowerCorrection=object$PowerCorrection)
            } else if(method==2){
                newBounds <- updateMethod2(rho_alpha = object$planned$rho_alpha,
                                           rho_beta = object$planned$rho_beta,
                                           alpha = object$alpha, alphaSpent = object$alphaSpent,
                                           beta = object$beta, betaSpent = object$betaSpent, 
                                           Kmax = object$kMax,
                                           Info.max = object$planned$Info.max,
                                           uk = object$uk,
                                           lk = object$lk,
                                           k = k, type.k = type.k, ImaxAnticipated = FALSE,
                                           InfoR.i = object$Info.i/object$planned$Info.max,
                                           InfoR.d = object$Info.d/object$planned$Info.max,
                                           delta = object$planned$delta, ## object$delta$estimate, 
                                           alternative = object$alternative,
                                           binding=bindingFutility,
                                           Trace = trace,
                                           cMin = object$cMin)
            } else if(method==3){
                newBounds <- updateMethod3(rho_alpha = object$planned$rho_alpha,
                                           rho_beta = object$planned$rho_beta,
                                           alpha = object$alpha, alphaSpent = object$alphaSpent,
                                           beta = object$beta, betaSpent = object$betaSpent, 
                                           Kmax = object$kMax,
                                           Info.max = object$planned$Info.max,
                                           uk = object$uk,
                                           lk = object$lk,
                                           k = k, type.k = type.k, ImaxAnticipated = FALSE,
                                           InfoR.i = object$Info.i/object$planned$Info.max,
                                           InfoR.d = object$Info.d/object$planned$Info.max,
                                           delta = object$planned$delta, ## object$delta$estimate, 
                                           alternative = object$alternative,
                                           binding=bindingFutility,
                                           Trace = trace)
            }
            
            object$lk[k]  <- newBounds$lk[k]
            object$uk[k]  <- newBounds$uk[k]
            object$ck[k]  <- newBounds$ck
            object$ck.unrestricted[k]  <- newBounds$ck.unrestricted
            object$alphaSpent[k] <- newBounds$alphaSpent[k]
            object$betaSpent[k] <- newBounds$betaSpent[k]
            
            object$lk[(k+1):kMax] <- NA
            object$uk[(k+1):kMax]  <- NA
            
            
            #if(k<kMax){
            #    object$lk[(k+1):kMax]  <- NA
            #    
            #    if(k<kMax-1){
            #        object$ck[(k+1):kMax]  <- NA 
            #        object$ck.unrestricted[(k+1):kMax]  <- NA 
            #    }
            #}
        } 

    }

    ## *** at decision
    if(type.k == "decision"){

        if(Info.d[k] < Info.i[k]){ ## if information has decreased since interim analysis
          
            object$conclusion["comment.decision",k] <- "decreasing information"
            warning("Information has decreased between interim and decision, replacing information at decision with information at interim + epsilon. \n")

            Info.d[k] <- Info.i[k]+Info.i[k]/10000
            
            ## Possible solution: when estimating ck by balancing the probability of reversal, add a term corresponding to the type 1 error we should have spent vs. the type 1 error we spent.
            ## usual evaluation
            if(method==1){
                newBounds <- updateMethod1(rho_alpha = object$planned$rho_alpha,
                                           rho_beta = object$planned$rho_beta,
                                           alpha = object$alpha, alphaSpent = object$alphaSpent,
                                           beta = object$beta, betaSpent = object$betaSpent, 
                                           Kmax = object$kMax,
                                           Info.max = object$planned$Info.max,
                                           uk = object$uk,
                                           lk = object$lk,
                                           k = k, type.k = type.k, ImaxAnticipated = (object$conclusion["reason.interim",k]=="Imax reached"),
                                           InfoR.i = object$Info.i/object$planned$Info.max,
                                           InfoR.d = Info.d/object$planned$Info.max,
                                           delta = object$planned$delta, ## object$delta$estimate, 
                                           alternative = object$alternative,
                                           binding=bindingFutility,
                                           Trace = trace,
                                           cMin = object$cMin,
                                           PowerCorrection=object$PowerCorrection)
                
            } else if(method==2){
                newBounds <- updateMethod2(rho_alpha = object$planned$rho_alpha,
                                           rho_beta = object$planned$rho_beta,
                                           alpha = object$alpha, alphaSpent = object$alphaSpent,
                                           beta = object$beta, betaSpent = object$betaSpent, 
                                           Kmax = object$kMax,
                                           Info.max = object$planned$Info.max,
                                           uk = object$uk,
                                           lk = object$lk,
                                           k = k, type.k = type.k, ImaxAnticipated = (object$conclusion["reason.interim",k]=="Imax reached"),
                                           InfoR.i = object$Info.i/object$planned$Info.max,
                                           InfoR.d = Info.d/object$planned$Info.max,
                                           delta = object$planned$delta, ## object$delta$estimate, 
                                           alternative = object$alternative,
                                           binding=bindingFutility,
                                           Trace = trace,
                                           cMin = object$cMin)
                
            } else if(method==3){
                newBounds <- updateMethod3(rho_alpha = object$planned$rho_alpha,
                                           rho_beta = object$planned$rho_beta,
                                           alpha = object$alpha, alphaSpent = object$alphaSpent,
                                           beta = object$beta, betaSpent = object$betaSpent, 
                                           Kmax = object$kMax,
                                           Info.max = object$planned$Info.max,
                                           uk = object$uk,
                                           lk = object$lk,
                                           k = k, type.k = type.k, ImaxAnticipated = (object$conclusion["reason.interim",k]=="Imax reached"),
                                           InfoR.i = object$Info.i/object$planned$Info.max,
                                           InfoR.d = Info.d/object$planned$Info.max,
                                           delta = object$planned$delta, ## object$delta$estimate, 
                                           alternative = object$alternative,
                                           binding=bindingFutility,
                                           Trace = trace)                
            }
            

            object$ck[k]  <- newBounds$ck
            if(k<kMax){
                object$ck.unrestricted[k]  <- newBounds$ck.unrestricted
            }
            if(k<kMax-1){
                object$ck[(k+1):kMax]  <- NA                    
            }
            
        }else{

            ## NOTE: include special case Info.d[k] >= Info.max which is treated as usual - conservative procedure

            ## usual evaluation
            if(method==1){
                newBounds <- updateMethod1(rho_alpha = object$planned$rho_alpha,
                                           rho_beta = object$planned$rho_beta,
                                           alpha = object$alpha, alphaSpent = object$alphaSpent,
                                           beta = object$beta, betaSpent = object$betaSpent, 
                                           Kmax = object$kMax,
                                           Info.max = object$planned$Info.max,
                                           uk = object$uk,
                                           lk = object$lk,
                                           k = k, type.k = type.k, ImaxAnticipated = (object$conclusion["reason.interim",k]=="Imax reached"),
                                           InfoR.i = object$Info.i/object$planned$Info.max,
                                           InfoR.d = object$Info.d/object$planned$Info.max,
                                           delta = object$planned$delta, ## object$delta$estimate, 
                                           alternative = object$alternative,
                                           binding=bindingFutility,
                                           Trace = trace,
                                           cMin = object$cMin,
                                           PowerCorrection=object$PowerCorrection)
                
            } else if(method==2){
                newBounds <- updateMethod2(rho_alpha = object$planned$rho_alpha,
                                           rho_beta = object$planned$rho_beta,
                                           alpha = object$alpha, alphaSpent = object$alphaSpent,
                                           beta = object$beta, betaSpent = object$betaSpent, 
                                           Kmax = object$kMax,
                                           Info.max = object$planned$Info.max,
                                           uk = object$uk,
                                           lk = object$lk,
                                           k = k, type.k = type.k, ImaxAnticipated = (object$conclusion["reason.interim",k]=="Imax reached"),
                                           InfoR.i = object$Info.i/object$planned$Info.max,
                                           InfoR.d = object$Info.d/object$planned$Info.max,
                                           delta = object$planned$delta, ## object$delta$estimate, 
                                           alternative = object$alternative,
                                           binding=bindingFutility,
                                           Trace = trace,
                                           cMin = object$cMin)
                
            } else if(method==3){
                newBounds <- updateMethod3(rho_alpha = object$planned$rho_alpha,
                                           rho_beta = object$planned$rho_beta,
                                           alpha = object$alpha, alphaSpent = object$alphaSpent,
                                           beta = object$beta, betaSpent = object$betaSpent, 
                                           Kmax = object$kMax,
                                           Info.max = object$planned$Info.max,
                                           uk = object$uk,
                                           lk = object$lk,
                                           k = k, type.k = type.k, ImaxAnticipated = (object$conclusion["reason.interim",k]=="Imax reached"),
                                           InfoR.i = object$Info.i/object$planned$Info.max,
                                           InfoR.d = object$Info.d/object$planned$Info.max,
                                           delta = object$planned$delta, ## object$delta$estimate, 
                                           alternative = object$alternative,
                                           binding = bindingFutility,
                                           Trace = trace)
            }
            object$ck[k]  <- newBounds$ck
            object$ck.unrestricted[k]  <- newBounds$ck.unrestricted
            
            if(k<kMax-1){
                object$ck[(k+1):kMax]  <- NA                    
                object$ck.unrestricted[(k+1):kMax]  <- NA
            }
        }
    }

    
    ## *** at final
    if(type.k == "final"){ ##  that could be replace by a call to CalcBoundaries

        if(Info.d[k] < Info.i[k-1]){
            object$conclusion["comment.decision",k] <- "decreasing information"
            warning("boundaries will not be computed correctly")
        }

        if(method==1){
            newBounds <- updateMethod1(rho_alpha = object$planned$rho_alpha,
                                       rho_beta = object$planned$rho_beta,
                                       alpha = object$alpha, alphaSpent = object$alphaSpent,
                                       beta = object$beta, betaSpent = object$betaSpent, 
                                       Kmax = object$kMax,
                                       Info.max = object$planned$Info.max,
                                       uk = object$uk,
                                       lk = object$lk,
                                       k = k, type.k = type.k, ImaxAnticipated = FALSE,
                                       InfoR.i = object$Info.i/object$planned$Info.max,
                                       InfoR.d = object$Info.d/object$planned$Info.max,
                                       delta = object$planned$delta, ## object$delta$estimate, 
                                       alternative = object$alternative,
                                       binding=bindingFutility,
                                       Trace = trace,
                                       cMin = object$cMin,
                                       PowerCorrection=object$PowerCorrection)
        } else if(method==2){
            newBounds <- updateMethod2(rho_alpha = object$planned$rho_alpha,
                                       rho_beta = object$planned$rho_beta,
                                       alpha = object$alpha, alphaSpent = object$alphaSpent,
                                       beta = object$beta, betaSpent = object$betaSpent, 
                                       Kmax = object$kMax,
                                       Info.max = object$planned$Info.max,
                                       uk = object$uk,
                                       lk = object$lk,
                                       k = k, type.k = type.k, ImaxAnticipated = FALSE,
                                       InfoR.i = object$Info.i/object$planned$Info.max,
                                       InfoR.d = object$Info.d/object$planned$Info.max,
                                       delta = object$planned$delta, ## object$delta$estimate, 
                                       alternative = object$alternative,
                                       binding=bindingFutility,
                                       Trace = trace,
                                       cMin = object$cMin)
        } else if(method==3){
            newBounds <- updateMethod3(rho_alpha = object$planned$rho_alpha,
                                       rho_beta = object$planned$rho_beta,
                                       alpha = object$alpha, alphaSpent = object$alphaSpent,
                                       beta = object$beta, betaSpent = object$betaSpent, 
                                       Kmax = object$kMax,
                                       Info.max = object$planned$Info.max,
                                       uk = object$uk,
                                       lk = object$lk,
                                       k = k, type.k = type.k, ImaxAnticipated = FALSE,
                                       InfoR.i = object$Info.i/object$planned$Info.max,
                                       InfoR.d = object$Info.d/object$planned$Info.max,
                                       delta = object$planned$delta, ## object$delta$estimate, 
                                       alternative = object$alternative,
                                       binding=bindingFutility,
                                       Trace = trace)
        }
        #object$uk[k]  <- newBounds$uk[k]
        #object$lk[k]  <- newBounds$lk[k]
        object$ck[k]  <- newBounds$ck
        object$ck.unrestricted[k]  <- newBounds$ck.unrestricted
        object$alphaSpent[k] <- newBounds$alphaSpent[k]
        object$betaSpent[k] <- newBounds$betaSpent[k]
    }

    ## ** export stage
    if(update.stage){
        object$stage$k <- k
        object$stage$type <- type.k
    }

    ## ** export object
    return(object)
}


