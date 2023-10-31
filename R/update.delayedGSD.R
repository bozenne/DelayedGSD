## * update.delayedGSD (documentation)
#' @title Perform Next Analysis Stage in a GSD
#' @description Use the newly observed data to perform the next stage of the analysis:
#' fit a mixed model, update the information, recompute the boundaries, take a decision, and possibly correct the estimated treatment effect.
#' 
#' @param object object of class \code{delayedGSD}, typically output from \code{\link{CalcBoundaries}}.
#' @param delta [data.frame or lmmGSD] estimated effect (\code{"estimate"}), standard error (\code{"se"}), test statistic (\code{"statistic"}), degrees of freedom (\code{"df"}), and p-value (\code{"p.value"}) all unadjusted for GSD. Alternatively output of \code{analyzeData}.
#' @param Info.i [numeric] information at the current stage when interim. Not used when argument \code{beta} is a lmmGSD object.
#' @param Info.d [numeric] information at the current stage when decision or final or at the coming decision when interim leading to early stop. Not used when argument \code{beta} is a lmmGSD object.
#' @param k [integer] index of the analysis.
#' @param type.k [character] type of analysis: \code{"interim"} (after continuing recruitment),
#' \code{"decision"} (after stopping recruitment for efficacy or futility),
#' or \code{"final"} (after reaching the last stage of the trial).
#' @param p.value [logical] should the p-value be computed at decision?
#' @param ci [logical] should the confidence intervalsbe computed at decision?
#' @param estimate [logical] should a de-biased estimate be computed at decision? WARNING: this is experiment and not reliable.
#' @param trace [logical] should the execution of the function be traced?
#' @param ... not used, for compatibility with the generic method.

## * update.delayedGSD (examples)
#' @examples
#' 
#' #### Planning #####
#' theAlpha <- 0.025
#' theBeta <- 0.2
#' theDelta <- 1.5
#' theK <- 2
#' theN <- 82
#'
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
#' myLmmI <- analyzeData(theInterimData)
#' myInterim1 <- update(myBound0, delta = myLmmI) ## k = 1, analysis = "interim"
#' print(myInterim1)
#' print(myInterim1, planned = FALSE)
#' print(myInterim1, planned = "only")
#'
#' par(mfrow = c(1,2))
#' plot(myInterim1, planned = "only")
#' plot(myInterim1)
#' 
#' #### Analyse data at the final stage ####
#' theFinalData <- SelectData(theData$d, t = 1e7) 
#' myLmmF <- analyzeData(theFinalData)
#' myFinal <- update(myInterim1, delta = myLmmF) ## k = 2, analysis = "final"
#' myFinal
#' print(myFinal, abreviated = FALSE)
#' plot(myFinal)

## * update.delayedGSD (code)
#' @export
update.delayedGSD <- function(object, delta, Info.i, Info.d, 
                              k = NULL, type.k = NULL, overrule.futility = FALSE,
                              p.value = TRUE, ci = TRUE, estimate = TRUE, trace = TRUE, ...){

    if(overrule.futility){        
        if(any(names(match.call()[-1]) %in% c("object","overrule.futility") == FALSE)){
            args.ignore <- names(match.call()[-1])[names(match.call()[-1]) %in% c("object","overrule.futility") == FALSE]
            warning("Arguments \"",paste(args.ignore, collapse="\", \""),"\" are ignored when overruling the futility bound. \n")
        }
            
        return(.overrule(object))
    }

    ## ** normalize user input
    kMax <- object$kMax
    object.type.k <- object$stage$type
    object.k <- object$stage$k

    if(object.type.k %in% "decision"){
        stop("No more boundary to update when a decision analysis has been performed. \n")
    }
    if(object.type.k %in% "final"){
        stop("No more boundary to update when the final analysis has been performed. \n")
    }
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    
    resStage <- .getStage(object.stage = object$stage,
                          object.conclusion = object$conclusion,
                          kMax = kMax,
                          k = k,
                          type.k = type.k,
                          nextStage = TRUE)

    k <- resStage$k
    type.k <- resStage$type.k
    type.update <- resStage$type.update

    if(trace>0){
        if(object.type.k=="planning"){
            cat("Group Sequential Trial at the ",object.type.k," stage. \n", sep = "")
        }else{
            cat("Group Sequential Trial at the ",object.type.k," analysis of stage ",object.k,". \n", sep = "")
        }
        if(type.update=="normal"){
            cat("Update for the ",type.k," analysis of stage ",k,". \n", sep = "")
        }else if(type.update=="information"){
            cat("Update information and bound relative to the (skipped) decision analysis of stage ",k,". \n", sep = "")
        }
    }
    
    #browser()

    ## ** update object with information and estimate
    ## skipped
    if(type.update=="information"){
        if(inherits(delta,"lmmGSD")){
            object$Info.d[k] <- as.double(delta$information["decision"])
            object$n.obs[k] <- delta$sample.size["total"]
        }else{
            object$Info.d[k] <- as.double(Info.d)
        }
        
        object <- updateBoundaries(object, k = k, type.k = type.k, trace = trace-1, update.stage = FALSE)
        return(object)
    }

    ## normal
    if(inherits(delta,"lmmGSD")){
        if(trace>0){cat(" - extract information from mixed model: ", sep = "")}
        if(type.k %in% c("interim","final")){
            object$lmm[[k]] <- delta
        }else{
            object$lmm[[k+1]] <- delta
        }
        if(type.k %in% c("interim","final")){
            Info.i <- delta$information["interim"]
        }
        if(type.k %in% c("interim","decision")){
            Info.d <- delta$information["decision"]
        }
        if(is.null(delta$data.decision)){
            object$n.obs[k] <- delta$sample.size["total"]
        }else{
            object$n.obs[k] <- delta$sample.size["decision"]
        }
        delta <- delta$delta
        if(trace>0){cat("done \n", sep = "")}
    }else{
        if(is.numeric(delta) && length(delta)==1){
            delta <- data.frame(estimate = delta, se = NA, statistic = NA, df = NA, p.value = NA)
            if(type.k=="interim"){
                delta$se <- 1/sqrt(Info.i)
            }else if(type.k=="decision"){
                delta$se <- 1/sqrt(Info.d)
            }else if(type.k=="final"){
                if(missing(Info.d) && !missing(Info.i)){
                    delta$se <- 1/sqrt(Info.i)
                }else if(missing(Info.i) && !missing(Info.d)){
                    delta$se <- 1/sqrt(Info.d)                    
                }else{
                    stop("Could not guess the standard error based on the provided information. \n",
                         "At final, exactly only one of the two arguments: \'Info.d\' and  \'Info.i\' should be provided. \n")
                }
            }
            delta$statistic <- delta$estimate / delta$se
            delta$df <- Inf
            delta$p.value <- 2*(1-pnorm(abs(delta$statistic))) ## ignore object$alternative
        }else if(!is.data.frame(delta) || NROW(delta)!=1 || NCOL(delta)!=5){
            stop("Argument \'delta\' should be a data.frame with 1 line and 5 columns \n.")
        }else if(any(names(delta) != c("estimate","se","statistic","df","p.value"))){
            stop("Column names in argument \'delta\' should be \"estimate\", \"se\", \"statistic\", \"df\", \"p.value\" \n")
        }
        
    }
    object$delta <- rbind(object$delta,
                          data.frame(method = "ML", stage = k, type = type.k,
                                     delta,
                                     lower = delta$estimate + stats::qt((1-object$conf.level)/2, df = delta$df) * delta$se,
                                     upper = delta$estimate + stats::qt(1-(1-object$conf.level)/2, df = delta$df) * delta$se))
    if(type.k == "interim"){
        object$Info.i[k] <- as.double(Info.i)
    }else if(type.k == "final"){
        if(missing(Info.d) && !missing(Info.i)){
            object$Info.d[k] <- as.double(Info.i)
        }else if(missing(Info.i) && !missing(Info.d)){
            object$Info.d[k] <- as.double(Info.d)
        }else{
            stop("Could not guess the information based on the provided information. \n",
                 "At final, exactly only one of the two arguments: \'Info.d\' and  \'Info.i\' should be provided. \n")
        }
    }
    if(type.k %in% "interim" && !missing(Info.d)){ ## predicted information at decision is used by certain methods
        object$Info.d[k] <- as.double(Info.d)
    }else if(type.k %in% "decision"){
        object$Info.d[k] <- as.double(Info.d)
    }

    ## ** update boundaries
    ## also update k and type.k
    if(trace>0){cat(" - update boundaries: ", sep = "")}
    object <- updateBoundaries(object, k = k, type.k = type.k, trace = trace-1, update.stage = TRUE)
    if(trace>0){cat("done \n", sep = "")}

    ## ** decision
    if(trace>0){cat(" - update decision: ", sep = "")}
    if(type.k != "interim" || is.na(object$conclusion["interim",k])){
            object <- Decision(object, k = k, type.k = type.k, trace = trace-1)
    }
    
    if(trace>0){cat("done \n", sep = "")}

    ## ** estimate
    if(type.k %in% c("decision","final")){
        if(trace>0){cat(" - correct estimate: ", sep = "")}

        ## read default settings
        options <- DelayedGSD.options()
        continuity.correction <- options$continuity.correction
        tolerance <- options$tolerance
        FCT.p_value <- options$FCT.p_value

        ## prepare
        delta.MUE <- data.frame(method = "MUE", stage = k, type = type.k,
                                estimate = NA, se = NA, statistic = NA, df = NA, lower = NA, upper = NA, p.value = NA)

        delta <- stats::confint(object)
        Info.i <- object$Info.i
        Info.d <- object$Info.d
        ck <- object$ck
        ck.unrestricted <- object$ck.unrestricted
        lk <- object$lk
        uk <- object$uk

        ## arguments
        ls.args <- list(Info.d = Info.d[1:k],  
                        Info.i = Info.i[1:min(k,kMax-1)],
                        ck = ck[1:min(k,kMax)],
                        ck.unrestricted = ck.unrestricted[1:min(k,kMax)],
                        lk = lk[1:min(k,kMax-1)],  
                        uk = uk[1:min(k,kMax-1)],
                        reason.interim = object$conclusion["reason.interim",1:k],
                        kMax = kMax, 
                        estimate = delta[1,"estimate"],
                        method = object$method,
                        bindingFutility = object$bindingFutility,
                        cNotBelowFixedc=object$cNotBelowFixedc,
                        continuity.correction=continuity.correction)

        ## *** p.value
        if(p.value || ci){
            resP <- do.call(FCT.p_value, c(ls.args, list(delta=0))) 
            if(p.value){
                delta.MUE[1,"p.value"] <- as.double(resP)
                attr(delta.MUE,"error") <- c(p.value = unname(attr(resP,"error")))
            }
        }
            
        ## *** CI
        if(ci & as.double(resP) < 1){

            resCI <- do.call(FinalCI, c(ls.args, list(FCT.p_value = FCT.p_value, conf.level = object$conf.level, tolerance = tolerance, conclusion = object$conclusion["decision",k])))
            delta.MUE[1,"lower"] <- resCI["lower"]
            delta.MUE[1,"upper"] <- resCI["upper"]

            if(is.null(attr(delta.MUE,"error"))){
                attr(delta.MUE,"error") <- c(lower = unname(attr(resCI,"error")["lower"]), upper = unname(attr(resCI,"error")["upper"]))
            }else{
                attr(delta.MUE,"error") <- c(attr(delta.MUE,"error"), lower = unname(attr(resCI,"error")["lower"]), upper = unname(attr(resCI,"error")["upper"]))
            }
        }
        
        ## *** Estimate
        if(estimate){
            resMUE <- do.call(FinalEstimate, c(ls.args, list(FCT.p_value = FCT.p_value, tolerance = tolerance)))

            delta.MUE[1,"estimate"] <- resMUE
            if(is.null(attr(delta.MUE,"error"))){
                attr(delta.MUE,"error") <- c(estimate = unname(attr(resMUE,"error")))
            }else{
                attr(delta.MUE,"error") <- c(attr(delta.MUE,"error"), estimate = unname(attr(resMUE,"error")))
            }
        }
        object$delta <- rbind(object$delta, delta.MUE)
        attr(object$delta,"error") <- attr(delta.MUE,"error")
        if(trace>0){cat("done \n", sep = "")}
    }

    ## ** export
    return(object)
}

## * .getStage
## Identify the current stage or check that the values of k and type.k are compatible with the object
.getStage <- function(object.stage, object.conclusion, kMax, k, type.k, nextStage){

    ## ** extract information
    object.k <- object.stage$k
    object.type.k <- object.stage$type
    if(nextStage){
        type.update <- "normal"
    }else{
        type.update <- NULL
    }

    ## ** basic checks
    if(!is.null(type.k)){
        type.k <- match.arg(type.k, c("interim","decision","final"))
    }
    if(!is.null(k)){
        if(k<0 || k>kMax || k %% 1 != 0){
            stop("Argument \'k\' must be an integer between 1 and ",kMax,"")
        }
    }

    ## ** more precise check
    if(object.k==0){

        ## *** planning
        if(nextStage){
            if(!is.null(k)){
                if(k!=1){
                    stop("Argument \'k\' should be 1 just after the planning stage. \n")
                }
            }else{ 
                k <- 1
            }
            if(!is.null(type.k)){
                if(type.k!="interim"){
                    stop("Argument \'type.k\' should be \"interim\" just after the planning stage. \n")
                }
            }else{
                type.k <- "interim"
            }
        }else{
            if(!is.null(k)){
                if(k!=0){
                    stop("Argument \'k\' should be 0 at the planning stage. \n")
                }
            }else{ 
                k <- 0
            }
            if(!is.null(type.k)){
                if(type.k!="planning"){
                    stop("Argument \'type.k\' should be \"planning\" just at the planning stage. \n")
                }
            }else{
                type.k <- "planning"
            }
        }
    }else if(object.type.k=="interim"){

        if(nextStage){
            object.conclusionInterim <- object.conclusion["interim",object.k]

            if(object.conclusionInterim=="continue"){ ## at interim where we conclude to continue
                if(!is.null(k)){
                    if((k==object.k) && (type.k == "decision")){
                        type.update <- "information"
                    }else if(k!=(object.k+1)){
                        stop("Argument \'k\' should be ",object.k+1," after continuing recruitment following the interim of stage ",object.k,". \n")
                    }
                }else{
                    k <- object.k + 1
                }
                if(k==kMax){
                    if(!is.null(type.k)){
                        if(type.k!="final"){
                            stop("Argument \'type.k\' should be \"final\" after continuing recruitment following the interim of stage ",object.k,". \n")
                        }
                    }else{
                        type.k <- "final"
                    }
                }else{ ## k < kMax
                    if(!is.null(type.k)){
                        if((k==object.k) && (type.k == "decision")){
                            ## ok
                        }else if(type.k!="interim"){
                            stop("Argument \'type.k\' should be \"interim\" after continuing recruitment following the interim of stage ",object.k,". \n")
                        }
                    }else{
                        type.k <- "interim"
                    }
                }
            }else if(object.conclusionInterim=="stop"){ ## at interim where we conclude to stop
                if(!is.null(k)){
                    if(k!=object.k){
                        stop("Argument \'k\' should be ",object.k," after stopping recruitment following the interim of stage ",object.k,". \n")
                    }
                }else{
                    k <- object.k
                }
                if(!is.null(type.k)){
                    if(type.k!="decision"){
                        stop("Argument \'type.k\' should be \"decision\" after stopping recruitment following the interim of stage ",object.k,". \n")
                    }
                }else{
                    type.k <- "decision"
                }
            }

        }else{
            if(!is.null(k)){
                if(k %in% 1:object.k == FALSE){
                    stop("Argument \'k\' should be an integer between 1 and ",object.k,". \n")
                }
            }else{
                k <- object.k
            }
            if(!is.null(type.k)){
                if(type.k %in% "interim" == FALSE){
                    stop("Argument \'type.k\' can only be \"interim\". \n")
                }
            }else{
                type.k <- "interim"
            }
        }

    }else if(object.type.k=="decision"){
        if(nextStage){
            stop("The decision analysis has already been reached, no further analysis. \n")
        }
        if(!is.null(k)){
            if(k %in% 1:object.k == FALSE){
                stop("Argument \'k\' should be an integer between 1 and ",object.k,". \n")
            }
        }else{
            k <- object.k
        }
        if(!is.null(type.k)){
            if(k<object.k && type.k=="decision"){
                stop("Argument \'type.k\' cannot be \"decision\" as recruitement continued after interim ",k,". \n")
            }else if(type.k %in% c("decision","interim") == FALSE){
                stop("Argument \'type.k\' can either be \"interim\" or \"decision\". \n")
            }
        }else{
            type.k <- "decision"
        }

    }else if(object.type.k=="final"){
    
        if(nextStage){
            stop("The final analysis has already been reachedm no further analysis. \n")
        }
        if(!is.null(k)){
            if(k %in% 1:object.k == FALSE){
                stop("Argument \'k\' should be an integer between 1 and ",object.k,". \n")
            }
        }else{
            k <- object.k
        }
        if(!is.null(type.k)){
            if(k==object.k){
                if(type.k!="final"){
                    stop("Argument \'type.k\' must be \"final\" at stage ",k,". \n")
                }
            }else if(type.k %in% c("decision","interim") == FALSE){
                if(k>2){
                    stop("Argument \'type.k\' can either be \"interim\" or \"decision\" between stage 1 and ",k-1,". \n")
                }else{
                    stop("Argument \'type.k\' can either be \"interim\" or \"decision\" at stage 1. \n")
                }
            }
        }else{
            type.k <- "final"
        }

    }

    ## ** export
    return(list(k=k,
                type.k = type.k,
                type.update = type.update))
}

## * .overrule
.overrule <- function(object){
    kMax <- object$kMax
    object.type.k <- object$stage$type
    object.k <- object$stage$k

    if(object.type.k != "interim"){
        stop("Can only overrule a futility boundary at an interim analysis. \n")
    }
    if(object$conclusion["interim",object.k]=="continue"){
        message("No need to overrule the futility bound as it has been decided to continue. \n")
    }else{
        object$conclusion["interim",object.k] <- "continue"
        object$conclusion["reason.interim",object.k] <- "overrule futility"
    }

    return(object)
}
