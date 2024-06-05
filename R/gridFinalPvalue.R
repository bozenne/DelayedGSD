### plotFinalPvalue.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 21 2023 (09:42) 
## Version: 
## Last-Updated: jun  5 2024 (14:25) 
##           By: Brice Ozenne
##     Update #: 140
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * gridFinalPvalue (documentation)
##' @title P-value over the Ordering Space
##' 
##' @description Evaluate the p-value over stages and test statistics, following armitage ordering (Armitage 1957).
##' 
##' @param object 
##' @inheritParams FinalPvalue
##' 
##' @return a \code{data.frame} containing \itemize{
##' \item \code{k} the stage
##' \item \code{z} the test statistic value
##' \item \code{p.value} the p-value
##' }
##'
##' @references P. Armitage, Restricted sequential procedures. Biometrika 1957 (9-56)
##' 
##' @examples
##' ## 2 stages, no fixC
##' myBound2 <- CalcBoundaries(kMax = 2,  
##'                          alpha = 0.025, 
##'                          beta = 0.2,  
##'                          InfoR.i = c(0.55, 1.00),  
##'                          InfoR.d = c(0.62, 1),  
##'                          rho_alpha = 2,  
##'                          rho_beta = 2,  
##'                          method = 1,  
##'                          cNotBelowFixedc = FALSE,
##'                          bindingFutility= TRUE,
##'                          delta = 0.6)
##' grid.p2 <- gridFinalPvalue(myBound2)
##' grid.p2
##' 
##' ## 2 stages, fixC
##' myBound2.fixC <- CalcBoundaries(kMax = 2,  
##'                          alpha = 0.025, 
##'                          beta = 0.2,  
##'                          InfoR.i = c(0.55, 1.00),  
##'                          InfoR.d = c(0.62, 1),  
##'                          rho_alpha = 2,  
##'                          rho_beta = 2,  
##'                          method = 1,  
##'                          cNotBelowFixedc = TRUE,
##'                          bindingFutility= TRUE,
##'                          delta = 0.6)
##' gridFinalPvalue(myBound2.fixC,
##'                 continuity.correction = 0)
##' gridFinalPvalue(myBound2.fixC,
##'                 continuity.correction = 1)
##' gridFinalPvalue(myBound2.fixC,
##'                 continuity.correction = 2)
##'
##' ## 3 stages, no fixC
##' myBound3 <- CalcBoundaries(kMax = 3,  
##'                          alpha = 0.025, 
##'                          beta = 0.2,  
##'                          InfoR.i = c(0.45, 0.65, 1.00),  
##'                          InfoR.d = c(0.52, 0.72, 1),  
##'                          rho_alpha = 2,  
##'                          rho_beta = 2,  
##'                          method = 1,  
##'                          cNotBelowFixedc = FALSE,
##'                          bindingFutility= TRUE,
##'                          delta = 0.6)
##' grid.p3 <- gridFinalPvalue(myBound3)
##' grid.p3
##' 
##' ## 3 stages, fixC
##' myBound3.fixC <- CalcBoundaries(kMax = 3,  
##'                          alpha = 0.025, 
##'                          beta = 0.2,  
##'                          InfoR.i = c(0.45, 0.65, 1.00),  
##'                          InfoR.d = c(0.52, 0.72, 1),  
##'                          rho_alpha = 2,  
##'                          rho_beta = 2,  
##'                          method = 1,  
##'                          cNotBelowFixedc = TRUE,
##'                          bindingFutility= TRUE,
##'                          delta = 0.6)
##' gridFinalPvalue(myBound3.fixC, digits = 3,
##'                 continuity.correction = 0)
##' gridFinalPvalue(myBound3.fixC, digits = 3,
##'                 continuity.correction = 1)
##' gridFinalPvalue(myBound3.fixC, digits = 3,
##'                 continuity.correction = 2)
##'

## * gridFinalPvalue (code)
##' @export
gridFinalPvalue <- function(object,
                            continuity.correction,
                            seq.futility = c(-3,0),
                            seq.kMax = c(-3,1,3,6),
                            seq.efficacy = c(2.5,3.25,4,5)){

    ## ** normalize arguments
    options <- DelayedGSD.options()
    FCT.p_value <- options$FCT.p_value
    seq.futility <- sort(seq.futility)
    seq.kMax <- sort(seq.kMax)
    seq.efficacy <- sort(seq.efficacy)

    #browser()
    
    if(inherits(object,"delayedGSD")){

        kMax <- object$kMax
        method <- object$method
        bindingFutility <- object$bindingFutility
        cNotBelowFixedc <- object$cNotBelowFixedc

        if(object$stage$type=="planning"){
            kCurrent <- kMax
            Info.d <- object$planned$Info.d
            Info.i <- object$planned$Info.i
            ck <- object$planned$ck
            ck.unrestricted <- object$planned$ck.unrestricted
            lk <- object$planned$lk
            uk <- object$planned$uk
            reason.interim <- rep("",kMax)
            alphaSpent <- object$planned$alphaSpent
        }else{
            kCurrent <- object$stage["k"]
            if(kCurrent!=kMax){
                stop("Current stage (k=",kCurrent,") is not kMax=",kMax,". \n",
                     "Can only display p-values over the whole space. \n")
            }
            Info.d <- object$Info.d
            Info.i <- object$Info.i
            ck <- object$ck
            ck.unrestricted <- object$ck.unrestricted
            lk <- object$lk
            uk <- object$uk
            reason.interim <- object$conclusion["reason.interim",]
            alphaSpent <- object$alphaSpent
        }
    }else if(is.list(object)){

        valid.args <- c("Info.d", "Info.i", "ck", "ck.unrestricted", "lk", "uk",
                        "reason.interim", "kMax", "method", "bindingFutility", "cNotBelowFixedc")
        if(any(names(object) %in% valid.args == FALSE)){
            stop("Unknown argument(s) \"",paste(names(object)[names(object) %in% valid.args == FALSE], collapse = "\", \""),"\".\n")
        }
        if(any(valid.args %in% names(object) == FALSE)){
            stop("Missing argument(s) \"",paste(valid.args[valid.args %in% names(object) == FALSE], collapse = "\", \""),"\".\n")
        }

        for(iArg in 1:length(valid.args)){
            assign(valid.args[iArg], value = object[[valid.args[iArg]]], envir = environment())
        }
        alphaSpent <- rep(NA, kMax)
    }else{
        stop("Argument \'object\' should either be \n",
             " a delayedGSD object (output of calcBoundaries functions) \n",
             " a list containing the arguments of the function finalPvalue. \n")
    }
    if(kMax==1){
        stop("kMax should be strictly greater than 1. \n")
    }

    if(missing(continuity.correction)){
        continuity.correction <- options$continuity.correction
    }

    ## ** grid of values
    grid <- NULL
    for(k in 1:(kMax-1)){
        grid <- rbind(grid,
                      data.frame(name = paste0("fut.",k,".",1:sum(seq.futility<ck.unrestricted[k])),
                                 z = seq.futility[seq.futility<ck.unrestricted[k]],
                                 k = k,
                                 fixC = FALSE,
                                 alphaSpent = NA)
                      )

        if(cNotBelowFixedc & ck.unrestricted[k]!=ck[k]){
            grid <- rbind(grid,
                          data.frame(name = paste0("fut.",k,".Bl"), z = ck.unrestricted[k], k = k, fixC = TRUE, alphaSpent = NA),
                          data.frame(name = paste0("fut.",k,".Bu"), z = 0.1*ck.unrestricted[k]+0.9*ck[k], k = k, fixC = TRUE, alphaSpent = NA))
        }else{
            grid <- rbind(grid,
                          data.frame(name = paste0("fut.",k,".B"), z = ck[k]-1e-5, k = k, fixC = FALSE, alphaSpent = NA))
        }
    }

    grid <- rbind(grid,
                  data.frame(name = paste0("fut.",kMax,".",1:sum(seq.kMax<ck[kMax])),
                             z = seq.kMax[seq.kMax<ck[kMax]],
                             k = kMax,
                             fixC = FALSE,
                             alphaSpent = NA),
                  data.frame(name = paste0("eff.",kMax,".B"),
                             z = ck[kMax],
                             k = kMax,
                             fixC = FALSE,
                             alphaSpent = alphaSpent[kMax]),
                  data.frame(name = paste0("eff.",kMax,".",1:sum(seq.kMax>ck[kMax])),
                             z = seq.kMax[seq.kMax>ck[kMax]],
                             k = kMax,
                             fixC = FALSE,
                             alphaSpent = NA)
                  )
    for(k in (kMax-1):1){
        grid <- rbind(grid,
                      data.frame(name = paste0("eff.",k,".B"),
                                 z = ck[k]+1e-5,
                                 k = k,
                                 fixC = FALSE,
                                 alphaSpent = alphaSpent[k]),
                      data.frame(name = paste0("eff.",k,".",1:sum(seq.efficacy>ck[k])),
                                 z = seq.efficacy[seq.efficacy>ck[k]],
                                 k = k,
                                 fixC = FALSE,
                                 alphaSpent = NA)
                      )
        
    }
    n.grid <- NROW(grid)

    ## ** evaluate p-value over the domain
    calcP <- function(z, k){
        outP <- do.call(FCT.p_value, list(Info.d = Info.d[1:k],
                                         Info.i = Info.i[1:min(k,kMax-1)],
                                         ck = ck[1:min(k,kMax)],
                                         ck.unrestricted = ck.unrestricted[1:min(k,kMax)],
                                         lk = lk[1:min(k,kMax-1)],
                                         uk = uk[1:min(k,kMax-1)],
                                         kMax = kMax,
                                         delta = 0, 
                                         statistic = z,
                                         reason.interim = reason.interim[1:k],
                                         method = method,
                                         bindingFutility = bindingFutility, 
                                         cNotBelowFixedc = cNotBelowFixedc,
                                         continuity.correction = continuity.correction)
                        )
        out <- data.frame(z = z, k = k, p.value = outP, correction = continuity.correction)
        attr(out,"terms") <- attr(outP,"terms")
        return(out)
    }
    
    ls.out <- lapply(1:n.grid, function(iG){ ## iG <- 1
        outP <- calcP(z = grid[iG,"z"], k = grid[iG,"k"])
        out <- cbind(outP, fixC = grid[iG,"fixC"], alphaSpent = grid[iG,"alphaSpent"])
        attr(out,"terms") <- attr(outP,"terms")
        return(out)
    })
    out <- do.call(rbind,ls.out)
    attr(out,"terms") <- stats::setNames(lapply(ls.out, attr, "terms"), paste0("stage=",out$k,": z=",out$z))

    ## ** export
    attr(out,"class") <- append("gridDelayedGSD",attr(out,"class"))
    return(out)
}


##----------------------------------------------------------------------
### plotFinalPvalue.R ends here
