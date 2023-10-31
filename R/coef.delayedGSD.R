## * coef.delayedGSD (documentation)
#' @title Treatment Effect Estimates for a Group Sequential Design with Delayed Endpoints
#' @description Extract estimate relative to the treatment effect at a specific stage of a group sequential design with delayed endpoints.
#' By default extract value for the latest stage that has been performed.
#' 
#' @param object object of class \code{delayedGSD}, typically output from \code{\link{CalcBoundaries}}.
#' @param type [character] Shoudl the estimate effect (\code{"effect"}), boundaries (\code{"boundaries"}), information  (\code{"information"}), or decision  (\code{"deicsion"})
#' be output. The estimate is only displayed for the current stage while the other are displayed for all stages.
#' @param planned [logical] Should the planned or estimated boundaries/information be output.
#' @param predicted [logical] Should the predicted information/boundaries at decision based on interim data be output (when relevant).
#' @param k [integer] stage relative to which the estimates should be output.
#' @param type.k [character] type of stage relative to which the estimates should be output: \code{"interim"}, \code{"decision"}, or \code{"final"}.
#' @param method [character] type of estimate to output: can be \code{"ML"} or  \code{"corrected ML"}, the latter accounting for the group sequential design.
#' @param ... not used, for compatibility with the generic method.
#' 

## * coef.delayedGSD
##' @export
coef.delayedGSD <- function(object, type = "effect", planned = NULL, predicted = TRUE, k = NULL, type.k = NULL, method = NULL, ...){

    ## ** extract information from object
    type <- match.arg(type, c("effect","boundary","information","decision"))
    kMax <- object$kMax
    resStage <- .getStage(object.stage = object$stage,
                          object.conclusion = object$conclusion,
                          kMax = kMax,
                          k = k,
                          type.k = type.k,
                          nextStage = FALSE)
    if(is.null(k)){
        k <- resStage$k
    }
    if(is.null(type.k)){
        type.k <- resStage$type.k
    }

    if(is.null(planned)){
        planned <- (type.k=="planning")
    }else if(planned %in% c(TRUE,FALSE) == FALSE){
        stop("Argument \'planned\' should be TRUE or FALSE. \n")
    }

    Info.max <- object$planned$Info.max

    #browser()
    
    ## ** retrive quantities
    if(type=="effect"){
        if(planned==TRUE){
            out <- object$planned$delta
        }else{
            outConfint <- stats::confint(object, method = method, k = k, type.k = type.k)
            out <- stats::setNames(outConfint$estimate,outConfint$method)
        }
    }else if(type=="boundary"){
        if(planned==TRUE){

            out <- data.frame(stage = 1:kMax,
                              Fbound = c(object$planned$lk[1:(kMax-1)],NA),
                              Ebound = c(object$planned$uk[1:(kMax-1)],NA),
                              Cbound = object$planned$ck)

        }else{

            statistic <- unname(object$delta$statistic)
            statistic.interim <- c(statistic[0:k], rep(NA,kMax-k))
            statistic.interim[kMax] <- NA ## at final analysis no interim only decision
            statistic.decision <- rep(NA,kMax)
            if(type.k=="decision"){
                statistic.decision[k] <- statistic[k+1]
            }else if(type.k=="final"){
                statistic.decision[k] <- statistic[k]
            }

            ## assemble
            out <- data.frame(stage = 1:kMax, Fbound = c(object$lk[1:(kMax-1)],NA), Ebound = c(object$uk[1:(kMax-1)],NA), statistic.interim = statistic.interim,
                              Cbound = object$ck, statistic.decision = statistic.decision)
            rownames(out) <- NULL                

            ## remove critical boundaries when continuing at interim
            if(type.k == "interim" && predicted==FALSE){
                out$Cbound[k] <- NA
            }        
                
                
        }

    }else if(type=="information"){

        if(planned==TRUE){
            Info.i <- c(object$planned$Info.i,utils::tail(object$planned$Info.d,1)) ## interim and final
            Info.d <- object$planned$Info.d ## decision and final
            InfoR.i <- Info.i/Info.max
            InfoR.d <- Info.d/Info.max
            if(!is.null(object$n.obs)){
                n.obs <- object$n.obs*InfoR.d
            }else{
                n.obs <- NA
            }
        }else{
          #browser()
            Info.i <- c(object$Info.i, utils::tail(object$Info.d,1)) ## interim and final
            Info.d <- object$Info.d ## decision and final
            Info <- c(object$planned$Info.i,object$planned$Info.d[object$kMax])
            InfoR.i <- Info.i/Info.max
            InfoR.d <- Info.d/Info.max
            if(type.k=="decision"){
                index.lmm <- setdiff(1:(k+1),k)
            }else{
                index.lmm <- 1:k
            }            
            n.obs <- object$n.obs
            #browser()
            index.lastNNA <- utils::tail(which(!is.na(n.obs)),1)
            index.NA <- which(is.na(n.obs))
            if(length(index.lastNNA)>0 && length(index.NA)>0 && predicted && type.k == "interim" && Info[index.NA[1]]>Info.i[index.lastNNA]){
                n.obs[index.NA] <- n.obs[index.lastNNA] * Info[index.NA]/Info.i[index.lastNNA]
            }
            if(type.k == "interim" && predicted==FALSE){
                Info.d[k] <- NA
                InfoR.d[k] <- NA
            }
        }
        out <- data.frame(stage = 1:kMax,
                          Interim = Info.i,
                          Interim.pc = InfoR.i,
                          Decision = Info.d,
                          Decision.pc = InfoR.d,
                          n = n.obs)
        rownames(out) <- NULL
        attr(out,"Info.max") <- Info.max
    }else if(type == "decision"){

        if(resStage$k==0){
            return(NULL)
        }
        
        out <- object$conclusion[c("interim","reason.interim"),1:resStage$k,drop=FALSE]
        rownames(out) <- c("decision","comment")
        colnames(out) <- paste0("stage ",1:resStage$k)

        if(resStage$type.k == "decision"){
            index.decision <- which(!is.na(object$conclusion["decision",]))
            add.decision <- object$conclusion[c("decision","comment.decision"),index.decision,drop=FALSE]
            colnames(add.decision) <- paste0("stage ",resStage$k," decision")
            out <- cbind(out, add.decision)
        }else if(resStage$type.k == "final"){
            out[,paste0("stage ",kMax)] <- object$conclusion[c("decision","comment.decision"),kMax,drop=FALSE]            
        }
    }

    ## ** export
    return(out)
}
