## * confint.delayedGSD (documentation)
#' @title Treatment Effect Estimates for a Group Sequential Design with Delayed Endpoints
#' @description Extract estimate, confidence intervals, p-value, ... relative to the treatment effect at a specific stage of a group sequential design with delayed endpoints.
#' By default extract value for the latest stage that has been performed.
#' 
#' @param object object of class \code{delayedGSD}, typically output from \code{\link{CalcBoundaries}}.
#' @param parm not used, for compatibility with the generic method.
#' @param level not used, for compatibility with the generic method.
#' @param method [character] type of estimate to output: can be \code{"ML"} or  \code{"MUE"}, the latter accounting for the group sequential design (median unbiased estimates).
#' @param k [integer] stage relative to which the estimates should be output.
#' @param type.k [character] type  of stage relative to which the estimates should be output: \code{"interim"}, \code{"decision"}, or \code{"final"}.
#' @param ... not used, for compatibility with the generic method.
#' 

## * confint.delayedGSD (code)
#' @export
confint.delayedGSD <- function(object, parm = NULL, level = NULL, method = NULL, k = NULL, type.k = NULL, ...){


    ## ** check user input
    if(object$stage$type=="planning"){return(NULL)} ## no estimate at the planning stage
    if(!is.null(parm)){
        stop("Argument \'parm\' not used. ")
    }
    if(!is.null(level)){
        stop("Argument \'level\' not used. ")
    }
    dots <- list(...)
    if(length(dots)){
        stop("Argument(s) \'",paste(names(dots), collapse = "\' \'"),"\' not used. ")
    }

    if(identical(k,"all")){

        out <- object$delta

    }else{
    
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
    
        if(!is.null(method)){
            method <- toupper(method)
            if(any(method %in% c("ML","MUE") == FALSE)){
                stop("Incorrect argument \'method\': should take value \"ML\" and/or \"MUE\". \n")
            }
            ## if("MUE" %in% method &&  type.k == "interim"){
            ##     stop("Argument \'method\' cannot be \"MUE\" at interim. \n")
            ## }
        }

        ## ** extract information from object
        out <- object$delta[object$delta$type == type.k & object$delta$stage == k,]
        
    }



    ## ** export
    if(any(sapply(object$lmm, length)>0)){
        out$coef <- unlist(sapply(object$lmm,"[[","name.coef"))[1]
    }else{
        out$coef <- NA
    }
    out <- out[,c("method","stage","type","coef","estimate","se","statistic","df","p.value","lower","upper")]
    if(!is.null(method)){
        out <- out[out$method %in% method,,drop=FALSE]
    }
    rownames(out) <- NULL
    if("MUE" %in% method){
        attr(out,"error") <- attr(object$delta,"error")
        test <- na.omit(attr(out,"error"))
        if(!is.null(test) && any(abs(test)>1e-3)){
            warning("Possibly incorrect evaluation of the MUE ",paste(names(test[abs(test)>1e-3]),collapse = ", "),".\n",
                    "Mismatch in the optimization process in term of confidence level: ",paste(test[abs(test)>1e-3],collapse = ", "),". \n",sep="")
        }
    }

    return(out)

}
