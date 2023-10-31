#' @title Student's t-Test
#' @description Performs one and two sample t-tests on vectors of data. 
#' The only modification to the \code{stats::t.test} is that it stores the call and the arguments.
#'
#' @param ... arguments passed to \code{stats::t.test}
#'
#' @export
ttest <- function(...){
    out <- stats::t.test(...)

    ## ** add call and arguments to the output
    out$call <- match.call()
    out$args <- list(...)

    ## ** properly name arguments
    if("formula" %in% names(out$args) == FALSE & "x" %in% names(out$args) == FALSE){
        if(is.null(names(out$args))){
            if(inherits(out$args[[1]],"formula")){
                names(out$args)[1] <- "formula"
            }else{
                names(out$args)[1] <- "x"
            }
        }else{
            index.first <- which(names(out$args[1])=="")[1]
            if(inherits(out$args[[index.first]],"formula")){
                names(out$args)[index.first] <- "formula"
            }else{
                names(out$args)[index.first] <- "x"
            }
        }
    }

    if("formula" %in% names(out$args) & "data" %in% names(out$args) == FALSE){
        names(out$args)[min(c(which(names(out$args)==""),which(is.na(names(out$args)))))] <- "data"
    }
    if("x" %in% names(out$args) & out$method != "One Sample t-test" & "y" %in% names(out$args) == FALSE){
        names(out$args)[min(c(which(names(out$args)==""),which(is.na(names(out$args)))))] <- "y"
    }

    ## ** convert all t-test to (x,y)
    if("formula" %in% names(out$args)){
        name.outcome <- all.vars(out$args[["formula"]])[1]
        if(out$method == "One Sample t-test"){
            out$args$x <- out$args[["data"]][[name.outcome]]
        }else{
            name.group <- all.vars(out$args[["formula"]])[2]
            out$args[c("x","y")]  <- tapply(out$args[["data"]][[name.outcome]],
                                            out$args[["data"]][[name.group]],
                                            function(iVec){list(iVec)})
        }
    }

    ## ** set class
    class(out) <- append("ttest",class(out))
    return(out)
}
