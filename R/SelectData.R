#' @title Select available longitudnal data for analysis at a specific follow-up time.
#' @description yy
#' 
#' @param d data set generated with GenData function.
#' @param t time at which we want the available data.
#' 
#' @details xxx
#' @return input dataset  with less rows and NA where appropriate.
#' @author Paul Blanche
#' 
#' @examples
#' x <- GenData(n=50)
#' head(x$d,n=20)
#' tail(x$d)
#' y <- SelectData(x$d,t=4)
#' head(y)
#' tail(y)
#' PlotProgress(x$d,at=4)
#'
#' #----- select when half of the subjects have one follow-up measuement---
#' # which depends on accrual rate (ar) 
#'
#' x <- GenData(n=35)
#' thear <- 10
#' thet <- x$d$t2[ceiling(nrow(x$d)/2) + ceiling(thear)]
#' PlotProgress(x$d,at=thet)
#' y <- SelectData(x$d,t=thet)
#' y
#'
#' @export
SelectData <- function(d,t){

    ## find columns in the dataset d corresponding to the outcome and the time
    col.outcome <- grep("X",names(d), value = TRUE)
    col.times <- grep("t",names(d), value = TRUE)
    n.times <- length(col.times)

    ## add indicator of missing values
    col.missing <- gsub("^t","missing",col.times)
    test.na <- 1.0 * is.na(d[,col.outcome,drop=FALSE]) ## 1.0 * to convert to numeric
    colnames(test.na) <- col.missing
    d.augm <- cbind(d, test.na)

    ## subset according to time
    dd <- d.augm[which(d.augm$t1 <= t),,drop=FALSE]

    ## add missing value according to time as the baseline can be observed but some of the follow-up value may be not
    for(j in 2:n.times){
        whereNA <- which(dd[,paste0("t",j)] > t)
        dd[whereNA,col.missing[[j]]] <- -1
        dd[whereNA,paste0("t",j)] <- NA
        dd[whereNA,paste0("X",j)] <- NA
    }
    ## dd$missing <- any(dd[,paste0("missing",1:n.times)]==1)

    ## keep arguments
    attr(dd,"t") <- t
    return(dd)
}
