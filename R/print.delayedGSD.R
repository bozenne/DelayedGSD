## * print.delayedGSD (documentation)
#' @title Status of a Group Sequential Design with Delayed Endpoints
#' @description Display boundaries and estimated treatment effect up the current stage.
#' 
#' @param x object of class \code{delayedGSD}, typically output from \code{\link{CalcBoundaries}}.
#' @param ... arguments passed to the summary function.
#' 

## * print.delayedGSD
#' @export
print.delayedGSD <- function(x, ...){

    summary(x, space = NULL, abreviated = 2, ...)
}

