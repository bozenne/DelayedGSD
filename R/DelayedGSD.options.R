### DelayedGSD.options.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 11 2023 (11:47) 
## Version: 
## Last-Updated: okt 11 2023 (15:41) 
##           By: Brice Ozenne
##     Update #: 15
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * DelayedGSD.options (documentation) 
#' @title Global options for DelayedGSD package
#' @include 0-onload.R
#'
#' @description Update or select global options for the DelayedGSD package.
#'
#' @param ... options to be selected or updated
#' @param reinitialise should all the global parameters be set to their default value
#'
#' @details The options are: \itemize{
#' \item FCT.p_value [character]: function to evaluate p-values (FinalPvalue or FinalPvalue2)
#' \item continuity.correction [0,1,2]: correction used to evaluate the p-value when using \code{cNotBelowFixedc=TRUE}. Positive values ensures continuity of the p-value across stages.
#' \item tolerance [numeric,>0]: acceptable discrepancy to the objective level when evaluating the confidence intervals and median unbiased estimate.
#' }
#'
#' @return A list containing the default options.
#'
#' @keywords utilities
#' @examples
#' DelayedGSD.options()


 
## * DelayedGSD.options (code)
#' @export
DelayedGSD.options <- function(..., reinitialise = FALSE){

    default <- list(FCT.p_value = "FinalPvalue",
                    continuity.correction = 1,
                    tolerance = 1e-3)
    
    if (reinitialise == TRUE) {
        assign(".DelayedGSD-options", 
               default,
               envir = DelayedGSD.env)
    
    return(invisible(get(".DelayedGSD-options", envir = DelayedGSD.env)))
    
    }else{
        args <- list(...)

        if(!is.null(names(args))){
            object <- get(".DelayedGSD-options", envir = DelayedGSD.env)
        }else{
            object <- try(get(".DelayedGSD-options", envir = DelayedGSD.env))
            if(inherits(object,"try-error")){
                object <- default
            }
        }
        
        if(length(args)==0){ ## read all
            return(object)
        }else if (!is.null(names(args))) { ## write

          if(any(names(args) %in% names(object) == FALSE)){
              stop("Incorrect element selected: \"",paste0(names(args)[names(args) %in% names(object) == FALSE], collapse = "\" \""),"\"\n",
                   "Available elements: \"",paste0(setdiff(names(object),names(args)), collapse = "\" \""),"\"\n")
          }

            if("continuity.correction" %in% names(args) && (any(args$continuity.correction %in% 0:2 == FALSE) || length(args$continuity.correction)!=1)){
                stop("Argument \'continuity.correction\' must be a integer with values 0, 1, or 2.\n")
            }
            valid.FCT <- c("FinalPvalue","FinalPvalue2")
            if("FCT.p_value" %in% names(args) && (any(args$FCT.p_value %in% valid.FCT == FALSE) || length(args$FCT.p_value)!=1)){
                stop("Argument \'FCT.p_value\' must be a integer with values \"",paste(valid.FCT, collapse ="\", \""),"\".\n")
            }
            if("tolerance" %in% names(args) && (!is.numeric(args$tolerance) || any(args$tolerance<=0) || length(args$tolerance)!=1)){
                stop("Argument \'tolerance\' must be a strictly positive number.\n")
            }
            object[names(args)] <- args
      
          assign(".DelayedGSD-options", 
                 object, 
                 envir = DelayedGSD.env)
      
          return(invisible(object))
      
      } else {# read
          args <- unlist(args)
          if(any(args %in% names(object) == FALSE)){
              stop("Incorrect element selected: \"",paste0(args[args %in% names(object) == FALSE], collapse = "\" \""),"\"\n",
                   "Available elements: \"",paste0(setdiff(names(object),args), collapse = "\" \""),"\"\n")
          }
          return(object[args])
      }
    
  }
  
  
}


##----------------------------------------------------------------------
### DelayedGSD.options.R ends here
