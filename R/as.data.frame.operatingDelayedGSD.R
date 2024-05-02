### as.data.frame.operatingDelayedGSD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  2 2024 (15:19) 
## Version: 
## Last-Updated: maj  2 2024 (15:20) 
##           By: Brice Ozenne
##     Update #: 2
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * as.data.frame.operatingDelayedGSD
##' @export
as.data.frame.operatingDelayedGSD <- function(x, row.names = NULL, optional = FALSE, ...){
    return(x$results)
}

##----------------------------------------------------------------------
### as.data.frame.operatingDelayedGSD.R ends here
