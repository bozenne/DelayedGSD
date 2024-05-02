### print.operatingDelayedGSD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  1 2024 (10:49) 
## Version: 
## Last-Updated: maj  1 2024 (11:38) 
##           By: Brice Ozenne
##     Update #: 12
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

print.operatingDelayedGSD <- function(x, ...){

    ## ** extract from object
    N.fw <- x$args$N.fw ## number of follow-up values per individual
    n.sim <- x$args$n.sim
    n.obs <- x$args$n.obs
    kMax <- x$args$kMax

    ## ** summary statistics
    results.decision <- x$results[x$results$type.stage %in% c("decision","final"),] 
    by(results.decision, interaction(results.decision$method,results.decision$binding,results.decision$fixC),
       function(iDF){
           return(data.frame(n = NROW(iDF) rejection.rate = mean(iDF$decision=="efficacy", na.rm = TRUE)))
       })
    browser()

    ## ** display
    cat("\t\tSimulation study for Group Sequential trial with ",kMax-1," interim. \n",
        "\t\t(",n.sim," iterations with maximum sample size ",n.obs,")\n\n",sep="")

    
    return(invisible(NULL))
}

##----------------------------------------------------------------------
### print.operatingDelayedGSD.R ends here
