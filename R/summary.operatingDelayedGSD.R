### summary.operatingDelayedGSD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  2 2024 (15:17) 
## Version: 
## Last-Updated: maj  3 2024 (14:34) 
##           By: Brice Ozenne
##     Update #: 68
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summary.operatingDelayedGSD
##' @title Summary of the Operating Characteristics of a GSD
##' @description Display key simulation parameters (sample size, occurence of a stage, information) and the key operating characteristics of a GSD.
##' @param x output of operatingDelayedGSD.
##' @param print [logical] should the information be printed in the console.
##' @param digits [logical] should the information be printed in the console.
##' @param ... not used, for compatibility with the generic method
##' 
##' @export
summary.operatingDelayedGSD <- function(object, print = TRUE, digits = c(2,4), ...){

    ## ** extract from object
    N.fw <- object$args$N.fw ## number of follow-up values per individual
    n.sim <- object$args$n.sim
    n.obs <- object$args$n.obs
    kMax <- object$args$kMax
    delta <- object$args$delta[N.fw+1]

    test.name <- c(length(unique(object$results$method))>1,length(unique(object$results$binding))>1,length(unique(object$results$fixC))>1)
    col.method <- c("method","binding","fixC")[]
    ls.name <- list(paste0("method ", object$results$method),
                    c(" non-binding"," binding")[object$results$binding+1],
                    c(""," ck>=1.96")[(object$results$method==3)+(object$results$method!=3)*(object$results$fixC+1)])
    if(length(digits)!=2){
        stop("Argument \'digits\' should have length 2. \n",
             "The first element refers to percentages and the second to numbers. \n")
    }

    
    ## ** summary statistics
    if(all(test.name==FALSE)){
        x.results <- cbind(name.method = paste0("method ",object$results$method[1]),object$results)
        level.method <- x.results$name.method[1]
    }else{
        x.results <- cbind(name.method = interaction(as.data.frame(do.call(cbind,ls.name[test.name])), sep ="", drop = TRUE),object$results)
        level.method <- levels(x.results$name.method)
    }
    xd.results <- x.results[x.results$type.stage!="interim",]
    xd.results <- xd.results[order(xd.results$stage),]
    n.method <- length(level.method)

    x.print <- print(object, print = print, digits = digits, ...)

    ## rejection rate
    table.run <- do.call(rbind,by(xd.results, xd.results$name.method,
                                  function(iDF){
                                      return(data.frame("median duration" = median(iDF$time),
                                                        "median sample size" = median(iDF$n.patients),
                                                        check.names = FALSE))
                                  }))

    ## rejection rate
    table.rejection <- do.call(rbind,by(xd.results, xd.results$name.method,
                                        function(iDF){
                                            return(data.frame(n.NNA = sum(!is.na(iDF$decision=="efficacy")),
                                                              "rejection rate" = mean(iDF$decision=="efficacy", na.rm = TRUE),
                                                              check.names = FALSE))
                                            }))

    ## reversal
    vec.reversal <- colSums(do.call(rbind,by(xd.results, xd.results$stage,
                                             function(iDF){
                                                 return(c(tapply(iDF$reason.interim=="futility" & iDF$decision=="efficacy", iDF$name.method,sum, na.rm = TRUE),
                                                          tapply(iDF$reason.interim=="efficacy" & iDF$decision=="futility", iDF$name.method,sum, na.rm = TRUE)))
                                             })))/table.rejection$n.NNA
    table.reversal <- data.frame(n.NNA = table.rejection$n.NNA,
                                 "reversal F->E" = vec.reversal[1:n.method],
                                 "reversal E->F" = vec.reversal[(n.method+1):(2*n.method)],
                                 check.names = FALSE)
    

    ## rejection below 1.96
    table.lowrejection <- do.call(rbind,by(xd.results, xd.results$name.method,
                                           function(iDF){
                                               return(data.frame(n.NNA = sum(!is.na(iDF$decision=="efficacy") & !is.na(iDF$statistic)),
                                                                 "low rejection" = mean(iDF$decision=="efficacy" & iDF$statistic<1.96, na.rm = TRUE),
                                                                 check.names = FALSE))
                                           }))


    ## coverage
    table.coverage <- do.call(rbind,by(xd.results, xd.results$name.method,
                                       function(iDF){
                                           return(data.frame(n.NNA = sum(!is.na(iDF$lower_MUE) & !is.na(iDF$upper_MUE)),
                                                             "coverage" = mean(iDF$lower_MUE <= delta & delta <= iDF$upper_MUE, na.rm = TRUE),
                                                             check.names = FALSE))
                                       }))

    ## sum bias
    table.bias <- do.call(rbind,by(xd.results, xd.results$name.method,
                                       function(iDF){
                                           return(data.frame(n.NNA = sum(!is.na(iDF$estimate_MUE)),
                                                             "mean bias" = mean(iDF$estimate_MUE - delta, na.rm = TRUE),
                                                             "median bias" = mean(iDF$estimate_MUE > delta, na.rm = TRUE) - 0.5,
                                                             check.names = FALSE))
                                       }))


    table.operating <- cbind(rbind(table.run[["median duration"]],
                                   table.run[["median sample size"]],
                                   table.rejection[["rejection rate"]],
                                   table.reversal[["reversal F->E"]],
                                   table.reversal[["reversal E->F"]],
                                   table.lowrejection[["low rejection"]],
                                   table.coverage[["coverage"]],
                                   table.bias[["mean bias"]],
                                   table.bias[["median bias"]]
                                   ),
                             rbind(0,
                                   0,
                                   n.sim - table.rejection[["n.NNA"]],
                                   n.sim - table.reversal[["n.NNA"]],
                                   n.sim - table.reversal[["n.NNA"]],
                                   n.sim - table.lowrejection[["n.NNA"]],
                                   n.sim - table.coverage[["n.NNA"]],
                                   n.sim - table.bias[["n.NNA"]],
                                   n.sim - table.bias[["n.NNA"]]
                                   )
                             )
    colnames(table.operating) <- c(level.method, rep("NA", n.method))
    rownames(table.operating) <- c("median duration","median sample size","rejection rate", "reversal F->E", "reversal E->F", "low rejection", "coverage", "mean bias", "median bias")

    ## necessary because "  " white space will be added to the name latter on
    row.nb <- which(rownames(table.operating) %in% c("median duration","median sample size", "mean bias"))
    row.pc <- which(rownames(table.operating) %in% c("rejection rate", "reversal F->E", "reversal E->F", "low rejection", "coverage", "median bias"))
    row.coverage <- which(rownames(table.operating) == "coverage")
    row.meanbias <- which(rownames(table.operating) == "mean bias")
    row.medianbias <- which(rownames(table.operating) == "median bias")
    
    ## ** display
    if(print){
    cat("\n - operating characteristics (true effect=",delta,"):\n",sep="")
    print.operating <- table.operating[,level.method,drop=FALSE]
    rownames(print.operating) <- paste0("  ",rownames(print.operating))
    print.operating[row.pc,level.method] <- paste0(round(100*unlist(table.operating[row.pc,level.method]), digits = digits[1]),"%")
    print.operating[row.nb,level.method] <- as.character(round(unlist(table.operating[row.nb,level.method]), digits = digits[2]))
    if(any(table.operating[row.coverage,(n.method+1):(2*n.method)]>0)){        
        print.operating[row.coverage,1:n.method] <- paste0(print.operating[row.coverage,1:n.method],
                                                         " (NA: ",round(100*table.operating[row.coverage,(n.method+1):(2*n.method)]/n.sim, digits = digits[1]),"%)")
    }
    if(any(table.operating[row.meanbias,(n.method+1):(2*n.method)]>0)){
        print.operating[row.meanbias,1:n.method] <- paste0(print.operating[row.meanbias,1:n.method],
                                                          " (NA: ",round(100*table.operating[row.meanbias,(n.method+1):(2*n.method)]/n.sim, digits = digits[1]),"%)")
    }
    if(any(table.operating[row.medianbias,(n.method+1):(2*n.method)]>0)){
        print.operating[row.medianbias,1:n.method] <- paste0(print.operating[row.medianbias,1:n.method],
                                                            " (NA: ",round(100*table.operating[row.medianbias,(n.method+1):(2*n.method)]/n.sim, digits = digits[1]),"%)")
    }

    if(n.method == 1){
        print(as.data.frame(as.list(print.operating[,1]), check.names = FALSE), quotes = FALSE, row.names = FALSE)
    }else{
        print(as.data.frame(print.operating, check.rows = FALSE), quotes = FALSE)
    }
    cat("\n")
    }

    ## ** export
    return(invisible(table.operating))
}


##----------------------------------------------------------------------
### summary.operatingDelayedGSD.R ends here
