### print.operatingDelayedGSD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  1 2024 (10:49) 
## Version: 
## Last-Updated: maj  2 2024 (17:00) 
##           By: Brice Ozenne
##     Update #: 50
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * print.operatingDelayedGSD
##' @export
print.operatingDelayedGSD <- function(x, print = TRUE, ...){

    ## ** extract from object
    N.fw <- x$args$N.fw ## number of follow-up values per individual
    n.sim <- x$args$n.sim
    n.obs <- x$args$n.obs
    kMax <- x$args$kMax

    test.name <- c(length(unique(x$results$method))>1,length(unique(x$results$binding))>1,length(unique(x$results$fixC))>1)
    ls.name <- list(paste0("method ", x$results$method),
                        c(" non-binding"," binding")[x$results$binding+1],
                        c(""," ck>=1.96")[(x$results$method==3)+(x$results$method!=3)*(x$results$fixC+1)])

    ## ** summary statistics
    if(all(test.name==FALSE)){
        x.results <- cbind(name.method = paste0("method ",x$results$method[1]),x$results)
    }else{
        x.results <- cbind(name.method = interaction(as.data.frame(do.call(cbind,ls.name[test.name])), sep ="", drop = TRUE),x$results)
    }
    table.n <- do.call(rbind,by(x.results, interaction(x.results$stage,x.results$type.stage),
                                function(iDF){
                                    return(data.frame(
                                        stage = iDF$stage[1],
                                        type = iDF$type.stage[1],
                                        frequency = NROW(iDF),
                                        median.patients = median(iDF$n.patients), min.patients = min(iDF$n.patients), max.patients = max(iDF$n.patients),
                                        median.outcome = median(iDF$n.outcome), min.outcome = min(iDF$n.outcome), max.outcome = max(iDF$n.outcome),
                                        median.pipeline = median(iDF$n.pipeline), min.pipeline = min(iDF$n.pipeline), max.pipeline = max(iDF$n.pipeline),
                                        median.complete = median(iDF$n.complete), min.complete = min(iDF$n.complete), max.complete = max(iDF$n.complete))
                                        )
                                }))
    rownames(table.n) <- NULL
    
    table.n$type <- factor(table.n$type,c("interim","decision","final"))
    table.n <- table.n[order(table.n$stage,table.n$type),]

    
    table.visit <- do.call(rbind,by(x.results, interaction(x.results$stage,x.results$type.stage),
                                    function(iDF){
                                        return(cbind(data.frame(
                                            stage = iDF$stage[1],
                                            type = iDF$type.stage[1]),
                                            as.data.frame(as.list(table(iDF$name.method)), check.names = FALSE)
                                            ))
                                    }))
    table.visit$type <- factor(table.visit$type,c("interim","decision","final"))
    table.visit <- table.visit[order(table.visit$stage,table.visit$type),]

    ## ** prepare txt
    if(any(test.name==FALSE)){
            if(any(test.name[1]==FALSE)){
                txt.method <- paste0(" ",x$results$method[1])
            }else{
                txt.method <-  "s"
            }
            if(any(test.name[2]==FALSE)){
                if(all(x$results$binding==FALSE)){
                    txt.binding <- " with non-binding futility rules"
                }else if(all(x$results$binding==FALSE)){
                    txt.binding <-  " with binding futility rules"
                }
            }else{
                txt.binding <- ""
            }
            if(any(test.name[3]==FALSE)){
                if(all(x$results$fixC==FALSE)){
                    txt.fixC <- " with free decision boundaries"
                }else if(all(x$results$binding==FALSE)){
                    txt.fixC <-  " with decision boundaries above 1.96"
                }
            }else{
                txt.fixC <- ""
            }
            
    }

    ## ** display
    if(print){
        cat("\t\tSimulation study for Group Sequential trial with ",kMax-1," interim. \n",
            "\t\t(",n.sim," simulations with maximum sample size ",n.obs,")\n",sep="")
        if(any(test.name==FALSE)){
            cat("\t\t(method",txt.method,txt.binding,txt.fixC,")\n\n",sep="")
        }else{
            cat("\n")
        }

        cat(" - median [min;max] number of included participant, participant with incomplete follow-up, and participant with observed outcome. \n")
        cat(" - frequency (% w.r.t. number of simulations) of occurence of each stageper method. \n")
        print.n <- data.frame(stage = paste(table.n$type, table.n$stage, sep =" "),
                              included = paste(table.n$median.patients, " [", table.n$min.patients, ";",table.n$max.patients,"]", sep =""),
                              pipeline = paste(table.n$median.pipeline, " [", table.n$min.pipeline, ";",table.n$max.pipeline,"]", sep =""),
                              outcome = paste(table.n$median.outcome, " [", table.n$min.outcome, ";",table.n$max.outcome,"]", sep ="")
                              )

    print.visit <- data.frame(stage = paste(table.visit$type, table.visit$stage, sep =" "),
                              table.visit[, levels(x.results$name.method)], check.names = FALSE)
    print.visit[, levels(x.results$name.method)] <- paste0(unlist(print.visit[, levels(x.results$name.method)]), " (", round(100*unlist(print.visit[, levels(x.results$name.method)])/n.sim,2),"%)")
    print(cbind(print.n,print.visit[-1]), row.names = FALSE)
    }

    ## ** export
    return(invisible(cbind(table.n,table.visit[-(1:2)])))
}

##----------------------------------------------------------------------
### print.operatingDelayedGSD.R ends here
