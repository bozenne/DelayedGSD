### print.operatingDelayedGSD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  1 2024 (10:49) 
## Version: 
## Last-Updated: maj  3 2024 (11:34) 
##           By: Brice Ozenne
##     Update #: 74
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
print.operatingDelayedGSD <- function(x, index.method = 1, print = TRUE, digits = c(2,2), ...){

    ## ** extract from object
    N.fw <- x$args$N.fw ## number of follow-up values per individual
    n.sim <- x$args$n.sim
    n.obs <- x$args$n.obs
    kMax <- x$args$kMax

    test.name <- c(length(unique(x$results$method))>1,length(unique(x$results$binding))>1,length(unique(x$results$fixC))>1)
    ls.name <- list(paste0("method ", x$results$method),
                        c(" non-binding"," binding")[x$results$binding+1],
                        c(""," ck>=1.96")[(x$results$method==3)+(x$results$method!=3)*(x$results$fixC+1)])

    if(length(index.method)!=1){
        stop("Argument \'index.method\' should have length 1. \n")
    }
    if(length(digits)!=2){
        stop("Argument \'digits\' should have length 2. \n",
             "The first element refers to percentages and the second to numbers. \n")
    }
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** summary statistics
    if(all(test.name==FALSE)){
        x.results <- cbind(name.method = paste0("method ",x$results$method[1]),x$results)
        level.method <- x.results$name.method[1]
    }else{
        x.results <- cbind(name.method = interaction(as.data.frame(do.call(cbind,ls.name[test.name])), sep ="", drop = TRUE),x$results)
        level.method <- levels(x.results$name.method)
    }
    if(index.method <= 0 || index.method>length(level.method)){
        stop("Argument \'index.method\' should take integer values between 1 and ",length(level.method),". \n",sep="")
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
    
    table.info <- do.call(rbind,by(x.results, interaction(x.results$stage,x.results$type.stage),
                                   function(iDF){
                                       return(cbind(data.frame(
                                           stage = iDF$stage[1],
                                           type = iDF$type.stage[1]),
                                           as.data.frame(as.list(tapply(iDF$info,iDF$name.method,mean,na.rm=TRUE)), check.names = FALSE)
                                           ))
                                   }))
    table.info$type <- factor(table.info$type,c("interim","decision","final"))
    table.info <- table.info[order(table.info$stage,table.info$type),]

    table.infoPC <- do.call(rbind,by(x.results, interaction(x.results$stage,x.results$type.stage),
                                   function(iDF){
                                       return(cbind(data.frame(
                                           stage = iDF$stage[1],
                                           type = iDF$type.stage[1]),
                                           as.data.frame(as.list(tapply(iDF$infoPC,iDF$name.method,mean,na.rm=TRUE)), check.names = FALSE)
                                           ))
                                   }))
    table.infoPC$type <- factor(table.infoPC$type,c("interim","decision","final"))
    table.infoPC <- table.infoPC[order(table.infoPC$stage,table.infoPC$type),]

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
        cat(" - frequency (% w.r.t. number of simulations) for ",level.method[index.method],". \n", sep = "")
        cat(" - information (% w.r.t. to planned max information) for ",level.method[index.method],". \n", sep = "")
        print.n <- data.frame(stage = paste(table.n$type, table.n$stage, sep =" "),
                              included = paste(table.n$median.patients, " [", table.n$min.patients, ";",table.n$max.patients,"]", sep =""),
                              pipeline = paste(table.n$median.pipeline, " [", table.n$min.pipeline, ";",table.n$max.pipeline,"]", sep =""),
                              outcome = paste(table.n$median.outcome, " [", table.n$min.outcome, ";",table.n$max.outcome,"]", sep ="")
                              )

        print.visit <- data.frame(stage = paste(table.visit$type, table.visit$stage, sep =" "),
                                  table.visit[, levels(x.results$name.method)], check.names = FALSE)
        print.visit[, levels(x.results$name.method)] <- paste0(unlist(print.visit[, levels(x.results$name.method)]), " (", round(100*unlist(print.visit[, levels(x.results$name.method)])/n.sim,digits[1]),"%)")
        
        print.info <- data.frame(stage = paste(table.info$type, table.info$stage, sep =" "),
                                 table.info[, levels(x.results$name.method)], check.names = FALSE)
        print.info[, levels(x.results$name.method)] <- paste0(round(unlist(print.info[, levels(x.results$name.method)]),digits[2]), " (", round(100*unlist(table.infoPC[, levels(x.results$name.method)]),digits[1]),"%)")

        print(cbind(print.n, frequency = print.visit[[1+index.method]], information = print.info[[1+index.method]]), row.names = FALSE)
    }

    ## ** export
    return(invisible(list(n = table.n,
                          visit = table.visit,
                          info = table.info,
                          infoPC = table.infoPC)))
}

##----------------------------------------------------------------------
### print.operatingDelayedGSD.R ends here
