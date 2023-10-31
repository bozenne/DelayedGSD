## * summary.delayedGSD (documentation)
#' @title Status of a Group Sequential Design with Delayed Endpoints
#' @description Display boundaries and estimated treatment effect up the current stage.
#' 
#' @param x object of class \code{delayedGSD}, typically output from \code{\link{CalcBoundaries}}.
#' @param planned [logical] should the planned boundaries be displayed along with the boundaries that have been computed based on the data.
#' Can also be \code{"only"} to only display planned boundaries.
#' @param predicted [logical] Should the predicted information/boundaries at decision based on interim data be output (when relevant).
#' @param digits [interger,>0] integer indicating the number of decimal places to be used.
#' @param space [character] column separator.
#' @param abreviated [logical] abreviate some column names and action that should be taken
#' (e.g. C for continue, C-Info for continue due to decreasing information, S-F for stop for futility, and S-Imax for stop because maximum information has been reached).
#' @param ... not used, for compatibility with the generic method.
#' 

## * summary.delayedGSD
#' @export
summary.delayedGSD <- function(x, planned = NULL, predicted = TRUE, digits = 5, space = "|", abreviated = TRUE, ...){

    ## ** extract information from object
    call <- x$call
    kMax <- x$kMax
    k <- x$stage$k
    type.k <- x$stage$type
    
    if(is.null(planned)){
        if(type.k=="planning"){
            planned <- "only"
        }else{
            planned <- FALSE
        }
        ## if(type.k == "interim" && x$conclusion["interim",k]!="stop"){
        ##     planned <- TRUE            
        ## }else{
        ##     planned <- FALSE
        ## }
    }else if(planned %in% c(TRUE,FALSE) == FALSE && planned != "only"){
        stop("Argument \'planned\' should be TRUE, FALSE, or \"only\". \n")
    }

    ## ** Welcome message
    if(planned=="only"){

        cat("	    Planning for a GSD with repeated measurements \n")
        if(abreviated<=1){
            cat("\n * Call \n")
            print(call)
        }

    }else{

        if(type.k == "final"){
            cat("	    GSD with repeated measurements at the final analysis (stage ",k,") \n", sep = "")
        }else if(type.k == "decision"){
            cat("	    GSD with repeated measurements at the decision analysis of stage ",k," \n", sep = "")
        }else if(type.k == "interim"){
            cat("	    GSD with repeated measurements at the interim analysis of stage ",k," \n", sep = "")
        }

    }

    ## ** Extract information
    if(planned=="only"){
        outInfo <- coef(x, type = "information", planned = TRUE)
        outBound <- coef(x, type = "boundary", planned = TRUE)
        outDelta <- coef(x, type = "effect", planned = TRUE)
    }else{
        outInfo <- coef(x, type = "information", planned = FALSE, predicted = predicted)
        outBound <- coef(x, type = "boundary", planned = FALSE, predicted = predicted)
        outDelta <- stats::confint(x)
        if(planned && type.k == "interim"){ ## add planned values
            outInfo2 <- coef(x, type = "information", planned = TRUE)
            outInfo$Interim[(k+1):kMax] <- outInfo2$Interim[(k+1):kMax]
            outInfo$Interim.pc[(k+1):kMax] <- outInfo2$Interim.pc[(k+1):kMax]
            outInfo$Decision[(k+1):kMax] <- outInfo2$Decision[(k+1):kMax]
            outInfo$Decision.pc[(k+1):kMax] <- outInfo2$Decision.pc[(k+1):kMax]
            outBound2 <- coef(x, type = "boundary", planned = TRUE)
            outBound$Fbound[(k+1):kMax] <- outBound2$Fbound[(k+1):kMax]
            outBound$Ebound[(k+1):kMax] <- outBound2$Ebound[(k+1):kMax]
            outBound$Cbound[(k+1):kMax] <- outBound2$Cbound[(k+1):kMax]
        }
    }
    outInfo[NROW(outInfo),c("Interim","Interim.pc")] <- NA
    outBound[NROW(outBound),c("Fbound","Ebound")] <- NA

    ## ** Boundaries
    df.printBound <- outBound
    if(!is.null(digits)){        
        df.printBound$Fbound <- round(df.printBound$Fbound, digits)
        df.printBound$Ebound <- round(df.printBound$Ebound, digits)
        df.printBound$Cbound <- round(df.printBound$Cbound, digits)
    }


    if(planned=="only"){

        if(abreviated>0){
            colnames(df.printBound) <- c("stage","F-bound","E-bound","C-bound")
        }else{
            colnames(df.printBound) <- c("stage","Futility boundary","Efficacy boundary","Critical boundary")
        }
        df.printBound[["alpha-spent"]] <- round(x$planned$alphaSpent, digits)
        df.printBound[["beta-spent"]] <- round(x$planned$betaSpent, digits)
        df.printBound[is.na(df.printBound)] <- ""

        if(!is.null(space)){
            df.printBound <- cbind(df.printBound[,1,drop=FALSE],
                                    " " = space,
                                    df.printBound[,2:3,drop=FALSE],
                                    " " = space,
                                    df.printBound[,4,drop=FALSE],
                                    " " = space,
                                    df.printBound[,5:6,drop=FALSE]
                                   )
            names(df.printBound)[c(2,5,7)] <- space
        }
        cat("\n * Planned boundaries: \n")
        print(df.printBound, row.names = FALSE, quotes="")

    }else{
        df.printBound$statistic.interim <- round(df.printBound$statistic.interim, digits)
        df.printBound$statistic.decision <- round(df.printBound$statistic.decision, digits)
        df.printBound[["alpha-spent"]] <- round(x$alphaSpent, digits+2)
        df.printBound[["beta-spent"]] <- round(x$betaSpent, digits+2)
        if(planned){
            df.printBound[["alpha-spent"]][(k+1):kMax] <- x$planned$alphaSpent[(k+1):kMax]
            df.printBound[["beta-spent"]][(k+1):kMax] <- x$planned$betaSpent[(k+1):kMax]
        }
        df.printBound[is.na(df.printBound)] <- ""
        conclusion.interim.print <- .reformatInterimTxt(x$conclusion, kMax = kMax, abreviated = abreviated>0)
        conclusion.decision.print <- .reformatDecisionTxt(x$conclusion, kMax = kMax, abreviated = abreviated>0)
        if(any(conclusion.interim.print=="S-Imax")){ ## when stopping at interim because we anticipate to reach Imax at decision put NA to boundaries
            df.printBound$Fbound[conclusion.interim.print=="S-Imax"] <- as.character(NA)
            df.printBound$Ebound[conclusion.interim.print=="S-Imax"] <- as.character(NA)
            df.printBound$statistic.interim[conclusion.interim.print=="S-Imax"] <- as.character(NA)
        }
        df.printBound2 <- data.frame("stage" = df.printBound$stage,
                                     "Futility boundary" = df.printBound$Fbound,
                                     "Efficacy boundary" = df.printBound$Ebound,
                                     "Statistic" = df.printBound$statistic.interim,
                                     "Action" = conclusion.interim.print,
                                     "Critical boundary" = df.printBound$Cbound,
                                     "Statistic" = df.printBound$statistic.decision,
                                     "Action" = conclusion.decision.print,
                                     "alpha" = df.printBound[["alpha-spent"]],
                                     "beta" = df.printBound[["beta-spent"]],
                                     check.names = FALSE)
        if(abreviated>0){
            names(df.printBound2) <- c("stage", "F-bound", "E-bound", "Stat", "", "C-bound", "Stat", "","alpha","beta")
        }
        df.printBound3 <- rbind(names(df.printBound2),df.printBound2)
        df.printBound3[1,1] <- ""
        
        if(!is.null(space)){
            df.printBound3 <- cbind(df.printBound3[,1,drop=FALSE],
                                    " " = space,
                                    df.printBound3[,2:5,drop=FALSE],
                                    " " = space,
                                    df.printBound3[,6:8,drop=FALSE],
                                    " " = space,
                                    df.printBound3[,9:10,drop=FALSE]
                                    )
            names(df.printBound3)[-1] <- ""
            names(df.printBound3)[c(2,7,11)] <- space
            index.name <- c(4,8,12)
        }else{
            names(df.printBound3)[-1] <- ""
            index.name <- c(3,6,9)
        }
        if(abreviated>0){
            names(df.printBound3)[index.name[1]] <- "Interim"
            names(df.printBound3)[index.name[2]] <- "Decision"
            names(df.printBound3)[index.name[3]] <- "Spent"
            if(abreviated>1){
                df.printBound3[,index.name[3]:NCOL(df.printBound3)] <- NULL
            }
        }else{
            names(df.printBound3)[index.name[1]] <- "Interim analysis"
            names(df.printBound3)[index.name[2]] <- "Decision analysis"
            names(df.printBound3)[index.name[3]] <- "Error spent"
        }
        
        ## display
        cat("\n * Boundaries and observed statistics \n")
        print(df.printBound3, row.names = FALSE, quote = FALSE)
    }
    cat("\n")

    ## ** Information
    df.printInfo <- outInfo
    if(!is.null(digits)){        
        df.printInfo$Interim <- round(df.printInfo$Interim, digits)
        df.printInfo$Decision <- round(df.printInfo$Decision, digits)
        df.printInfo$Interim.pc <- round(df.printInfo$Interim.pc, digits)
        df.printInfo$Decision.pc <- round(df.printInfo$Decision.pc, digits)
    }

    if(planned == "only"){
        df.printInfo$n <- NULL
        names(df.printInfo) <- c("stage", "Interim", "(%)", "Decision", "(%)")
        cat(" * Planned information: \n",sep="")
    }else{
        names(df.printInfo) <- c("stage", "Interim", "(%)", "Decision", "(%)", "n")
        if(planned){
            cat(" * Observed and planned information: \n",sep="")
        }else{
            cat(" * Observed and predicted information: \n",sep="")
        }
    }
    df.printInfo[is.na(df.printInfo)] <- ""
    if(!is.null(space)){
        if(planned == "only" || all(df.printInfo$n=="")){
            df.printInfo <- cbind(df.printInfo[,1,drop=FALSE]," "=space,df.printInfo[,2:3]," "=space,df.printInfo[,4:5]," "=space)
        }else{
            df.printInfo <- cbind(df.printInfo[,1,drop=FALSE]," "=space,df.printInfo[,2:3]," "=space,df.printInfo[,4:5]," "=space,df.printInfo[,6,drop=FALSE])
        }
        names(df.printInfo)[c(2,5,8)] <- space
    }
    print(df.printInfo, row.names = FALSE)
   
    ## ** Estimates
    if(planned=="only" && abreviated<=1){
        cat("\n")
        cat(" * Expected effect : ",outDelta," \n",sep="")
    }else if(abreviated<=1){
        x.CI <- stats::confint(x, method = c("ML","MUE"))
        x.CI$estimate <- round(x.CI$estimate, digits)
        x.CI$se <- round(x.CI$se, digits)
        x.CI$lower <- round(x.CI$lower, digits)
        x.CI$upper <- round(x.CI$upper, digits)
        x.CI$statistic <- round(x.CI$statistic, digits)
        x.CI$df <- round(x.CI$df, digits)
        x.CI$p.value <- format.pval(x.CI$p.value, digits)
        x.CI$stage <- paste0(x.CI$stage," (",x.CI$type,")")
        x.CI$type <- NULL

        if(all(x.CI$method == "ML")){
            cat("\n")
            if(!is.na(x.CI$coef[1])){
                cat(" * Current ML-estimate of the treatment effect (",x.CI$coef[1],") \n",sep="")
            }else{
                cat(" * Current ML-estimate of the treatment effect \n",sep="")
            }
            print(x.CI[x.CI$method == "ML",c("estimate","se","lower","upper","statistic","df","p.value")], row.names = FALSE)
        }else if("MUE" %in% x.CI$method){
            cat("\n")
            if(!is.na(x.CI$coef[1])){
                cat(" * Current MUE-estimate of the treatment effect (",x.CI$coef[1],") \n",sep="")
            }else{
                cat(" * Current MUE-estimate of the treatment effect \n",sep="")
            }
            print(x.CI[x.CI$method == "MUE",c("estimate","lower","upper","p.value")], row.names = FALSE)
        }else{
            cat("\n")
            if(!is.na(x.CI$coef[1])){
                cat(" * Current estimates of the treatment effect (",x.CI$coef[1],") \n",sep="")
            }else{
                cat(" * Current estimates of the treatment effect \n",sep="")
            }
            x.CI$se <- as.character(round(x.CI$se, digits))
            x.CI$statistic <- as.character(round(x.CI$statistic, digits))
            x.CI$df <- as.character(round(x.CI$df))
            print(x.CI[,c("method","estimate","se","lower","upper","statistic","df","p.value")], row.names = FALSE, na.print = "", digits = digits)
        }
    }
  
    ## ** Sample size
    if(planned=="only" && abreviated<=1){
        cat("\n")
        cat(" * Inflation factor: ",x$planned$InflationFactor," \n",sep="")
        if(!is.na(x$planned$n.obs)){
            cat(" * Expected sample size: ",x$planned$n.obs," \n",sep="")
        }
        cat("\n")
    }else if(abreviated<=1 && !is.na(x$n.obs[k])){

        iLMM <- x$lmm[[k+(type.k=="decision")]]

        df.nobs <- data.frame(c1 = paste0(" * ",x$n.obs[k]," patients in the study"))
        if(!is.null(iLMM)){
            df.nobs$c1 <- paste0(df.nobs$c1,":")
            df.nobs$c2 <- paste0(iLMM$sample.size["interim"]," with at least one outcome value ")
            df.nobs <- rbind(df.nobs, data.frame(c1 = "", c2 = paste0(iLMM$sample.size["interim.cc"]," with complete data")))
        }
        colnames(df.nobs) <- NULL
        print(df.nobs, quote = FALSE, right = FALSE, row.names = FALSE)
        cat("\n")
    }else{
        cat("\n")
    }
    return(invisible(NULL))
}


## * reformat text
## ** interim
.reformatInterimTxt <- function(conclusion, kMax, abreviated){

    if(all(is.na(conclusion["interim",]))){return(rep("",kMax))}
    
    out <- sapply(1:kMax,function(iK){
        if(is.na(conclusion["interim",iK])){
            return("")
        }else if(conclusion["interim",iK]=="continue"){
            if(conclusion["reason.interim",iK]=="no boundary crossed"){
                if(abreviated){
                    return("C")
                }else{
                    return("continue")
                }
            }else if(conclusion["reason.interim",iK]=="decreasing information"){
                if(abreviated){
                    return("C-Info")
                }else{
                    return("continue (decreasing information)")
                }
            }else if(conclusion["reason.interim",iK]=="overrule futility"){
                if(abreviated){
                    return("C-OF")
                }else{
                    return("continue (overrule futility)")
                }
            }
        }else if(conclusion["interim",iK]=="stop"){
            if(conclusion["reason.interim",iK]=="efficacy"){
                if(abreviated){
                    return("S-E")
                }else{
                    return("stop for efficacy")
                }
            }else if(conclusion["reason.interim",iK]=="futility"){
                if(abreviated){
                    return("S-F")
                }else{
                    return("stop for futility")
                }
            }else if(conclusion["reason.interim",iK]=="Imax reached"){
                if(abreviated){
                    return("S-Imax")
                }else{
                    return("stop (maximum information reached)")
                }
            }
        }
    })
    return(out)
}

## ** decision
.reformatDecisionTxt <- function(conclusion, kMax, abreviated){

    out <- rep("",kMax)

    if(all(is.na(conclusion["decision",]))){return(out)}

    out <- sapply(1:kMax,function(iK){
        if(is.na(conclusion["decision",iK])){
            return("")
        }else if(conclusion["decision",iK]=="efficacy"){
            if(abreviated){
                return("E")
            }else{
                return("conclude efficacy")
            }
            return(out)
        }else if(conclusion["decision",iK]=="futility"){
            if(abreviated){
                return("F")
            }else{
                return("conclude futility")
            }
        }
    })

    return(out)
}
