## * analyzeData (documentation)
#' @title Fit a Linear Mixed Model
#' @description Fit a linear mixed model and extract the information relative to the parameter of interest.
#' 
#' @param data [data.frame] dataset for the current analysis.
#' @param ddf [character] method used to compute the degrees of freedom of the Wald statistic.
#' Can be \code{"satterthwaite"} or \code{"nlme"}.
#' @param getinfo [logical] should the information be computed at interim and decision? Otherwise no information is computed.
#' @param data.decision [integer or data.frame] data or total number of patients for the future decision analysis (only relevant at interim).
#' Used to compute the predicted total information at decision where only current (or future) drop-out value are excluded but future observations are kept.
#' @param trace [logical] should the execution of the function be traced?
#' 
#' @examples
#' ## simulate data in the wide format
#' set.seed(1)
#' dataW <- GenData(n=104)$d
#'
#' ## estimate mixed model
#' analyzeData(dataW, getinfo = FALSE)
#' analyzeData(dataW, getinfo = TRUE)
#' analyzeData(dataW, ddf = "satterthwaite", getinfo = TRUE)
#'


## * analyzeData (code)
#' @export
analyzeData <- function(data, ddf = "nlme", getinfo = TRUE, data.decision = NULL, trace = TRUE){

    requireNamespace("nlme")

    ## ** normalize arguments
    ddf <- match.arg(ddf, choices = c("nlme","satterthwaite"))
    long <- wide2long(data, rm.na = TRUE) ## remove missing value but not pipeline value from the data.set
    if(!is.null(data.decision) && inherits(data.decision, what = "data.frame")){
        data.decision <- wide2long(data.decision, rm.na = TRUE) ## remove missing value
    }
    ## head(long[order(long$id),],n=10)
    ## summary(long)

    ## ** fit gls model
    ctrl <- nlme::glsControl(opt='optim')

    m <- nlme::gls(X ~ baseline*visit + Z*visit,
                   data = long,
                   correlation = nlme::corSymm(form=~visit.num|id), 
                   weights = nlme::varIdent(form=~1|visit),
                   method = "REML",
                   na.action = stats::na.exclude)
    m$variable <- c(cluster = "id", time = "visit")

    ## ** store estimate from the gls model
    out <- list(data=long,
                fit=m,
                name.coef = "Z1",
                delta = data.frame(estimate = as.double(stats::coef(m)["Z1"]),
                                   se = as.double(sqrt(m$varBeta["Z1","Z1"])),
                                   statistic = NA,
                                   df = NA,
                                   p.value = NA)
                )
    #out$delta$statistic <- abs(out$delta$estimate)/out$delta$se
    out$delta$statistic <- out$delta$estimate/out$delta$se
    
    ## computation of the p-value  
    if(ddf=="nlme"){ ## with weird estimator of the degree of freedom from nlme::gls 
        out$delta$df <- as.double(m$dims$N - m$dims$p)
        out$delta$p.value <- as.double(1-stats::pt(out$delta$statistic, df = out$delta$df))  
        #out$delta$p.value <- as.double(2*(1-stats::pt(out$delta$statistic, df = out$delta$df)))#TO DO: this looks like a two-sided p-value, shouldn't it be one-sided?
        ## summary(m)$tTable["Z1","p-value"]
    }else{ ## or using satterthwaite approximation
        requireNamespace("emmeans")
        
        groupTest <- suppressMessages(emmeans::emmeans(m, specs = ~Z|visit, data = long[!is.na(long$X),,drop=FALSE]))
        e.satterthwaite <- summary(graphics::pairs(groupTest, reverse = TRUE), by = NULL, infer = TRUE, adjust = "none", side = ">")
        index <- which(e.satterthwaite$visit==levels(long$visit)[1])

        if(abs(out$delta$estimate - e.satterthwaite$estimate[index])>1e-10){
            warning("Something went wrong when extracting the coefficient from emmeans \n",
                    "discrepancy estimate emmeans vs. gls: ",out$delta$estimate-e.satterthwaite$estimate[index])
        }
        if(abs(out$delta$se - e.satterthwaite$SE[index])>1e-10){
            warning("Something went wrong when extracting the standard error in emeans \n",
                    "discrepancy se emmeans vs. gls: ",out$delta$se-e.satterthwaite$SE[index])
        }
        out$delta$df <- e.satterthwaite$df[index]
        out$delta$p.value <- e.satterthwaite$p.value[index]
    }

    ## ** Estimate the information
    if(getinfo){
        out$data.decision <- data.decision
        if(is.null(data.decision) || inherits(data.decision, what = "data.frame")){
            ## current information and information at decision
            out <- c(out,getInformation(m, name.coef = "Z1", data = long, newdata = data.decision, details = TRUE))
            if(!is.null(data.decision)){
                out$sample.size["decision"] <- length(unique(data.decision[["id"]]))
            }
            ## details: extract information at interim (complete case) interim(full information)
            ## to get information at decision (i.e. including rows for which the covariates are missing) upweight the rows with observed values.
        }else if(inherits(data.decision,"numeric") || inherits(data.decision,"integer")){
            out <- c(out,getInformation(m, name.coef = "Z1", data = long, newdata = NULL, details = TRUE))
            if(data.decision < out$sample.size["total"]){
                stop("Argument \'data.decision\' should be larger than the number of patients (=",out$sample.size["total"],") in the dataset. \n")
            }
            out$information["decision"] <- (data.decision/out$sample.size["total"])*as.double(out$information["decision"])
            out$sample.size["decision"] <- data.decision
        }else{
            stop("When not NULL, argument \'data.decision\' should be a data.frame or an integer. \ns")
        }
        
    
        if(any(abs(out$Info$interim$vcov - stats::vcov(m))>1e-10)){
            warning("Something went wrong when extracting the variance covariance matrix from gls \n",
                    "largest discrepancy between getInformation and gls: ",max(abs(out$Info$interim$vcov - stats::vcov(m))))
     
        }
    }

    ## ** Export
    class(out) <- "lmmGSD"
    return(out)
}


## * wide2long
wide2long <- function(d, rm.na = FALSE, id.na = NULL, Z.id.na = TRUE){

    dW <- as.data.frame(d)

    ## get the number of timepoints
    col.time <- grep("^t([[:digit:]]+)$",names(dW), value = TRUE)
    col.X <- grep("^X([[:digit:]]+)$*",names(dW), value = TRUE)
    col.missing <- grep("^missing([[:digit:]]+)$*",names(dW), value = TRUE)
    if(length(col.missing)==0){
        for(iCol in col.X){ ## iCol <- col.X[1]
            dW[[paste0("missing",gsub("^X","",iCol))]] <- as.numeric(is.na(dW[[iCol]]))
        }
        col.missing <- grep("^missing([[:digit:]]+)$*",names(dW), value = TRUE)
    }
    n.times <- length(grep("t",names(dW)))

    ## convert to long format
    dL <- stats::reshape(dW,
                         idvar = c("id","Z","missing1","t1","X1"),
                         timevar = "visit",
                         v.names = c("missing","time","X"),
                         varying = list(c(col.missing[-1]),
                                        c(col.time[-1]),
                                        c(col.X[-1])),
                         direction = "long")
    rownames(dL) <- NULL

    names(dL)[names(dL)=="X1"] <- "baseline"
    dL$id <- factor(dL$id)
    dL$Z <- factor(dL$Z)
    dL$visit.num <- as.numeric(dL$visit)
    dL$visit <- stats::relevel(as.factor(dL$visit),ref=n.times-1)## relevel to make the treatment coefficient have the interpretation of threatment effect at end of follow-up
    ## add missing values for not yet observed ids
    if(!is.null(id.na)){
        index.notObs <- which(dL$id %in% id.na)
        dL[index.notObs, c("missing1","missing")] <- -1
        dL[index.notObs, c("t1","baseline","time","X")] <- NA
        if(Z.id.na){
            dL[index.notObs, c("Z")] <- NA
        }
    }
    dL.save <- dL

    ## remove lines corresponding to missing data
    if(rm.na){
        dL <- dL[(dL$missing1<=0),] ## at baseline 
        dL <- dL[(dL$missing<=0),] ## or at follow-up
        dL$id <- droplevels(dL$id)
        dL$Z <- droplevels(dL$Z)
        dL$visit <- droplevels(dL$visit)
    }

    ## export
    attr(dL,"df.allobs") <- dL.save
    return(dL)
}
