### getInformation.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 11 2020 (10:18) 
## Version: 
## Last-Updated: feb 18 2022 (10:50) 
##           By: Brice Ozenne
##     Update #: 1098
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * getInformation (documentation)
#' @title Extract information relative to a parameter.
#' @description Extract information relative to a parameter.
#' @name getInformation
#'
#' @param object a \code{ttest} object or a \code{gls} object.
#' @param name.coef [character] For which coefficient the information should be computed (only relevant for gls models).
#' @param type [character] Should the information be estimated only for the observed data (\code{"estimation"}, excluding missing values)
#' or for the full dataset  (\code{"prediction"}, including missing values).
#' @param variance [NULL, list] If \code{NULL} then the estimated variance-covariance matrix from the object is used.
#' Otherwise it should be a list containing the theoretical value of the variance-covariance matrix for each treatment group.
#' @param newdata,data [data.frame] The dataset relative to which the information should be computed. See details section.
#' @param details [logical] Should intermediate results be output. See details section.
#' @param weighting [logical] Should weight be used to handle missing values in the covariates.
#' @param ... not used. For compatibility with the generic details.
#' 
#' @details Argument \bold{data}: the dataset may contain missing value in the outcome but no in the covariates. Missing value in the outcome indicates that the information is not available at the interim analysis but will be come available at decision.
#'
#'
#' Argument \bold{details}: when using gls models, an attribute detail is added to the output which contain a list:
#' \itemize{
#' \item decision: information at decision analysis
#' \item interim: information at the interim analysis using all available data
#' \item interim.cc: information at the interim analysis using a complete case analysis
#' \item n: sample size at decision, interim with complete observation, interim with only partial observations
#' }
#' 
#' @export
`getInformation` <- function(object, ...) UseMethod("getInformation")

## * getInformation (examples)
#' @examples
#' library(nlme)
#' n <- 1e2
#'
#' #### Single endpoint ####
#' ## simulate data
#' set.seed(10)
#' X <- rnorm(n)
#' Y <- rnorm(n)
#' df <- rbind(data.frame(group=0,value=X),
#'             data.frame(group=1,value=Y))
#'
#' ## t-test
#' getInformation(ttest(value~1, df[df$group==0,])) ## only work for R>=4.0
#' getInformation(ttest(x = X)) ##  based on the estimate variance
#' getInformation(ttest(X), variance = 1) ## based on a variance of 1
#' getInformation(ttest(X), variance = 2) ## based on a variance of 2
#' getInformation(ttest(rnorm(length(X))), variance = 2) ## note: the X values do not matter here
#' 
#' getInformation(ttest(value~group, data = df))
#' getInformation(ttest(x = X, Y))
#' getInformation(ttest(X,Y), variance = 1:2) ## information with a variance of 1 in one group and 2 in the other group
#' getInformation(ttest(X,Y), variance = c(1,3)) ## information with a variance of 1 in one group and 3 in the other group
#'
#' ## gls
#' library(nlme)
#'
#' e1.gls <- gls(value~1, data = df[df$group==0,])
#' getInformation(e1.gls, name.coef = "(Intercept)")
#' getInformation(e1.gls, name.coef = "(Intercept)", variance = matrix(1,1,1))
#' 
#' e2.gls <- gls(value~group, data = df, weights = varIdent(form=~1|group))
#' getInformation(e2.gls, name.coef = "group")
#' getInformation(e2.gls, name.coef = "group", variance = list(matrix(1,1,1),matrix(2,1,1)))
#' 
#' #### Two endpoints ####
#' ## simulate data
#' library(mvtnorm)
#' 
#' set.seed(10)
#' X <- rmvnorm(n, sigma = 0.5 + diag(0.5,2))
#' Y <- rmvnorm(n, sigma = 0.5 + diag(0.5,2))
#' df <- rbind(data.frame(id = paste0("id",1:n),group=0,time=0,value=X[,1]),
#'             data.frame(id = paste0("id",1:n),group=0,time=1,value=X[,2]),
#'             data.frame(id = paste0("id",n+1:n),group=1,time=0,value=Y[,1]),
#'             data.frame(id = paste0("id",n+1:n),group=1,time=1,value=Y[,2]))
#'
#'
#' ## gls
#' e.gls <- gls(value~time-1, data = df[df$group==0,],
#'              correlation = corSymm(form=~1|id),
#'              weights = varIdent(form=~1|time))
#' getInformation(e.gls, name.coef = "time")  ## 1/vcov(e.gls)
#' getInformation(e.gls, name.coef = "time", variance = list(diag(1:2)))
#'
#' ## with interaction
#' e.gls <- gls(value~time*group, data = df,
#'              correlation = corSymm(form=~1|id),
#'              weights = varIdent(form=~1|time*group))
#' getInformation(e.gls, name.coef = "time:group") ## 1/vcov(e.gls)
#' getInformation(e.gls, name.coef = "time", variance = list(diag(1:2),diag(1:2)))
#'
#' ## with random ordering
#' df2 <- df[sample.int(NROW(df)),,drop=FALSE]
#' e2.gls <- gls(value~time*group, data = df2,
#'              correlation = corSymm(form=~1|id),
#'              weights = varIdent(form=~1|time*group))
#' getInformation(e2.gls, name.coef = "time:group") ## 1/vcov(e2.gls)
#' getInformation(e.gls, name.coef = "time", variance = list(diag(1:2),diag(1:2)))

## * getInformation.matrix
#' @rdname getInformation
#' @export
getInformation.matrix <- function(object, variance, ...){
    ## object is the design matrix

    ## ** get hidden arguments
    dots <- list(...)
    index.variance <- dots$index.variance
    index.cluster <- dots$index.cluster
    weight <- dots$weight
    
    ## ** prepare
    precision <- lapply(variance,solve)
    n.cluster <- length(index.variance)
    n.allcoef <- NCOL(object)
    name.allcoef <- colnames(object)
    Info <- matrix(0, nrow = n.allcoef, ncol = n.allcoef,
                   dimnames = list(name.allcoef, name.allcoef))

    ## ** compute information
    for(iId in 1:n.cluster){ ## iId <- 200
        iX <- object[index.cluster==iId,,drop=FALSE]
        iOmega <- precision[[index.variance[iId]]]
        if(!is.null(weight)){
            Info <- Info + weight[iId] * (t(iX) %*% iOmega %*% iX)
        }else{
            Info <- Info + (t(iX) %*% iOmega %*% iX)
        }
    }

    ## ** export
    return(Info)
}

## * getInformation.ttest
#' @rdname getInformation
#' @export
getInformation.ttest <- function(object, type = "estimation", variance = NULL, ...){

    ## ** normalize arguments
    type <- match.arg(type, c("estimation","prediction"))
    match.arg(object$method, c("One Sample t-test","Welch Two Sample t-test"))

    ## ** compute standard error for the mean/difference in mean
    if(object$method=="One Sample t-test"){
        z <- list(object$args$x)
    }else if(object$method=="Welch Two Sample t-test"){
        z <- list(x=object$args$x, y=object$args$y)
    }

    if(is.null(variance)){
        
        if(type=="estimation"){ ## output the current estimated information
        
            se <- object$stderr
        
        }else if(type=="prediction"){ ## output the predicted information if there was no missing value

            se <- sqrt(Reduce("+",lapply(z, function(iZ){stats::var(iZ, na.rm = TRUE)/length(iZ)})))

        }
        
    }else{

        ## get n
        if(type=="estimation"){
            n <- sapply(z, function(iZ){sum(!is.na(iZ))})
        }else if(type == "prediction"){
            n <- sapply(z, function(iZ){length(iZ)})
        }
        
        ## get variance
        if(is.list(variance)){variance <- unlist(variance)}
        if(object$method=="One Sample t-test"){
            if(length(variance) != 1){
                stop("The number of variance parameters contained in the \'method\' argument do not match the number of groups.\n")
            }
        }else if(object$method=="Welch Two Sample t-test"){
            if(length(variance) != 2){
                stop("The number of variance parameters contained in the \'method\' argument do not match the number of groups.\n")
            }
        }
        
        se <- sqrt(sum(variance/n))        
    }
        
    ## ** export information
    return(as.double(1/se^2))
}

## * getInformation.gls
#' @rdname getInformation
#' @export
getInformation.gls <- function(object, name.coef, data = NULL, details = FALSE,
                               newdata = NULL, variance = NULL, ...){

    requireNamespace("nlme")
    
    ## ** normalize arguments
    if(is.null(data)){
        data <- try(nlme::getData(object), silent = TRUE)
        if(inherits(data,"try-error")){
            stop("Could not retrieve the data used to fit the gls model. \n",
                 "Consider passing the data via the \'data\' argument. \n")
        }
    }
    data <- as.data.frame(data)
    
    if(!is.null(object$modelStruct$varStruct) && !inherits(object$modelStruct$varStruct, what = "varIdent")){
        stop("Only handles \"varIdent\" variance structures \n.")
    }
    if(!is.null(object$modelStruct$corStruct) && !inherits(object$modelStruct$corStruct, what = c("corSymm","corCompSymm"))){
        stop("Only handles \"corSymm\" and \"corCompSymm\" correlation structures \n.")
    }

    if(!is.null(variance) && details){
        stop("Cannot compute information for user-defined covariance values when argument \'details\' is TRUE \n")
    }
    
    ## ** extract elements from model
    ## *** formula
    f.gls <- stats::formula(object)
    outcome.var <- all.vars(f.gls)[1]
    cluster.var <- utils::tail(all.vars(stats::formula(object$modelStruct$corStruct)),1)

    if(length(cluster.var)>0){
        if(is.factor(data[[cluster.var]])){
            data[[cluster.var]] <- droplevels(data[[cluster.var]])
            level.cluster <- levels(data[[cluster.var]])
        }else{
            level.cluster <- sort(unique(data[[cluster.var]]))
        }
    }else{
        level.cluster <- 1:NROW(data)
    }
    n.cluster <- length(level.cluster)
    
    ## *** coefficients
    name.allcoef <- names(stats::coef(object))
    n.allcoef <- length(name.allcoef)
    name.coef <- match.arg(name.coef, name.allcoef)

    ## ** design matrix and covariance pattern
    ## *** at decision
    ## design matrix
    X.decision <- stats::model.matrix(f.gls, data = stats::model.frame(formula = f.gls, data = data, na.action = stats::na.pass))
    ## instead of X <- model.matrix(stats::formula(object), data = data) to keep rows with missing data        

    ## covariance pattern
    resPattern.decision <- .getPattern(object, data = data, variance = NULL)

    ## number of observation per cluster at decision
    resPattern.decision$nobs.vargroup <- stats::setNames(sapply(resPattern.decision$variance.vargroup,NCOL)[resPattern.decision$index.vargroup], names(resPattern.decision$index.vargroup))

    ## *** at interim: full information approach
    ## dataset
    index.interim <- which(rowSums(is.na(data[,all.vars(f.gls)]))==0)
    data.interim <- data[index.interim,,drop=FALSE]

    ## design matrix
    X.interim <- stats::model.matrix(f.gls, data = data)

    ## covariance pattern
    resPattern.interim <- .getPattern(object, data = data.interim, variance = NULL)

    ## number of observation per cluster at interim
    resPattern.interim$nobs.vargroup <- stats::setNames(sapply(resPattern.interim$variance.vargroup,NCOL)[resPattern.interim$index.vargroup], names(resPattern.interim$index.vargroup))

    ## clusters
    if(length(cluster.var)>0){
        if(is.factor(data.interim[[cluster.var]])){
            data.interim[[cluster.var]] <- droplevels(data.interim[[cluster.var]])
            level.cluster.interim <- levels(data.interim[[cluster.var]])
        }else{
            level.cluster.interim <- sort(unique(data.interim[[cluster.var]]))
        }
    }else{
        level.cluster.interim <- 1:NROW(data.interim)
    }
    n.cluster.interim <- length(level.cluster.interim)

    ## *** at interim: complete case approach
    if(resPattern.decision$rep.full==1){
        X.interim.cc <- X.interim
        resPattern.interim.cc <- resPattern.interim
        data.interim.cc <- data.interim
    }else{
        ## dataset (keep all observations for which the cluster has no missing values in the response)
        cluster.cc <- names(resPattern.interim$nobs.vargroup[resPattern.interim$nobs.vargroup==resPattern.decision$rep.full])
        index.cc <- which(data[[cluster.var]] %in% cluster.cc)
        data.interim.cc <- data[index.cc,,drop=FALSE]

        ## design matrix
        X.interim.cc <- stats::model.matrix(f.gls, data = data.interim.cc)

        ## covariance pattern
        resPattern.interim.cc <- .getPattern(object, data = data.interim.cc, variance = NULL)

        ## number of observation per cluster in the complete case
        resPattern.interim.cc$nobs.vargroup <- stats::setNames(sapply(resPattern.interim.cc$variance.vargroup,NCOL)[resPattern.interim.cc$index.vargroup], names(resPattern.interim.cc$index.vargroup))
    }    

    ## *** sanity checks
    if(any(abs(X.decision[index.interim,,drop=FALSE]-X.interim)>1e-10)){ ## Sanity check
        warning("Something went wrong when selecting the data at interim. \n",
                "Could be due to missing values in the regressors. \n")
    }
    if(abs(NROW(X.interim)-NROW(data.interim))>1e-10){ ## Sanity check
        warning("Something went wrong when selecting the data at interim. \n",
                "Could be due to missing values in the regressors. \n")
    }
    if(n.cluster != resPattern.decision$n.cluster){
        warning("Something went wrong when identifying the clusters. \n")
    }
    if(n.cluster.interim != resPattern.interim$n.cluster){
        warning("Something went wrong when identifying the clusters at interim. \n")
    }

    
    ## ** compute information
    ## *** at decision
    ## computed as if all outcome values would be observed
    ## do not handle missing values in covariates
    info.decision <- getInformation(X.decision,
                                    variance = resPattern.decision$variance.vargroup,
                                    index.variance = resPattern.decision$index.vargroup,
                                    index.cluster = resPattern.decision$index.cluster)
    var.decision <- solve(info.decision)[name.coef,name.coef]

    ## *** at interim: full information approach
    info.interim <- getInformation(X.interim,
                                   variance = resPattern.interim$variance.vargroup,
                                   index.variance = resPattern.interim$index.vargroup,
                                   index.cluster = resPattern.interim$index.cluster)
    var.interim <- solve(info.interim)[name.coef,name.coef]
    
    ## *** at interim: complete case approach
    info.interim.cc <- getInformation(X.interim.cc,
                                      variance = resPattern.interim.cc$variance.vargroup,
                                      index.variance = resPattern.interim.cc$index.vargroup,
                                      index.cluster = resPattern.interim.cc$index.cluster)
    var.interim.cc <- solve(info.interim.cc)[name.coef,name.coef]
    
    ## ** sample size
    n.decision <- resPattern.decision$n.cluster ## number of patients at decision (with at least one observable proxy or outcome).
    n.interim <- resPattern.interim$n.cluster ## number of patients at interim (with at least one observation proxy or outcome). Same as decision except in the case of pipeline patients
    n.interim.cc <- resPattern.interim.cc$n.cluster ## number of patients at interim with complete data

    ## ** export
    if(!is.null(attr(data,"df.allobs"))){
        n.all <- as.double(length(unique(attr(data,"df.allobs")[[cluster.var]])))
    }else{
        n.all <- n.decision
    }
    out <- list(decision = list(data = data,
                                X = X.decision,
                                pattern = resPattern.decision,
                                vcov = solve(info.decision)),
                interim = list(data = data.interim,
                               X = X.interim,
                               pattern = resPattern.interim,
                               vcov = solve(info.interim)),
                interim.cc = list(data = data.interim.cc,
                                  X = X.interim.cc,
                                  pattern = resPattern.interim.cc,
                                  vcov = solve(info.interim.cc)),
                sample.size = c(total = n.all, ## all patients in the original dataset: all patient including those who dropped out early and have no observable outcome
                                decision = n.decision,  ## no missing data analysis: all patients with at least one observable outcome (i.e. including those with not yet observed values)
                                interim = n.interim, ## full information analysis: all patients with at least one observed outcome
                                interim.cc = n.interim.cc), ## complete case analysis: all patients with all observed outcomes
                information = c(decision = 1/solve(info.decision)[name.coef,name.coef],
                                interim = 1/solve(info.interim)[name.coef,name.coef],
                                interim.cc = 1/solve(info.interim.cc)[name.coef,name.coef]),
                formula.mean = f.gls, name.coef = name.coef, cluster.var = cluster.var
                )
    if(details){
        if(!is.null(newdata)){
            out$information["decision"] <- getInformation.lmmGSD(c(out, list(fit = object)), newdata = newdata, variance = variance, ...)
        }else if(any(is.na(data[,all.vars(f.gls)]))){ ## better handling of missing value compared to initial calculation
            out$information["decision"] <- getInformation.lmmGSD(c(out, list(fit = object)), newdata = data, variance = variance, ...)
        }
        return(out)
    }else if(is.null(newdata) && is.null(variance)){        
        return(stats::setNames(1/solve(info.interim)[name.coef,name.coef], name.coef))
    }else{
        return(getInformation.lmmGSD(c(out, list(fit = object)), newdata = newdata, variance = variance, ...))
    }
    return(out)
}

## * getInformation.lmmGSD
#' @rdname getInformation
#' @export
getInformation.lmmGSD <- function(object, newdata = NULL, variance = NULL, timevar = NULL, ...){
    
    ## ** normalize arguments
    if(is.null(newdata)){return(object$info)}
    if(!is.null(variance)){
        if(!is.list(variance)){
            variance <- list(variance)
        }
        if(any(sapply(variance,is.matrix)==FALSE)){
            stop("Argument \"variance\" must contain a list of matrices. \n")
        }
        
    }

    cluster.var <- object$cluster.var
    f.gls <- object$formula.mean
    name.regressor <- all.vars(stats::update(f.gls,0~.))
    name.coef <- object$name.coef
    name.outcome <- setdiff(all.vars(f.gls), name.regressor)

    ## ** design matrix
    index.newdata <- which(rowSums(is.na(newdata[,all.vars(f.gls)]))==0)
    newdataRed <- newdata[index.newdata,,drop=FALSE]
    X.newdata <- stats::model.matrix(f.gls, data = newdataRed)
    if(any(abs(X.newdata - stats::model.matrix(f.gls, data = newdata))>1e-10)){
        warning("Something when wrong when subsetting the design matrix.")
    }
    
    ## ** covariance pattern
    resPattern.newdata <- .getPattern(object$fit, data = newdataRed, variance = variance)
    ## number of observation per cluster at decision
    resPattern.newdata$nobs.vargroup <- stats::setNames(sapply(resPattern.newdata$variance.vargroup,NCOL)[resPattern.newdata$index.vargroup], names(resPattern.newdata$index.vargroup))

    ## ** weights
    if(NROW(X.newdata) != NROW(newdata)){

        ## find time variable
        if(is.null(timevar)){
            if(!is.null(object$fit$variable)){
                timevar <- object$fit$variable["time"]
            }else{
                stop("Unknown time variable - consider specifying the argument \'timevar\'. \n")
            }
        }
        ## regressors where there are missing values
        name.regressorNA <- name.regressor[sapply(name.regressor, function(iR){any(is.na(newdata[[iR]]))})]
        ## regressors where there are no missing values
        name.regressorNNA <- setdiff(name.regressor,name.regressorNA)
        ## regressors where there are no missing values and that take constant values within clusters (used to reshape)
        test.cstRegressorCluster <- sapply(name.regressorNNA, function(iReg){ all(tapply(newdata[[iReg]], newdata[[cluster.var]], function(iX){length(unique(iX))})==1) })
        name.regressorNNA2 <- name.regressorNNA[test.cstRegressorCluster]
        ## move dataset to wide format
        newdataW <- stats::reshape(newdata[,c(cluster.var, all.vars(f.gls))], direction = "wide", idvar = c(cluster.var,name.regressorNNA2),
                                   timevar = timevar, v.names = setdiff(all.vars(f.gls),name.regressorNNA2))
        newdataW <- newdataW[order(newdataW[[cluster.var]]),]
        for(iCol in paste(timevar,unique(newdata[[timevar]]),sep=".")){
            if(length(stats::na.omit(unique(newdataW[[iCol]])))>1){
                stop("Time column should only take a single value. \n")
            }
            newdataW[[iCol]][is.na(newdataW[[iCol]])] <- unique(stats::na.omit(newdataW[[iCol]]))
        }
        ## individual with complete follow-up or who had the opportunity to have complete follow-up
        Uid.NA <- as.character(unique(newdata[rowSums(is.na(newdata[,all.vars(f.gls)]))>0,cluster.var]))
        Uid.NNA <- setdiff(as.character(newdata[[cluster.var]]),Uid.NA)
        ## identify reference dataset (no missing value)
        newdataW.NNA <- newdataW[newdataW[[cluster.var]] %in% Uid.NNA,,drop=FALSE]
        ## identify columns with missing values
        name.colNA <- names(which(colSums(is.na(newdataW))>0))
        ## identify missing data patterns to be added to the information
        newdataW.NA <- newdataW[newdataW[[cluster.var]] %in% Uid.NA,,drop=FALSE] ## only keep individuals with incomplete follow-up
        test.newpatternNA <- do.call(cbind,lapply(name.colNA, function(iX){is.na(newdataW.NA[[iX]])}))
        colnames(test.newpatternNA) <-  name.colNA
        vec.newpatternNA <- interaction(as.data.frame(do.call(cbind,lapply(name.colNA, function(iX){ifelse(is.na(newdataW.NA[[iX]]),iX,"")}))),sep="_",drop=TRUE)
        vec.newpatternNA <- factor(vec.newpatternNA, levels(vec.newpatternNA), gsub("_$","",gsub("^_","",levels(vec.newpatternNA))))
        vec.UnewpatternNA <- as.character(unique(vec.newpatternNA))
        n.newpatternNA <- length(vec.UnewpatternNA)
        
        ## generate weights
        weight <- rep(1, NROW(X.newdata))
        
        for(iP in 1:n.newpatternNA){ ## iP <- 1
            iPattern <- vec.UnewpatternNA[iP]
            ## covariates with observed values in this pattern
            iCov.pattern <- names(which(!colSums(test.newpatternNA[vec.newpatternNA==iPattern,,drop=FALSE])>0))
                
            ## stratify on randomization group
            for(iZ in unique(newdataW.NA$Z)){ ## iZ <- 1
                ## lines corresponding to the missing values
                iId <- newdataW.NA[which( (vec.newpatternNA==iPattern) * (newdataW.NA$Z==iZ) == 1),cluster.var]
                iNewdata <- newdata[newdata[[cluster.var]] %in% iId & rowSums(is.na(newdata[,c(name.regressor,name.outcome)]))>0,,drop=FALSE] ## lines with missing value for this pattern
                iTime <- unique(iNewdata[[timevar]])
                ## control lines matching the one with missing values
                if(length(iCov.pattern)>0){
                    iKeep.cluster <- newdataW.NNA[[cluster.var]][which(rowSums(is.na(newdataW.NNA[,iCov.pattern,drop=FALSE]))==0)]
                    iIndex <- which((newdataRed$Z==iZ) * (newdataRed$id %in% iKeep.cluster) * (newdataRed[[timevar]] %in% iTime) == 1)
                }else{
                    iIndex <- intersect(which(newdataRed$Z==iZ), which(newdataRed[[timevar]] %in% iTime))
                }
                ## iDataW <- newdataRed[iIndex,]
                ## update weigths
                weight[iIndex] <- weight[iIndex] + NROW(iNewdata)/length(iIndex)
                ## cat(sum(weight), " ", NROW(iNewdata),"\n")
                ## technically numerator and denominator should be time specific, but taking the (non time specific) ratio leads to same result
                
            }
        }

        if(abs(sum(weight)-NROW(newdata)) > 1e-10){
            stop("Something went wrong when attributing weights to account for not yet observed values in the information.")
        }
    }else{
        weight <- NULL
    }

    ## ** predict information
    info.newdata <- getInformation(X.newdata,
                                   variance = resPattern.newdata$variance.vargroup,
                                   index.variance = resPattern.newdata$index.vargroup,
                                   index.cluster = resPattern.newdata$index.cluster,
                                   weight = weight)
    var.newdata <- solve(info.newdata)[name.coef,name.coef]

    ## ** export
    out <- 1/var.newdata 
    return(out)

}


## * .getPattern
.getPattern <- function(object, data, variance = NULL){
    requireNamespace("nlme")
    
    ## ** reinitialize correlation/variance structure according to data
    if(!is.null(object$modelStruct$corStruct)){
        corStructNew <- do.call(utils::getFromNamespace(class(object$modelStruct$corStruct)[1], ns = "nlme"),
                                args = list(form = stats::formula(object$modelStruct$corStruct),
                                            value = stats::coef(object$modelStruct$corStruct, unconstrained = FALSE))
                                )
        corStructNew <- nlme::Initialize(corStructNew, data = data)
        corgroups <- nlme::getGroups(corStructNew)
    }
    if(!is.null(object$modelStruct$varStruct)){
        varStructNew <- do.call(utils::getFromNamespace(class(object$modelStruct$varStruct)[1], ns = "nlme"),
                                args = list(form = stats::formula(object$modelStruct$varStruct),
                                            value = stats::coef(object$modelStruct$varStruct, unconstrained = FALSE))
                                )
        varStructNew <- nlme::Initialize(varStructNew, data = data)
        vargroups <- nlme::getGroups(varStructNew)
    }

    ## ** number and position of the clusters among the observations
    if(is.null(object$modelStruct$corStruct)){ ## no correlation structure
        index.cluster <- 1:NROW(data)
        n.cluster <- NROW(data)

    }else{ ## correlation structure
        index.cluster <- stats::setNames(as.numeric(corgroups), corgroups)
        n.cluster <- attr(corStructNew,"Dim")$M

    }

    ## ** extract number and position of unique residual variance-covariance structures
    if(is.null(object$modelStruct$varStruct)){ ## no variance structure
        index.vargroup <- stats::setNames(rep(1,n.cluster), sort(unique(names(index.cluster))))
        n.vargroup <-  1
        Sigma.pattern <- list("sigma") ## not needed

    }else if(is.null(object$modelStruct$corStruct)){ ## no correlation structure
        index.vargroup <- as.numeric(as.factor(vargroups))
        n.vargroup <-  max(index.vargroup)
        Sigma.pattern <- as.list(levels(as.factor(vargroups))) ## not needeed
        
    }else{ ## variance and correlation structure
        ## variance parameter within cluster
        variance.per.cluster <- tapply(X = vargroups, INDEX = corgroups, FUN = function(iVec){list(iVec)}) ## WARNING MAY MESS UP THE ORDER
        variance.per.cluster <- variance.per.cluster[levels(corgroups)] ## PROPERLY REORDER
            
        ## unique variance patterns
        Sigma.pattern <- unique(variance.per.cluster)
        names(Sigma.pattern) <- sapply(Sigma.pattern, paste, collapse = "|")
        n.vargroup <- length(Sigma.pattern)
    
        ## associate each cluster to a variance structure
        index.vargroup <- sapply(variance.per.cluster, function(x){
            as.double(which(unlist(lapply(Sigma.pattern, identical, x))))
        })

    }
    
    ## ** correlation coefficient
    if(!is.null(variance)){
        if(NCOL(variance[[1]])>1 && NROW(variance[[1]])>1){
            rho <- variance[[1]][lower.tri(variance[[1]])]
        }else{
            rho <- 0
        }
        
    }else if(is.null(object$modelStruct$corStruct)){
        rho <- 0
    }else{
        rho <- stats::coef(object$modelStruct$corStruct, unconstrained = FALSE)
    }
    
    ## ** full residual variance-covariance matrix
    rep.full <- max(sapply(Sigma.pattern, length))
    Sigma.pattern.full <- Sigma.pattern[sapply(Sigma.pattern, length) == rep.full]
    names(Sigma.pattern.full) <- sapply(Sigma.pattern.full, paste, collapse = "|")
    variance <- .getResVcov(object, variance = variance, rho = rho, Sigma.pattern = Sigma.pattern.full, rep.full = rep.full)

    ## ** (missing data) residual variance-covariance matrix    
    if(is.null(object$modelStruct$varStruct) || is.null(object$modelStruct$corStruct)){ ## no variance or correlation structure
        variance.vargroup <- variance
    }else{ ## variance and correlation structure
        name.vargroup <- names(Sigma.pattern)
        variance.vargroup <- stats::setNames(vector(mode = "list", length = n.vargroup), name.vargroup)
        MSigma.pattern.full <- do.call(rbind,Sigma.pattern.full)

        for(i.vargroup in 1:n.vargroup){ ## i.vargroup <- 2
            iname.vargroup <- name.vargroup[i.vargroup]
            if(iname.vargroup %in% names(Sigma.pattern.full)){
                variance.vargroup[[i.vargroup]] <- variance[[iname.vargroup]]
            }else{
                test <- t(apply(MSigma.pattern.full, 1, `%in%`, Sigma.pattern[[iname.vargroup]])) ## t() because apply(,1,) returns the transposed version
                index <- which.max(rowSums(test))
                variance.vargroup[[i.vargroup]] <- variance[[names(index)]][test[index,],test[index,],drop=FALSE]
            }
        }
    }

    ## ** export
    return(list(index.cluster = index.cluster, ## index of the cluster per observation
                index.vargroup = index.vargroup, ## index of the variance pattern per cluster
                n.cluster = n.cluster, ## number of clusters
                n.vargroup = n.vargroup, ## number of variance patterns
                variance.vargroup = variance.vargroup, ## variance-covariance matrix associated to each variance pattern
                rho =  as.double(rho),
                rep.full = rep.full)) ## number of observation per cluster when no missing data
}

## * .getResVcov
.getResVcov <- function(object, variance, rho, Sigma.pattern, rep.full){

    ## ** build or check variance covariance patterns for fully observed cluster
    if(is.null(variance)){
        if(is.null(object$modelStruct$varStruct) && is.null(object$modelStruct$corStruct)){
            variance <- list(matrix(stats::sigma(object)^2, nrow = 1, ncol = 1))
        }else if(is.null(object$modelStruct$varStruct)){
            variance <- unclass(nlme::getVarCov(object))
        }else if(is.null(object$modelStruct$corStruct)){
            variance <- as.list((stats::sigma(object)*stats::coef(object$modelStruct$varStruct, allCoef = TRUE, unconstrained = FALSE))^2)
        }else{
            vec.sigma <- stats::sigma(object)*stats::coef(object$modelStruct$varStruct, allCoef = TRUE, unconstrained = FALSE)
            Mcor <- matrix(1, nrow = rep.full, ncol = rep.full)
            Mcor[lower.tri(Mcor)] <- rho
            Mcor[upper.tri(Mcor)] <- t(Mcor)[upper.tri(Mcor)]
            variance <- stats::setNames(lapply(Sigma.pattern, function(x){
                tcrossprod(vec.sigma[x]) * Mcor
            }), names(Sigma.pattern))
        }
    }else{

        ## check format and names
        if(length(variance) != length(Sigma.pattern)){
            stop("Incorrect \'variance\' argument (length = ",length(variance),"). \n",
                 length(Sigma.pattern)," variance-covariance patterns detected (levels: \"",paste0(names(Sigma.pattern), collapse="\" \""),"\"). \n")
        }
        if(!is.null(names(variance))){
            if(any(names(variance) %in% names(Sigma.pattern) == FALSE)){
                stop("Incorrect \'variance\' argument (names:  \"",paste0(names(variance), collapse="\" \""),"\"). \n",
                     "variance-covariance patterns detected (levels: \"",paste0(names(Sigma.pattern), collapse="\" \""),"\"). \n")
            }else{
                variance <- variance[names(Sigma.pattern)]
            }
        }else{
            names(variance) <- names(Sigma.pattern)
        }

        
        ## check correlation
        ls.rho <- lapply(variance[sapply(variance,length)>1], function(x){stats::cov2cor(x)[lower.tri(x)]})
        if(length(ls.rho)>1){
            M.rho <- do.call(rbind, ls.rho)
            test <- apply(M.rho, MARGIN = 2, FUN = function(x){max(abs(diff(x)))})
            if(test>1e-10){
                stop("When specificy the argument \"variance\", the correlation should be the same for all matrices. \n")
            }
        }
    }

    return(variance)
}


## * .nProxy
## get the number of observations missing at interim for the parameter of interest but not missing for the proxy
.nProxy <- function(n.interim, n.interim.cc,
                    X.decision, X.interim, name.coef, index.interim,
                    data, data.interim, cluster.var){

    ## ** if only one parameter, then any missing value is missing for the parameter of interest
    if(NCOL(X.decision)==1){
        return(n.interim - n.interim.cc)
    }
    ## else it could be that some "lines" of the design matrix are not useful to estimate the parameter of interest
    ## and are only useful to estimate nuisance parameters.
    
    ## ** do we have full information about the patients included at interim?
    ## if yes we are comparing at the cluster level so we can indeed substract the sample size
    index.interim.full <- which(data[[cluster.var]] %in% unique(data.interim[[cluster.var]]))
    if(length(setdiff(index.interim.full, index.interim))==0){
        return(n.interim - n.interim.cc)
    }
    ## else could be that the added lines when doing full information do not benefit the parameter of interest
    ## and therefore the sample size is the same as complete case analysis.
    
    ## ** Is having full data for each cluster at interim better than the observed data at interim for estimating the parameter of interest?
    ## i.e. does score of the parameter of interest contains a linear combination of the score of other parameters, so that their contribution to the score of the parameter of interest will always be 0.
    ## example Score = [Score_alpha & Score_beta] where Score_alpha = [S_1 \\ S_2] and Score_beta = [S_2 \\ 0] so sum(S_2) must be 0 and therefore the last observations do not contribute to Score_alpha
    X.new <- X.decision[setdiff(index.interim.full, index.interim),,drop=FALSE]
    scorefit <- stats::lm.fit(x = X.new[,setdiff(colnames(X.new),name.coef),drop=FALSE], y = X.new[,name.coef])
    if(any(abs(scorefit$residuals)>1e-10)){ ## not a linear combination
        return(n.interim - n.interim.cc)
    }
    
    ## check whether,  at interim (full information), the score of the parameter of interest contains the score of the linear combination, via the design matrix
    combin <- stats::na.omit(stats::coef(scorefit))
    X.combin <- Reduce("+",lapply(1:length(combin), function(x){
        X.interim[,names(combin)[x],drop=FALSE]*combin[x]
    }))
    test.0 <- (X.combin-X.interim[, name.coef,drop=FALSE])[which(abs(X.combin)>1e-10)]

    if((all(abs(test.0)<1e-10))){
        return(0) ## number of patients without the outcome measurment and only the proxy measurement 
    }else{
        return(n.interim - n.interim.cc) ## number of patients without the outcome measurment and only the proxy measurement
    }
}

######################################################################
### getInformation.R ends here
