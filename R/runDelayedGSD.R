### runDelayedGSD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 29 2024 (13:23) 
## Version: 
## Last-Updated: sep 25 2024 (14:29) 
##           By: Brice Ozenne
##     Update #: 96
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * runDelayedGSD (documentation)
##' @description Run Delayed GSD
##'
##' @param data output of GenData (dataset containing all patients).
##' @param boundaries [delayedGSD] output of CalcBoundaries
##' @param N.fw [integer] number of follow-up values per individual.
##' @param PropForInterim [numeric vector of length Kmax-1] percentage of all subjects, once they had the change to have complete follow-up, triggering interim and final analyses.
##' By default equal to \code{InfoR.i}.
##' @param lag [numeric] time lag between stop of recruitment and decision to stop recruitment.
##' @param export.GSD [logical] should a list of \code{delayedGSD} objects used to calculate boundaries and estimates be output as an attribute.
##' @param ... not used
##' 
##' @noRd
##'
##' @examples
##' e.bound <- CalcBoundaries(kMax = 2, alpha = 0.025, beta = 0.1, InfoR.i = c(0.5,1), InfoR.d = c(0.6,1),
##'                           rho_alpha = 2, rho_beta = 2,
##'                           method = 1, cNotBelowFixedc = FALSE, bindingFutility = TRUE, delta = 1)
##' set.seed(19)
##' df.GSD <- nGenData()$d
##' runDelayedGSD(df.GSD, boundaries = e.bound, N.fw = 2, PropForInterim = c(0.5), lag = 0, overrule.futility = FALSE)
##' 

## * runDelayedGSD (code)
##' @export
runDelayedGSD <- function(data, boundaries, N.fw,
                          PropForInterim, lag,
                          overrule.futility,
                          export.GSD = FALSE){

    ## ** check user input
    ## *** boundaries
    if(!inherits(boundaries,"delayedGSD")){
        stop("Incorrect argument \'boundaries\': should inherit from \"delayedGSD\". \n",
             "Consider using the function CalcBoundaries.\n")
    }
    kMax <- boundaries$kMax

    ## *** data
    if(!is.data.frame(data)){
        stop("Argument \'boundaries\' should inherit from \"data.frame\". \n")
    }
    if("id" %in% names(data) == FALSE){
        stop("Argument \'boundaries\' should contain a column \"id\". \n")
    }
    if("Z" %in% names(data) == FALSE){
        stop("Argument \'boundaries\' should contain a column \"Z\". \n")
    }
    if(any(paste0("X",1:kMax) %in% names(data) == FALSE)){
        stop("Argument \'boundaries\' should contain the columns \"X1\" to \"X",kMax,"\". \n")
    }
    if(any(paste0("t",1:kMax) %in% names(data) == FALSE)){
        stop("Argument \'boundaries\' should contain the columns \"t1\" to \"t",kMax,"\". \n")
    }

    ## *** PropForInterim
    if(length(PropForInterim)==kMax && PropForInterim[kMax]==1){
        PropForInterim <- PropForInterim[-kMax]
    }
    if(!is.numeric(PropForInterim) || length(PropForInterim)!=kMax-1){
        stop("Argument \'PropForInterim\' should have length ",kMax-1,". \n")
    }

    ## *** lag
    if(!is.numeric(lag) || length(lag)!=1){
        stop("Argument \'lag\' should a numeric value. \n")
    }

    ## ** prepare
    n.obs <- NROW(data)
    Mn.obs <- matrix(as.numeric(NA), nrow = kMax, ncol = 4,
                     dimnames = list(NULL, c("patients","outcome","pipeline","complete")))
    
    ## ** interim & decision
    time.interim <- data[[paste0("t",N.fw+1)]][ceiling(n.obs*PropForInterim)]
    data.interim <- vector(mode = "list", length = kMax-1)
    lmm.interim <- vector(mode = "list", length = kMax-1)
    gsd.interim <- vector(mode = "list", length = kMax-1)

    start.time <- Sys.time()

    for(iK in 1:(kMax-1)){ ## iK <- 1

        ## *** lmm at interim
        data.interim[[iK]] <- SelectData(data, t = time.interim[iK])
        Mn.obs[iK,c("patients","outcome","pipeline","complete")] <- c(NROW(data.interim[[iK]]),
                                                                      sum(!is.na(data.interim[[iK]][[paste0("t",N.fw+1)]])),
                                                                      sum(is.na(data.interim[[iK]][[paste0("t",N.fw+1)]])),
                                                                      sum(rowSums(is.na(data.interim[[iK]]))==0))
        lmm.interim[[iK]] <- analyzeData(data.interim[[iK]], ddf = "nlme",
                                         data.decision = sum(data$t1 <= time.interim[iK] + lag),
                                         getinfo = TRUE, trace = TRUE)

        ## *** boundary at interim
        if(iK==1){
            gsd.interim[[iK]] <- update(boundaries, delta = lmm.interim[[iK]], trace = FALSE)
        }else{
            gsd.interim[[iK]] <- update(gsd.interim[[iK-1]], delta = lmm.interim[[iK]], trace = FALSE)
        }

        ## *** stop recruitment
        iDecision <- coef(gsd.interim[[iK]], type = "decision")[,paste0("stage ",iK)]
        data.decision <- data[which(data$t1 <= time.interim[iK] + lag),]
        lmm.decision <- analyzeData(data.decision, ddf = "nlme", getinfo = TRUE, trace = TRUE)

        if(iDecision["decision"] == "stop"){

            if((iDecision["comment"]=="futility") && (boundaries$binding == FALSE) && overrule.futility){ 
                ## overrule futility in case of stopping for futility with non-binding bounds
                gsd.interim[[iK]] <- update(gsd.interim[[iK]], overrule.futility = TRUE)
                ## update decision
                iDecision <- coef(gsd.interim[[iK]], type = "decision")[,paste0("stage ",iK)] ## from stop to continue
                ## update information
                gsd.interim[[iK]] <- update(gsd.interim[[iK]], delta = lmm.decision, k = iK, type.k = "decision", trace = FALSE)
            }else{
                Mn.obs[iK+1,c("patients","outcome","pipeline","complete")] <- c(NROW(data.decision),
                                                                                sum(!is.na(data.decision[[paste0("t",N.fw+1)]])),
                                                                                sum(is.na(data.decision[[paste0("t",N.fw+1)]])),
                                                                                sum(rowSums(is.na(data.decision))==0))
                gsd.decision <- update(gsd.interim[[iK]], delta = lmm.decision, trace = FALSE)
                gsd.interim <- gsd.interim[1:iK] ## 'remove' future interim analyses from the results
                break
            }
            
        }else{
            ## update information
            gsd.interim[[iK]] <- update(gsd.interim[[iK]], delta = lmm.decision, k = iK, type.k = "decision", trace = FALSE)
        }
    }

    ## ** final
    if(iDecision["decision"] != "stop"){
        Mn.obs[iK+1,c("patients","outcome","pipeline","complete")] <- c(NROW(data),
                                                                        sum(!is.na(data[[paste0("t",N.fw+1)]])),
                                                                        sum(is.na(data[[paste0("t",N.fw+1)]])),
                                                                        sum(rowSums(is.na(data))==0))                
        lmm.final <- analyzeData(data, ddf = "nlme", getinfo = TRUE, trace = TRUE)
        gsd.decision <- update(gsd.interim[[kMax-1]], delta = lmm.final, trace = FALSE)
    }
    stop.time <- Sys.time()

    ## ** gather results
    ls.out <- lapply(c(gsd.interim[1:iK],list(gsd.decision)), function(iGSD){ ## iGSD <- gsd.interim[[1]]

        iK <- iGSD$stage[,"k"]
        iType <- iGSD$stage[,"type"]
        iBoundary <- coef(iGSD, type = "boundary")

        if(iType == "interim"){
            iConfintAll <- confint(iGSD)
            iConfint <- data.frame(method = c("ML","MUE"),
                                   stage = iK,
                                   type = iType,
                                   coef =  iConfintAll$coef[1],
                                   estimate = c(iConfintAll[iConfintAll$stage == iK & iConfintAll$type == iType,"statistic"],NA),
                                   se = c(iConfintAll[iConfintAll$stage == iK & iConfintAll$type == iType,"se"],NA),
                                   statistic = c(iConfintAll[iConfintAll$stage == iK & iConfintAll$type == iType,"statistic"],NA),
                                   df = c(iConfintAll[iConfintAll$stage == iK & iConfintAll$type == iType,"df"],NA),
                                   p.value = c(iConfintAll[iConfintAll$stage == iK & iConfintAll$type == iType,"p.value"],NA),
                                   lower = c(iConfintAll[iConfintAll$stage == iK & iConfintAll$type == iType,"lower"],NA),
                                   upper = c(iConfintAll[iConfintAll$stage == iK & iConfintAll$type == iType,"upper"],NA))

            iInfoAll <- coef(iGSD, type = "information")
            iInfo <- stats::setNames(iInfoAll[iInfoAll$stage==iK,c("Interim","Interim.pc","Decision","Decision.pc")],
                                     c("current","current.pc","prediction","prediction.pc"))
        }else{
            iConfint <- confint(iGSD, method = c("ML","MUE"))
            iConfint <- iConfint[iConfint$stage == iK & iConfint$type == iType,,drop=FALSE]

            iInfoAll <- coef(iGSD, type = "information")
            iInfo <- stats::setNames(c(unlist(iInfoAll[iInfoAll$stage==iK,c("Decision","Decision.pc")]),rep(as.numeric(NA),2)),
                                     c("current","current.pc","prediction","prediction.pc"))
        }

        if(iType == "interim"){
            iDecision <- coef(iGSD, type = "decision")[,paste0("stage ",iK)]
            iDecisionM1 <- c("decision" = NA, comment = NA)
            iN <- Mn.obs[iK,]
        }else if(iType == "decision"){
            iDecision <- coef(iGSD, type = "decision")[,paste0("stage ",iK," decision")]
            iDecisionM1 <- coef(gsd.interim[[iK]], type = "decision")[,paste0("stage ",iK)]
            iN <- Mn.obs[iK+1,]
        }else if(iType == "final"){
            iDecision <- coef(iGSD, type = "decision")[,paste0("stage ",iK)]
            iDecisionM1 <- coef(gsd.interim[[iK-1]], type = "decision")[,paste0("stage ",iK-1)]
            iN <- Mn.obs[iK,]
        }
        
        iResGSD <- data.frame(time = switch(iType,
                                            "interim" = time.interim[iK],
                                            "decision" = max(data.decision[[paste0("t",N.fw+1)]], na.rm = TRUE),
                                            "final" = max(data[[paste0("t",N.fw+1)]], na.rm = TRUE)),
                              stage = iK,
                              type.stage = iType,
                              n.patients = iN["patients"],
                              n.outcome = iN["outcome"],
                              n.pipeline = iN["pipeline"],
                              n.complete = iN["complete"],
                              uk = switch(iType,
                                          "interim" = iBoundary[iK,"Ebound"],
                                          "decision" = iBoundary[iK,"Cbound"],
                                          "final" = iBoundary[iK,"Cbound"]),
                              lk = switch(iType,
                                          "interim" = iBoundary[iK,"Fbound"],
                                          "decision" = iBoundary[iK,"Cbound"],
                                          "final" = iBoundary[iK,"Cbound"]),
                              statistic = iConfint[1,"statistic"],
                              estimate_ML = iConfint[iConfint$method=="ML","estimate"],
                              p.value_ML = iConfint[iConfint$method=="ML","p.value"],
                              lower_ML = iConfint[iConfint$method=="ML","lower"],
                              upper_ML = iConfint[iConfint$method=="ML","upper"],
                              estimate_MUE = iConfint[iConfint$method=="MUE","estimate"],
                              p.value_MUE = iConfint[iConfint$method=="MUE","p.value"],
                              lower_MUE = iConfint[iConfint$method=="MUE","lower"],
                              upper_MUE = iConfint[iConfint$method=="MUE","upper"],
                              info = unname(iInfo["current"]),
                              infoPC = unname(iInfo["current.pc"]),
                              info.pred = unname(iInfo["prediction"]),
                              infoPC.pred = unname(iInfo["prediction.pc"]),
                              decision = iDecision["decision"],
                              reason = iDecision["comment"],
                              reason.interim = iDecisionM1["comment"])

        return(iResGSD)
    })

    out <- cbind(method = boundaries$method,
                 binding = boundaries$bindingFutility,
                 fixC = boundaries$cNotBelowFixedc,
                 n.max = n.obs,
                 do.call(rbind, ls.out),
                 computation.time = as.numeric(stop.time - start.time))
                 
    rownames(out) <- NULL

    ## ## ** export
    if(export.GSD){
        attr(out, "delayedGSD") <- c(gsd.interim[1:iK],list(gsd.decision))
    }
    return(out)

}

            
        

##----------------------------------------------------------------------
### runDelayedGSD.R ends here
