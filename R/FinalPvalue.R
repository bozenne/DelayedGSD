## * FinalPvalue (documentation)
#' @title P-value at Decision
#' @description  Compute the p-value at the end of the trial.
#' @param Info.d Information at all decision analyses up to stage where study was stopped (including the final analysis if stopped at final analysis)
#' @param Info.i Information at all interim analyses up to stage where study was stopped (excluding final analysis if stopped at final analysis)
#' @param ck,ck.unrestricted decision boundaries for all decision analyses up to stage where study was stopped (should include final boundary if stopped at final analysis).
#' ck is possibly with restriction (when cNotBelowFixedc=TRUE) and ck.unrestricted always without.
#' @param lk lower bounds up to stage where study was stopped (should include interim boundaries only, not the boundary of the final analysis if study continued to the end)
#' @param uk upper bounds up to stage where study was stopped (should include interim boundaries only, not the boundary of the final analysis if study continued to the end)
#' @param reason.interim motivation for stopping or continuing at interim. Use to handle special cases (skipped interim because reach Imax, ...)
#' @param kMax maximum number of analyses
#' @param delta true effect under which to calculate the probability (should always be 0 for p-value, only change for calculation of CI)
#' @param estimate naive estimate (e.g. using  ML or REML).
#' @param method  method 1, 2 or 3
#' @param bindingFutility [logical] whether the futility stopping rule is binding.
#' @param cNotBelowFixedc [logical] whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
#' @param continuity.correction [logical] whether to add the probability of stopping between ck and ck.uncorrected to ensure continuity of the p-value across stages.
#' When used the p-value will always be greater than this probability of stopping bettwen ck and ck.uncorrected.
#'
#' @examples
#' library(mvtnorm)
#' 
#' if(require(rpact)){
#'
#' #### example 1 ####
#' ## simulate data for a given design
#' kMax <- 3 
#' design <- getDesignGroupSequential(kMax=kMax, sided = 1, alpha = 0.025, beta = 0.2,
#'                                   informationRates = c(0.5,0.75,1),
#'                                   typeOfDesign="asKD",
#'                                   typeBetaSpending="bsKD",gammaA=2,gammaB=2)
#' efficacyBound <- design$criticalValues
#' futilityBound <- design$futilityBounds
#' 
#' res <- getDataset(sampleSizes1 = c(40,20),sampleSizes2=c(40,20),means1=c(1.5,3),means2=c(0.5,-1),stDevs1=c(2,1.65),stDevs2=c(2,1.65))
#'
#' ## extract key informations
#' vec.stage <- res$stages ##  1 1 2 2
#' vec.mean <- res$overallMeans ##   1.5 0.5 2.0 0.0
#' vec.std <- res$overallStDevs ##  2.000000 2.000000 2.007307 2.007307
#' vec.n <- res$overallSampleSizes ##  40 40 60 60
#' vec.z <- vec.mean/(vec.std/sqrt(vec.n)) ## 4.743416 1.581139 7.717771 0.000000
#' estimate <- abs(diff(vec.mean[vec.stage==res$.kMax])) ##  2
#' 
#' Info.var <- as.double(1/tapply(vec.std^2/vec.n,vec.stage,sum)) ## 5.0000 7.4455 
#' statistic <- estimate*sqrt(Info.var[2]) ## 5.457289 
#' 
#' Info.cor <- (rep(1,res$.kMax) %o% sqrt(Info.var)) / t(rep(1,res$.kMax) %o% sqrt(Info.var))
#' Info.cor[upper.tri(Info.cor)] <- 0
#' Info.cor <- Info.cor + t(Info.cor)
#' diag(Info.cor) <- 1
#' 
#' ## p-value by hand
#' ## - more extreme 1: stop at first interim for efficacy
#' pval1 <- pmvnorm(lower = efficacyBound[1], upper = Inf, mean=0, sigma= Info.cor[1,1,drop=FALSE]) 
#' ## - more extreme 1: stop at second interim for efficacy with a larger effect
#' pval2 <- pmvnorm(lower = c(futilityBound[1], statistic), upper = c(efficacyBound[1],Inf), mean=c(0,0), sigma = Info.cor) 
#' ## global
#' pval1 + pval2 ## 0.00625 
#' 
#' ## p-value using FinalPvalue
#' FinalPvalue(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 3, estimate = estimate)
#' 
#' ## p-value using rpact
#' 
#' ests <- getAnalysisResults(design, res, normalApproximation=TRUE)
#'
#' ## Final p-value                                : NA, 0.00625, NA 
#' ##  Final CIs (lower)                            : NA, 0.2414, NA 
#' ##  Final CIs (upper)                            : NA, 2.001, NA 
#' ##  Median unbiased estimate                     : NA, 1.121, NA
#'
#' #### example 2 (see Wassmer 2016, Group sequential and confirmatory adaptive designs in clinical trials, section 4.1, page 87-88) ####
#' vec.stage <- 1:2 ##  1 1 2 2
#' vec.std <- c(1,1)
#' vec.z <- c(1,3)
#' vec.n <- c(22,44)
#' estimate <- vec.z/sqrt(vec.n[2])
#'
#' efficacyBound <- c(2.361,2.361)
#' futilityBound <- c(-2.361,-2.361)
#' 
#' Info.var <- vec.n
#' statistic <- estimate*sqrt(Info.var[2]) ## 5.457289 
#' 
#' Info.cor <- (rep(1,length(vec.stage)) %o% sqrt(Info.var)) / t(rep(1,length(vec.stage)) %o% sqrt(Info.var))
#' Info.cor[upper.tri(Info.cor)] <- 0
#' Info.cor <- Info.cor + t(Info.cor)
#' diag(Info.cor) <- 1
#' 
#' ## p-value by hand
#' ## - more extreme 1: stop at first interim for efficacy
#' pval1 <- pmvnorm(lower = efficacyBound[1], upper = Inf, mean=0, sigma= Info.cor[1,1,drop=FALSE]) 
#' ## - more extreme 1: stop at second interim for efficacy with a larger effect
#' pval2 <- pmvnorm(lower = c(futilityBound[1], statistic[2]), upper = c(efficacyBound[1],Inf), mean=c(0,0), sigma = Info.cor) 
#' ## global
#' pval1 + pval2 ## 0.009819342 
#' 
#' ## p-value using FinalPvalue
#' FinalPvalue(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 4, estimate = estimate)
#' 
#' ## confidence interval using FinalPvalue
#' FinalCI(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 4, estimate = estimate) ## [1] 0.07371176 0.72920245
#'
#' ## library(microbenchmark)
#' ## microbenchmark(optimise = FinalCI(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 4, estimate = estimate, optimizer = "optimise"),
#' ##               uniroot = FinalCI(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 4, estimate = estimate, optimizer = "uniroot") 
#' ## )
#' ## Unit: milliseconds
#' ##      expr      min       lq     mean   median       uq      max neval cld
#' ##  optimise 27.58699 30.60128 32.20530 32.24507 33.47103 44.99762   100   b
#' ##   uniroot 21.37910 23.32777 24.45567 24.29751 25.16585 29.36898   100  a
#' 
#' }
#' 
#' 
#' ##Example to show that the p-value will be 0.025 if the final test statistic is exactly on the boundary (requires cNotBelowFixedc=F and in case of non-binding futility that the info at decision is same as at interim)
#' myBound <- CalcBoundaries(kMax=3,
#'                           sided=1,
#'                           alpha=0.025,  
#'                           beta=0.2,  
#'                           InfoR.i=c(1/3,2/3,1),
#'                           rho_alpha=2,
#'                           rho_beta=2,
#'                           method=1,       #works both for method=1 and method=2
#'                           cNotBelowFixedc=FALSE,
#'                           bindingFutility=FALSE,
#'                           delta=1,
#'                           InfoR.d=c(1/3,2/3))
#' 
#' FinalPvalue(myBound$Info.d,
#'             myBound$Info.i,
#'             myBound$ck,
#'             myBound$lk,
#'             myBound$uk,
#'             myBound$kMax,
#'             estimate = myBound$uk[myBound$kMax]/sqrt(myBound$Info.max),
#'             bindingFutility=F,
#'             futility2efficacy=T) 

## * FinalPvalue (code)
#' @export
FinalPvalue <- function(Info.d,  
                        Info.i,  
                        ck,
                        ck.unrestricted,   
                        lk,  
                        uk,
                        reason.interim,
                        kMax, 
                        delta=0,  
                        estimate,
                        method,  
                        bindingFutility,
                        cNotBelowFixedc,
                        continuity.correction){
    
  
    requireNamespace("mvtnorm")
  
    ## ** reconstruct test statistic
    k <- length(Info.d)
    statistic <- estimate * sqrt(Info.d[k])
                                   
    if(continuity.correction == 3){
        statistic <- statistic - (ck.unrestricted[k] - ck[k])
        statistic_nocor <- statistic
    }else{
        statistic_nocor <- statistic
        if(cNotBelowFixedc  && (statistic < ck[k]) && (statistic > ck.unrestricted[k])){
            statistic <- ck.unrestricted[k] ## correaction to ensure consistency between rejection and p-value (we do not reject in this interval)
        }
    }
                                        #}
    ## ** modify information matrix to deal with decreasing information (special case)
    ## For an interim where we do not stop (nor skip): decreasing information between interim and hypothetical decision
    Info.d2 <- Info.d
    Info.i2 <- Info.i
    if(k>1 && any( (Info.i[1:(k-1)]>Info.d[1:(k-1)]) & reason.interim[1:(k-1)]!="decreasing information") ){
        index.pb <- which( (Info.i[1:(k-1)]>Info.d[1:(k-1)]) & reason.interim!="decreasing information")
        ## keep the correlation interim(iState),decision(iStage), rescale the information at decision accordingly
        Info.d2[index.pb] <- Info.i[index.pb] * Info.i[index.pb]/Info.d[index.pb]        
    }
    ## For an interim where we stop: decreasing information between interim and decision (or final)
    if(k<kMax && Info.i[k]>Info.d[k]){
        ## keep the correlation interim(iState),decision(iStage), rescale the information at decision accordingly
        ## note that we will never skip this interim since we stopped and went for decision, i.e. it must have greater information than previous interim
        Info.d2[k] <- Info.i[k] * Info.i[k]/Info.d[k]        
    }else if(k==kMax && any(Info.i > Info.d[kMax])){
        ## keep the correlation interim(iState),decision(iStage), rescale the information at decision accordingly
        Info.i2[kMax] <- max(Info.i[-kMax]) * max(Info.i[-kMax])/Info.i[kMax]
    }
    ## Previous version:
    ## we will ignore the decision analysis and simply calculate the probability of stopping for efficacy.
    ## This is correct for Methods 1 and 2 if c can be chosen freely (under the null).
    ## It is conservative otherwise as it is more likely to flip from efficacy to futility, hence the p-value will be larger than it should be (under the null)

    ## ** reconstruct information matrix
    if(k==kMax){
        I_all <- c(as.vector(rbind(Info.i, Info.d[-kMax])),Info.d[kMax])  ## vectorize in the following order: |interim 1| |decision 1| |interim 2| ... |interim k| |decision k| |decision final|
        I_all2 <- c(as.vector(rbind(Info.i2, Info.d2[-kMax])),Info.d2[kMax])  ## same but non-decreasing information
        index_tempo <- as.vector(rbind(rep(1,length(Info.i)), rep(2,length(Info.d[-kMax]))))
        index_interim <- which(index_tempo==1) ## 1, 3, ...
        index_decision <- which(index_tempo==2) ## 2, 4, ...
        index_final <- length(I_all)
        critval <- ck[kMax]
    }else{
        I_all <- as.vector(rbind(Info.i, Info.d))  ## vectorize in the following order: |interim 1| |decision 1| |interim 2| ... |interim k| |decision k|
        I_all2 <- as.vector(rbind(Info.i, Info.d2))  ## same but non-decreasing information
        index_tempo <- as.vector(rbind(rep(1,length(Info.i)), rep(2,length(Info.d))))
        index_interim <- which(index_tempo==1) ## 1, 3, ...
        index_decision <- which(index_tempo==2) ## 2, 4, ...
        index_final <- NA
        critval <- ck[k]
    }

    m <- length(I_all)
    sigmaZm <- diag(1,m)
    for(i in 1:m){
        for(j in i:m){
            sigmaZm[i,j] <- sqrt(I_all2[i]/I_all2[j])
            sigmaZm[j,i] <- sigmaZm[i,j]
        }
    }

  
    ## ** compute p-value 
    pval <- 0

    ## need to store original bounds, even in non-binding case to account for probability of switching from futility to efficacy at decision analysis k
    lk_orig <- lk  
    if(!bindingFutility){
        lk[1:min(k,kMax-1)] <- -Inf
    }

    ## under H0 or H1
    theta <- delta * sqrt(I_all)

    ## *** previous stages
    if(k>1){
        for(iStage in 1:(k-1)){ ## iStage <- 1
           
            if(reason.interim[iStage]=="decreasing information"){ ## special case: decreasing information from on interim analysis to another interim analysis
                next ## skip interim
            }else{ ## normal case and special case decreasing information between interim and decision analysis
                      
                if(iStage==1){
                    iSeq_interimM1 <- NULL
                }else{ ## index of previous interim that are not skipped
                    iSeq_interimM1 <- intersect(1:(iStage-1),which(reason.interim[1:(iStage-1)]!="decreasing information"))
                }
                iLk <- lk[iSeq_interimM1]
                iUk <- uk[iSeq_interimM1]
                iIndex <- c(index_interim[iSeq_interimM1], index_interim[iStage], index_decision[iStage])
                iSigma <- sigmaZm[iIndex,iIndex,drop=FALSE]
                
                ## probability to stop for efficacy at previous interim analysis and conclude efficacy
                pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,uk[iStage],ck.unrestricted[iStage]),  
                                                upper = c(iUk,Inf,Inf),
                                                mean = theta[iIndex],
                                                sigma= iSigma)
                
                if(method%in%c(1,2)){
                    ## probability to stop for futility at previous interim analysis and conclude efficacy
                    pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,-Inf,ck.unrestricted[iStage]),   
                                                    upper = c(iUk,lk_orig[iStage],Inf),
                                                    mean = theta[iIndex],
                                                    sigma= iSigma)
                    
                }
            }
        }
    }
    #continuity correction option 2:
    if((statistic >= critval) & (ck[k]>ck.unrestricted[k]) & continuity.correction==2){                
      shift <- max(0,ck[k]-ck.unrestricted[k])
    }else{
      shift <- 0
    }

    ## ** current stage (decision or final)
    if(k==1){
        iSeq_interimM1 <- NULL
    }else{ ## index of the previous interim that are not skipped
        iSeq_interimM1 <- intersect(1:(k-1),which(reason.interim[1:(k-1)]!="decreasing information")) 
    }
    iLk <- lk[iSeq_interimM1]
    iUk <- uk[iSeq_interimM1]

    if(k==kMax || reason.interim[k]=="Imax reached"){ ## is it the final analysis (interim 1:(k-1), final) or  (interim 1:(k-2), skip interim, decision) 
      
        if(k<kMax){
            iIndex <- c(index_interim[iSeq_interimM1], index_decision[k])
        }else if(k==kMax){
            iIndex <- c(index_interim[iSeq_interimM1], index_final)
        }
      
      if((statistic >= critval) & (ck[k]>ck.unrestricted[k]) & continuity.correction==1){
        ## continuity correction 1 when concluding efficacy at final
        correction <- mvtnorm::pmvnorm(lower = c(iLk,ck.unrestricted[k]),  
                                       upper = c(iUk,ck[k]),
                                       mean = theta[iIndex],
                                       sigma = sigmaZm[iIndex,iIndex,drop=FALSE])
    
        pval <- pval + correction
        attr(pval,"correction") <- correction
        shift <- 0
      }
        ## prob to continue until final analysis and obtain a higher result than the observed
        pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,statistic-shift),  
                                        upper = c(iUk,Inf),
                                        mean = theta[iIndex],
                                        sigma= sigmaZm[iIndex,iIndex,drop=FALSE])

    }else if(method == 3 && reason.interim[k]=="futility"){ ## is it a decision analysis after stopping at interim for futility (method 3)
        
      iIndex_interim <- c(index_interim[iSeq_interimM1], index_interim[k])
      iIndex <- c(iIndex_interim, index_decision[k])
      
      ## prob to continue (and therefore have a more extreme result than concluding futility now)
      pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,lk[k]),  
                                      upper = c(iUk,uk[k]),
                                      mean = theta[iIndex_interim],
                                      sigma = sigmaZm[iIndex_interim,iIndex_interim,drop=FALSE])
      
      ## prob to stop for futility at interim and conclude a more extreme result 
      pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,-Inf,statistic_nocor),   
                                      upper = c(iUk,lk_orig[k],Inf),
                                      mean = theta[iIndex],
                                      sigma = sigmaZm[iIndex,iIndex,drop=FALSE])
      
      
      ## prob to stop for efficacy at interim and conclude more extreme result (use min of critval and statistic as stopping for eff and z_k>c_k will result in efficacy conclusion (always more extreme than concluding futility, which is the case when stopping for futility at interim for Method 3 as flips from fut to eff not allowed))
      pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,uk[k],min(ck.unrestricted[k],statistic_nocor)), #OLD -shift here, note that shift=0 if statistic < critval, so shift is only subtracted for more extreme results that would have resulted in efficacy
                                      upper = c(iUk,Inf,Inf),
                                      mean = theta[iIndex],
                                      sigma = sigmaZm[iIndex,iIndex,drop=FALSE])
      
      #if((statistic >= critval) & (ck[k]>ck.unrestricted[k]) & continuity.correction==1){
      #  ## continuity correction when stopping for efficacy
      #  correction <- mvtnorm::pmvnorm(lower = c(iLk,uk[k],ck.unrestricted[k]),  
      #                                 upper = c(iUk,Inf,ck[k]),
      #                                 mean = theta[iIndex],
      #                                 sigma = sigmaZm[iIndex,iIndex,drop=FALSE])
      #  pval <- pval + correction
      #  attr(pval,"correction") <- correction
      #  shift <- 0
      #}
      
      
    }else { ## is it a decision analysis after stopping at interim (interim 1:1, decision i)
        iIndex_interim <- c(index_interim[iSeq_interimM1], index_interim[k])
        iIndex <- c(iIndex_interim, index_decision[k])

        if((statistic >= critval) & (ck[k]>ck.unrestricted[k]) & continuity.correction==1){
            ## continuity correction when stopping for efficacy
            correction <- mvtnorm::pmvnorm(lower = c(iLk,uk[k],ck.unrestricted[k]),  
                                           upper = c(iUk,Inf,ck[k]),
                                           mean = theta[iIndex],
                                           sigma = sigmaZm[iIndex,iIndex,drop=FALSE])
            if(method%in%c(1,2)){
                correction <- correction + mvtnorm::pmvnorm(lower = c(iLk,-Inf,ck.unrestricted[k]),  
                                                            upper = c(iUk,lk_orig[k],ck[k]),
                                                            mean = theta[iIndex],
                                                            sigma = sigmaZm[iIndex,iIndex,drop=FALSE])
            }
            pval <- pval + correction
            attr(pval,"correction") <- correction
            shift <- 0
        }

        if((statistic < critval)){
            ## prob to continue (and therefore have a more extreme result than stopping now for futility)
            pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,lk[k]),  
                                            upper = c(iUk,uk[k]),
                                            mean = theta[iIndex_interim],
                                            sigma = sigmaZm[iIndex_interim,iIndex_interim,drop=FALSE])
            
            if(method%in%3){
              ## prob to stop for futility at interim and conclude a more extreme result (this should only be added in case of futility for method 3, as for method 3 flips from fut to eff are not allowed and a futile result should not be considered more extreme than a positive result)
              pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,-Inf,statistic),
                                              upper = c(iUk,lk_orig[k],Inf),
                                              mean = theta[iIndex],
                                              sigma = sigmaZm[iIndex,iIndex,drop=FALSE])
            }
        }

        ## prob to stop for efficacy at interim and conclude more extreme result
        pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,uk[k],statistic - shift),  
                                        upper = c(iUk,Inf,Inf),
                                        mean = theta[iIndex],
                                        sigma = sigmaZm[iIndex,iIndex,drop=FALSE])

        if(method%in%c(1,2)){
             
            ## prob to stop for futility at interim and conclude a more extreme result
            pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,-Inf,statistic - shift),   
                                            upper = c(iUk,lk_orig[k],Inf),
                                            mean = theta[iIndex],
                                            sigma = sigmaZm[iIndex,iIndex,drop=FALSE])
                        
        }
    }

    ## ** Give another try if the pvalue is NA
    if(is.na(pval)){
        pval <- FinalPvalue2(Info.d = Info.d,  
                             Info.i = Info.i,  
                             ck = ck,
                             ck.unrestricted = ck,   
                             lk = lk,  
                             uk = uk,
                             reason.interim = reason.interim,
                             kMax = kMax, 
                             delta = delta,  
                             estimate = estimate,
                             method = method,  
                             bindingFutility = bindingFutility,
                             cNotBelowFixedc = cNotBelowFixedc,
                             continuity.correction = continuity.correction)
    }
    
    ## ** export
    return(unname(min(pval,1)))
}

## * FinalPvalue2 (code)
FinalPvalue2 <- function(Info.d,  
                         Info.i,  
                         ck,
                         ck.unrestricted,   
                         lk,  
                         uk,
                         reason.interim,
                         kMax, 
                         delta=0,  
                         estimate,
                         method,  
                         bindingFutility,
                         cNotBelowFixedc,
                         continuity.correction){

    requireNamespace("mvtnorm")

    ## ** extract 
    stage <- length(Info.d)
    reason.interim[is.na(reason.interim)] <- "NA"

    ## ** prepare

    ## statistic
    statistic <- estimate * sqrt(Info.d[stage])
    if(continuity.correction){
        Fstatistic <- statistic - (ck[stage] - ck.unrestricted[stage])
    }else{
        Fstatistic <- statistic
    }

    ## information
    Info.vec <- c(Info.i,Info.d)
    Info.type <- c(rep("interim",length(Info.i)),
                   rep("decision",min(stage,kMax-1)),
                   rep("final",1))[1:length(Info.vec)]
    Info.indexInterim <- which(Info.type=="interim")
    Info.indexDecision <- which(Info.type=="decision")
    Info.indexFinal <- which(Info.type=="final")
        
    Info.matrix <- diag(1, length(Info.vec))
    Info.matrix[lower.tri(Info.matrix)] <- sqrt((1/Info.vec) %*% t(Info.vec))[lower.tri(Info.matrix)]
    Info.matrix[upper.tri(Info.matrix)] <- t(Info.matrix)[upper.tri(Info.matrix)]

    ## futility bounds
    if(bindingFutility){
        lk.continue <- lk
    }else{
        lk.continue <- rep(-Inf,length(lk))
    }

    ## ** evaluate p-value at each stage
    ls.pvalue <- lapply(1:stage, function(iStage){ ## iStage <- 2
        if(reason.interim[iStage]=="decreasing information"){ ## SPECIAL CASE: decreasing information from on interim analysis to another interim analysis

            iOut <- 0
            attr(iOut,"reason") <- "decreasing information"

        }else{

            ## previous interim analyses
            if(iStage==1){
                iSeq_interimM1 <- NULL
            }else{ ## index of previous interim that are not skipped
                iSeq_interimM1 <- intersect(1:(iStage-1),which(reason.interim[1:(iStage-1)]!="decreasing information"))
            }

            if(iStage==kMax || reason.interim[iStage]=="Imax reached"){ ## FINAL OR SPECIAL CASE (decision after skipping interim)
                iIndex <- c(Info.indexInterim[iSeq_interimM1], c(Info.indexDecision,Info.indexFinal)[iStage])

                iOut <- pmvnorm2(lower = c(lk.continue[iSeq_interimM1], Fstatistic),  
                                 upper = c(uk[iSeq_interimM1],                 Inf),
                                 mean = delta * sqrt(Info.vec[iIndex]),
                                 sigma = Info.matrix[iIndex,iIndex,drop=FALSE])

                attr(iOut,"terms") <- stats::setNames(c(as.numeric(iOut),0,0), c("efficacy","reversal","continue"))

            }else{ ## DECISION OR SKIPPED DECISION (possibly after skipped interim due to decreasing information)

                iIndex_interim <- c(Info.indexInterim[iSeq_interimM1], Info.indexInterim[iStage])
                iIndex <- c(iIndex_interim, Info.indexDecision[iStage])

                if(iStage != stage){ ## at prevous stages
                    ## 1- probability to stop for efficacy and conclude efficacy
                    iTerm1 <- pmvnorm2(lower = c(lk.continue[iSeq_interimM1], uk[iStage], ck.unrestricted[iStage]),  
                                       upper = c(uk[iSeq_interimM1],                 Inf,        Inf),
                                       mean = delta * sqrt(Info.vec[iIndex]),
                                       sigma = Info.matrix[iIndex,iIndex,drop=FALSE])

                    ## 2- probability to stop for futility and conclude efficacy
                    if(method %in% 1:2){
                        iTerm2 <- pmvnorm2(lower = c(lk.continue[iSeq_interimM1], -Inf,       ck.unrestricted[iStage]),   
                                           upper = c(uk[iSeq_interimM1],          lk[iStage], Inf),
                                           mean = delta * sqrt(Info.vec[iIndex]),
                                           sigma = Info.matrix[iIndex,iIndex,drop=FALSE])
                    }else{
                        iTerm2 <- 0
                    }

                    ## 3- 
                    iTerm3 <- 0

                }else if(iStage == stage){ ## at the current stage (after interim)

                    ## 1- probability to stop for efficacy and have a more extreme test statistic
                    iTerm1 <- pmvnorm2(lower = c(lk.continue[iSeq_interimM1], uk[iStage], Fstatistic),  
                                       upper = c(uk[iSeq_interimM1],                 Inf,        Inf),
                                       mean = delta * sqrt(Info.vec[iIndex]),
                                       sigma = Info.matrix[iIndex,iIndex,drop=FALSE])

                    ## 2- probability to stop for futility and have a more extreme test statistic
                    ## NOTE: concluding for efficacy with method 3 is more extreme than stopping for futility (since that will lead to concluding futility).
                    ##       this case is therefore excluded
                    
                    if(!(method == 3 & statistic > ck[stage] & reason.interim[stage]!="futility")){
                        ## NOTE: notice the negation ("!") before the condition, i.e. either method 1 or 2
                        ##                                                                or method 3 where we conclude efficacy
                        iTerm2 <- pmvnorm2(lower = c(lk.continue[iSeq_interimM1], -Inf,       Fstatistic),   
                                           upper = c(uk[iSeq_interimM1],          lk[iStage], Inf),
                                           mean = delta * sqrt(Info.vec[iIndex]),
                                           sigma = Info.matrix[iIndex,iIndex,drop=FALSE])
                    }else{
                        iTerm2 <- 0
                    }

                    ## 3- probability to continue instead of stopping now for futility
                    ## NOTE: concluding for futility with method 3 can happen when the test statistic is large because we stopped for futility.
                    ##       this case is therefore added
                    if((statistic < ck[stage]) || (method == 3 && reason.interim[stage]=="futility")){
                        iTerm3 <- pmvnorm2(lower = c(lk.continue[iSeq_interimM1], lk.continue[stage]),
                                           upper = c(uk[iSeq_interimM1],          uk[stage]),
                                           mean = delta * sqrt(Info.vec[iIndex_interim]),
                                           sigma = Info.matrix[iIndex_interim,iIndex_interim,drop=FALSE])
                    }else{
                        iTerm3 <- 0
                    }
                }

                iOut <- unname(iTerm1 + iTerm2 + iTerm3)
                attr(iOut, "terms") <- stats::setNames(c(iTerm1,iTerm2,iTerm3), c("efficacy","reversal","continue"))
            }
        }
        return(iOut)

    })

    ## ** export
    out <- pmin(1,Reduce("+",ls.pvalue))
    attr(out,"error") <- try(sapply(ls.pvalue,attr,"error"), silent = TRUE)
    attr(out,"msg") <- try(sapply(ls.pvalue,attr,"msg"), silent = TRUE)
    attr(out,"terms") <- try(do.call(rbind,lapply(ls.pvalue,attr,"terms")), silent = TRUE)
    return(out)

}

