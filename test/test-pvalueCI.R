### test-pvalueCI.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 27 2022 (09:32) 
## Version: 
## Last-Updated: mar 22 2023 (14:45) 
##           By: Brice Ozenne
##     Update #: 66
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(testthat)
library(rpact)
library(DelayedGSD)
library(mvtnorm)
set.seed(10)

context("Computing p-value/CIs/MUestimate after stopping. \n")

## * rpact
test_that("Compare p-value/CIs/MUestimate to rpact in the non-longitudinal case",{

    kMax <- 3 
    design <- rpact::getDesignGroupSequential(kMax=kMax, sided = 1, alpha = 0.025, beta = 0.2,
                                              informationRates = c(0.5,0.75,1),
                                              typeOfDesign="asKD",
                                              typeBetaSpending="bsKD",gammaA=2,gammaB=2)
    efficacyBound <- design$criticalValues
    futilityBound <- design$futilityBounds
 
    res <- rpact::getDataset(sampleSizes1 = c(40,20),sampleSizes2=c(40,20),means1=c(1.5,3),means2=c(0.5,-1),stDevs1=c(2,1.65),stDevs2=c(2,1.65))
    ests <- rpact::getAnalysisResults(design, res, normalApproximation=TRUE)

    ## extract key informations
    vec.stage <- res$stages ##  1 1 2 2
    vec.mean <- res$overallMeans ##   1.5 0.5 2.0 0.0
    vec.std <- res$overallStDevs ##  2.000000 2.000000 2.007307 2.007307
    vec.n <- res$overallSampleSizes ##  40 40 60 60
    vec.z <- vec.mean/(vec.std/sqrt(vec.n)) ## 4.743416 1.581139 7.717771 0.000000
    estimate <- abs(diff(vec.mean[vec.stage==res$.kMax])) ##  2
 
    Info.var <- as.double(1/tapply(vec.std^2/vec.n,vec.stage,sum)) ## 5.0000 7.4455 
    statistic <- estimate*sqrt(Info.var[2]) ## 5.457289 
 
    Info.cor <- (rep(1,res$.kMax) %o% sqrt(Info.var)) / t(rep(1,res$.kMax) %o% sqrt(Info.var))
    Info.cor[upper.tri(Info.cor)] <- 0
    Info.cor <- Info.cor + t(Info.cor)
    diag(Info.cor) <- 1
 
    ## p-value by hand
    ## - more extreme 1: stop at first interim for efficacy
    pval1 <- pmvnorm(lower = efficacyBound[1], upper = Inf, mean=0, sigma= Info.cor[1,1,drop=FALSE]) 
    ## - more extreme 1: stop at second interim for efficacy with a larger effect
    pval2 <- pmvnorm(lower = c(futilityBound[1], statistic), upper = c(efficacyBound[1],Inf), mean=c(0,0), sigma = Info.cor) 
    ## global
    expect_equal(as.double(pval1 + pval2), 0.00625, tol = 1e-6)
 
    ## p-value using FinalPvalue
    test.p <- FinalPvalue(Info.d = Info.var,
                          Info.i = Info.var,
                          ck = efficacyBound,
                          ck.unrestricted = efficacyBound,
                          lk = futilityBound,
                          uk = efficacyBound[1:2],
                          kMax = 3,
                          estimate = estimate,
                          method = 1,
                          reason.interim = c("",""),
                          bindingFutility = TRUE,
                          cNotBelowFixedc = FALSE,
                          continuity.correction = TRUE)
    expect_equal(as.double(test.p), 0.00625)
    expect_equal(as.double(test.p), as.data.frame(ests)[2,"finalPValues"])
    
    ## confidence interval using FinalPvalue
    test.ci <- FinalCI(Info.d = Info.var,
                       Info.i = Info.var,
                       ck = efficacyBound,
                       ck.unrestricted = efficacyBound,
                       lk = futilityBound,
                       uk = efficacyBound[1:2],
                       kMax = 3,
                       conf.level = 0.95,
                       estimate = estimate,
                       reason.interim = c("",""),
                       method = 1,
                       bindingFutility = TRUE,
                       cNotBelowFixedc = FALSE,
                       continuity.correction = TRUE)
    expect_equal(as.double(test.ci), c(0.53413706, 1.99336202), tol = 1e-4)
    as.data.frame(ests)[2,"finalConfidenceIntervalLowerBounds"] ## 0.2413634
    as.data.frame(ests)[2,"finalConfidenceIntervalUpperBounds"] ## 2.000625

    ## median unbiased estimate
    test.beta <- FinalEstimate(Info.d = Info.var,
                               Info.i = Info.var,
                               ck = efficacyBound,
                               ck.unrestricted = efficacyBound,
                               lk = futilityBound,
                               uk = efficacyBound[1:2],
                               kMax = 3,
                               estimate = estimate,
                               reason.interim = c("",""),
                               method = 1,
                               bindingFutility = TRUE,
                               cNotBelowFixedc = FALSE,
                               continuity.correction = TRUE)
    expect_equal(as.double(test.beta), 1.26709147, tol = 1e-6)
    as.data.frame(ests)[2,"medianUnbiasedEstimates"] ## 1.121087
}) 

## * Wassmer 2016 (Group sequential and confirmatory adaptive designs in clinical trials, section 4.1, page 87-88)

test_that("p-value/CIs/MUestimate with four stage design stopped at interim 2",{

 vec.stage <- 1:2 ##  1 1 2 2
 vec.std <- c(1,1)
 vec.z <- c(1,3)
 vec.n <- c(22,44)
 estimate <- vec.z/sqrt(vec.n[2])

 Info.var <- vec.n
 statistic <- estimate*sqrt(Info.var[2]) ## 5.457289 

 Info.cor <- (rep(1,length(vec.stage)) %o% sqrt(Info.var)) / t(rep(1,length(vec.stage)) %o% sqrt(Info.var))
 Info.cor[upper.tri(Info.cor)] <- 0
 Info.cor <- Info.cor + t(Info.cor)
 diag(Info.cor) <- 1

 ## ** O Brien and Fleming
 ## alpha=0.05, power=0.8, delta = 0.31
 
 uk.Obrien <- c(4.049/sqrt(1:2))
 lk.Obrien <- c(-4.049/sqrt(1:2))
 ck.Obrien <- uk.Obrien
 
 ## p-value by hand
 ## - more extreme 1: stop at first interim for efficacy
 pval1.Obrien <- pmvnorm(lower = uk.Obrien[1], upper = Inf, mean=0, sigma= Info.cor[1,1,drop=FALSE]) 
 ## - more extreme 1: stop at second interim for efficacy with a larger effect
 pval2.Obrien <- pmvnorm(lower = c(lk.Obrien[1], statistic[2]), upper = c(uk.Obrien[1],Inf), mean=c(0,0), sigma = Info.cor) 
 ## global
 expect_equal(round(as.double(pval1.Obrien + pval2.Obrien),4), 0.0014)
 
 ## p-value using FinalPvalue
 pval.Obrien <- FinalPvalue(Info.d = Info.var,
                            Info.i = Info.var,
                            ck = ck.Obrien,
                            ck.unrestricted = ck.Obrien,
                            lk = lk.Obrien,
                            uk = uk.Obrien,
                            kMax = 4,
                            estimate = estimate[2],
                            reason.interim = c("",""),
                            method = 1,
                            bindingFutility = TRUE,
                            cNotBelowFixedc = FALSE,
                            continuity.correction = TRUE)
 expect_equal(as.double(pval.Obrien), 0.001362486, tol = 1e-6)
 
 ## confidence interval using FinalPvalue
 CI.Obrien <- FinalCI(Info.d = Info.var,
                      Info.i = Info.var,
                      ck = ck.Obrien,
                      ck.unrestricted = ck.Obrien,
                      lk = lk.Obrien,
                      uk = uk.Obrien,
                      kMax = 4,
                      estimate = estimate[2],
                      reason.interim = c("",""),
                      method = 1,
                      conf.level = 0.95,
                      bindingFutility = TRUE,
                      cNotBelowFixedc = FALSE,
                      continuity.correction = TRUE)
 expect_equal(as.double(CI.Obrien), c(0.1565499,0.7470000), tol = 1e-3) ## FROM THE BOOK: c(0.157, 0.748)

 ## median unbiased estimate
 estimate.Obrien <- FinalEstimate(Info.d = Info.var,
                                  Info.i = Info.var,
                                  ck = ck.Obrien,
                                  ck.unrestricted = ck.Obrien,
                                  lk = lk.Obrien,
                                  uk = uk.Obrien,
                                  kMax = 4,
                                  estimate = estimate[2],
                                  reason.interim = c("",""),
                                  method = 1,
                                  bindingFutility = TRUE,
                                  cNotBelowFixedc = FALSE,
                                  continuity.correction = TRUE)
 expect_equal(round(as.double(estimate.Obrien),3), 0.452, tol = 1e-3)

 ## ** Pocock case
 uk.Pocock <- rep(2.361,2)
 lk.Pocock <- rep(-2.361,2)
 ck.Pocock <- rep(2.361,2)

 ## p-value by hand
 ## - more extreme 1: stop at first interim for efficacy
 pval1.Pocock <- pmvnorm(lower = uk.Pocock[1], upper = Inf, mean=0, sigma= Info.cor[1,1,drop=FALSE]) 
 ## - more extreme 1: stop at second interim for efficacy with a larger effect
 pval2.Pocock <- pmvnorm(lower = c(lk.Pocock[1], statistic[2]), upper = c(uk.Pocock[1],Inf), mean=c(0,0), sigma = Info.cor) 
 ## global
 expect_equal(round(as.double(pval1.Pocock + pval2.Pocock),4), 0.0098)
 
 ## p-value using FinalPvalue
 pval.Pocock <- FinalPvalue(Info.d = Info.var,
                            Info.i = Info.var,
                            ck = ck.Pocock,
                            ck.unrestricted = ck.Pocock,
                            lk = lk.Pocock,
                            uk = uk.Pocock,
                            kMax = 4,
                            estimate = estimate[2],
                            reason.interim = c("",""),
                            method = 1,
                            bindingFutility = TRUE,
                            cNotBelowFixedc = FALSE,
                            continuity.correction = TRUE)
 expect_equal(as.double(pval.Pocock), 0.009819342, tol = 1e-6) 
 
 ## confidence interval using FinalPvalue
 CI.Pocock <- FinalCI(Info.d = Info.var,
                      Info.i = Info.var,
                      ck = ck.Pocock,
                      ck.unrestricted = ck.Pocock,
                      lk = lk.Pocock,
                      uk = uk.Pocock,
                      kMax = 4,
                      estimate = estimate[2],
                      reason.interim = c("",""),
                      method = 1,
                      bindingFutility = TRUE,
                      cNotBelowFixedc = FALSE,
                      conf.level = 0.95,
                      continuity.correction = TRUE)
 expect_equal(round(as.double(CI.Pocock),3), c(0.074, 0.729), tol = 1e-3)

 ## median unbiased estimate
 estimate.Pocock <- FinalEstimate(Info.d = Info.var,
                                  Info.i = Info.var,
                                  ck = ck.Pocock,
                                  ck.unrestricted = ck.Pocock,
                                  lk = lk.Pocock,
                                  uk = uk.Pocock,
                                  kMax = 4,
                                  estimate = estimate[2],
                                  reason.interim = c("",""),
                                  method = 1,
                                  bindingFutility = TRUE,
                                  cNotBelowFixedc = FALSE,
                                  continuity.correction = TRUE)
 expect_equal(round(as.double(estimate.Pocock),3), 0.419, tol = 1e-3)

})

## * P-value when at boundary value 
test_that("Check consistency between p-value and boundary (1 interim analysis)",{

    ## ** Stop at interim and statistic precisely at boundary at decision
    ## method 1 or 2 (the calculation of the p-value is the same for both methods)
    test <- FinalPvalue(Info.d = 2.6923896504288,
                        Info.i = 2.45494676691245,
                        ck = 1.69415016027832, ck.unrestricted = 1.69415016027832, 
                        lk = 0.992173622977366,
                        uk = 2.27907351772183,
                        kMax = 2, 
                        delta = 0,
                        estimate = 1.69415016027832 / sqrt(2.6923896504288),## 1.64370753556795,
                        reason.interim = c(""),
                        method = 2,
                        bindingFutility = TRUE, 
                        cNotBelowFixedc = FALSE,
                        continuity.correction = TRUE)
    GS <- ErrorSpend(I=2.45494676691245,rho=2,beta_or_alpha=0.025,Info.max=3.646459)
    expect_equal(as.double(test), GS, tol = 1e-5)
    
    ## method 3
    test <- FinalPvalue(Info.d = 12.58797814,  
                        Info.i = 10.60271846,
                        ck = 1.897834, ck.unrestricted = 1.897834,
                        lk = 0.22424864,  
                        uk = 2.52429353,  
                        kMax = 2, 
                        delta = 0,  
                        estimate = 1.897834 / sqrt(12.58797814),
                        reason.interim = c(""),
                        method = 3,
                        bindingFutility = TRUE,
                        cNotBelowFixedc = FALSE,
                        continuity.correction = TRUE)
    GS <- ErrorSpend(I=10.60271846,rho=2,beta_or_alpha=0.025,Info.max=22.71478)
    expect_equal(as.double(test), GS, tol = 1e-5)
    
    ## ** Continue at interim and and statistic precisely at boundary at final (method 1, cNotBelowFixedc)
    ## method 1 or 2 (the calculation of the p-value is the same for both methods)
    test <- FinalPvalue(Info.d = c(2.75046695782793, 4.30612228613852),
                        Info.i = c(2.5589317424577),
                        ck = c(1.7253424172699,2.01974814127831), ck.unrestricted = c(1.7253424172699,2.01974814127831),
                        lk = c(1.11127641794951),
                        uk = c(2.2459290489148),
                        kMax = 2,
                        delta = 0, 
                        estimate = 2.01974814127831 / sqrt(4.30612228613852),
                        reason.interim = c("",""),
                        method = 1,
                        bindingFutility = TRUE, 
                        cNotBelowFixedc = FALSE,
                        continuity.correction = TRUE)
    expect_equal(as.double(test), 0.025, tol = 1e-5)
    
    ## method 3
    test <- FinalPvalue(Info.d = c(12.58797814, 24.67092824),  
                        Info.i = c(10.60271846),
                        ck = c(1.897834, 1.99362293), ck.unrestricted = c(1.897834, 1.99362293),
                        lk = c(0.22424864),  
                        uk = c(2.52429353),  
                        kMax = 2, 
                        delta = 0,  
                        estimate = 1.99362293 / sqrt(24.67092824),
                        reason.interim = c("",""),
                        method = 3,
                        bindingFutility = TRUE,
                        cNotBelowFixedc = FALSE,
                        continuity.correction = TRUE)
    expect_equal(as.double(test), 0.025, tol = 1e-5)

})

test_that("Check consistency between p-value and boundary (2 interim analysis)",{

  
    for(iMethod in 1:3){ ## iMethod <- 1
      
      if(iMethod==3){
        cNotBelowFixedc <- TRUE
      } else {
        cNotBelowFixedc <- FALSE
      } 
      
      ## for binding futility rule
        myBound <- suppressWarnings(CalcBoundaries(kMax = 3,
                                                   alpha = 0.025, 
                                                   beta = 0.1, 
                                                   InfoR.i = c(3.5,6.75)/12,
                                                   rho_alpha = 1.345,
                                                   rho_beta = 1.345,
                                                   method = iMethod, ## has been changed from 2 to 1
                                                   cNotBelowFixedc = cNotBelowFixedc,
                                                   bindingFutility = TRUE,
                                                   delta = 1,
                                                   InfoR.d = c(5.5,8.75,12)/12))

        test <- FinalPvalue(Info.d = myBound$planned$Info.d,
                            Info.i = myBound$planned$Info.i,
                            ck = myBound$planned$ck,
                            ck.unrestricted = myBound$planned$ck,
                            lk = myBound$planned$lk,
                            uk = myBound$planned$uk,
                            kMax = myBound$kMax,
                            estimate = myBound$planned$ck[myBound$kMax]/sqrt(myBound$planned$Info.max),
                            reason.interim = c("","",""),
                            method = iMethod,
                            bindingFutility = TRUE,
                            cNotBelowFixedc = cNotBelowFixedc,
                            continuity.correction = TRUE)

        expect_equal(as.double(test), 0.025, tol = 1e-3)

        ## for non-binding futility rule
        myBound <- suppressWarnings(CalcBoundaries(kMax = 3,
                                                   alpha = 0.025, 
                                                   beta = 0.1, 
                                                   InfoR.i = c(3.5,6.75)/12,
                                                   rho_alpha = 1.345,
                                                   rho_beta = 1.345,
                                                   method = iMethod, 
                                                   cNotBelowFixedc = cNotBelowFixedc,
                                                   bindingFutility = FALSE,
                                                   delta = 1,
                                                   InfoR.d = c(5.5,8.75,12)/12))

            
        test <- FinalPvalue(Info.d = myBound$planned$Info.d,
                            Info.i = myBound$planned$Info.i,
                            ck = myBound$planned$ck,
                            ck.unrestricted = myBound$planned$ck,
                            lk = myBound$planned$lk,
                            uk = myBound$planned$uk,
                            kMax = myBound$kMax,
                            estimate = myBound$planned$ck[myBound$kMax]/sqrt(myBound$planned$Info.max),
                            reason.interim = c("","",""),
                            method = iMethod,
                            bindingFutility = FALSE,
                            cNotBelowFixedc = cNotBelowFixedc,
                            continuity.correction = TRUE)

        expect_equal(as.double(test), 0.025, tol = 1e-3)
    }

})

## * P-value over the space
test_that("Check ordering and continuity of p-values (1 interim analysis)",{

    ## ** binding and no fix C
    calcP_2stage <- function(z, k){
        if(k==1){
            FinalPvalue(Info.d = 2.75046695782793,
                        Info.i = 2.5589317424577,
                        ck = 1.7253424172699, ck.unrestricted = 1.7253424172699,
                        lk = 1.11127641794951,
                        uk = 2.2459290489148,
                        kMax = 2, 
                        delta = 0,
                        estimate = z / sqrt(2.75046695782793),
                        reason.interim = c(""),
                        method = 1,
                        bindingFutility = TRUE, 
                        cNotBelowFixedc = FALSE,
                        continuity.correction = TRUE)
        }else if(k==2){
            FinalPvalue(Info.d = c(2.75046695782793, 4.30612228613852),
                        Info.i = c(2.5589317424577),
                        ck = c(1.7253424172699,2.01974814127831), ck.unrestricted = c(1.7253424172699, 2.01974814127831),
                        lk = c(1.11127641794951),
                        uk = c(2.2459290489148),
                        kMax = 2,
                        delta = 0, 
                        estimate = z / sqrt(4.30612228613852),
                        reason.interim = c("",""),
                        method = 1,
                        bindingFutility = TRUE, 
                        cNotBelowFixedc = FALSE,
                        continuity.correction = TRUE)
        }
    }

    ## type 1 spent at interim: 1-pnorm(2.2459290489148)

    vecP <- c(fut.1_1 = calcP_2stage(z = -3, k = 1),
              fut.1_2 = calcP_2stage(z = 1, k = 1),
              fut.1_B = calcP_2stage(z = 1.725, k = 1),
              fut.2_1 = calcP_2stage(z = -5, k = 2),
              fut.2_2 = calcP_2stage(z = 1, k = 2),
              fut.2_3 = calcP_2stage(z = 2, k = 2),
              eff.2_B = calcP_2stage(z = 2.01974814127831, k = 2),
              eff.2_1 = calcP_2stage(z = 3, k = 2),
              eff.2_2 = calcP_2stage(z = 10, k = 2),
              eff.1_B = calcP_2stage(z = 1.72534241727, k = 1),
              eff.1_1 = calcP_2stage(z = 1.73, k = 1),
              eff.1_2 = calcP_2stage(z = 5, k = 1)
              )
    expect_true(all(round(diff(vecP),5)<=0))

    ## ** binding and fix C
    calcP_2stage <- function(z, k, continuity.correction = TRUE){
        if(k==1){
            FinalPvalue(Info.d = 12.58797814,
                        Info.i = 10.60271846,
                        ck = 1.959964, ck.unrestricted = 1.49772992,
                        lk = 0.24335814,
                        uk = 2.5458844,
                        kMax = 2, 
                        delta = 0,
                        estimate = z / sqrt(12.58797814),
                        reason.interim = c(""),
                        method = 1,
                        bindingFutility = TRUE, 
                        cNotBelowFixedc = TRUE,
                        continuity.correction = continuity.correction)
        }else if(k==2){
            FinalPvalue(Info.d = c(12.58797814, 20.74866926),
                        Info.i = c(10.60271846),
                        ck = c(1.959964, 1.9961649), ck.unrestricted = c(1.49772992, 1.9961649),
                        lk = c(0.24335814),
                        uk = c(2.5458844),
                        kMax = 2,
                        delta = 0, 
                        estimate = z / sqrt(20.74866926),
                        reason.interim = c("",""),
                        method = 1,
                        bindingFutility = TRUE, 
                        cNotBelowFixedc = TRUE,
                        continuity.correction = continuity.correction)
        }
    }

     ## type 1 spent at interim: 1-pnorm(2.5458844) = 0.005450064
    
    vecP <- c(fut.1_1 = calcP_2stage(z = -3, k = 1),
              fut.1_2 = calcP_2stage(z = 1, k = 1),
              fut.1_Bl = calcP_2stage(z = 1.49772992, k = 1),
              fut.1_Bu = calcP_2stage(z = 1.95, k = 1),
              fut.2_1 = calcP_2stage(z = -3, k = 2),
              fut.2_2 = calcP_2stage(z = 1, k = 2),
              fut.2_3 = calcP_2stage(z = 1.99, k = 2),
              eff.2_B = calcP_2stage(z = 1.9961650, k = 2),
              eff.2_1 = calcP_2stage(z = 3, k = 2),
              eff.2_2 = calcP_2stage(z = 10, k = 2),
              eff.1_B = calcP_2stage(z = 1.96, k = 1),
              eff.1_1 = calcP_2stage(z = 2, k = 1),
              eff.1_2 = calcP_2stage(z = 5, k = 1)
              )

    expect_true(all(round(diff(vecP),5)<=0))

})

##----------------------------------------------------------------------
### test-pvalueCI.R ends here

