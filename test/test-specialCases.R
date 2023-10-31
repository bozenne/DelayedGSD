### test-specialCases.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  7 2022 (15:00) 
## Version: 
## Last-Updated: mar 24 2023 (18:10) 
##           By: Brice Ozenne
##     Update #: 34
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(testthat)

## * Exceeding maximum information
plannedB <- CalcBoundaries(kMax=3,
                           alpha=0.025,
                           beta=0.2,
                           InfoR.i=c(0.5,0.75),
                           InfoR.d=c(0.55,0.8,1),
                           rho_alpha=2,
                           rho_beta=2,
                           method=1,
                           cNotBelowFixedc=FALSE,
                           bindingFutility=TRUE,
                           delta=0.8)

df.sim <- GenData(n=1000,
                  N.fw=2,
                  rand.block=c(1,1,0,0),
                  allsd=c(2.5, 2.1, 2.4),
                  mean0=c(10, 0, 0),
                  delta=c(0, 0.6, 0.8),
                  ar=3.44,
                  cor.01.1=-0.15,
                  cor.ij.1=0.68,
                  cor.0j.1=-0.27,
                  seed=19519,
                  MissProb=matrix(c(0.04807692, 0.05769231, 0.00961538, 0.88461538), nrow = 2, ncol = 2,
                                  dimnames = list(c("V1 missing", "V1 not missing"),c("V2 missing", "V2 not missing"))),
                  DigitsOutcome=2,
                  TimeFactor=14,
                  DigitsTime=0
                  )

## ** At interim and at decision
test_that("I(1st interim)>Imax, I(1st decison)>Imax",{

    ## interim
    df.interim1 <- df.sim$d[1:500,]
    lmm.interim1 <- analyzeData(df.interim1, ddf = "nlme", getinfo = TRUE, trace = TRUE)
    GSD.interim1 <- update(plannedB, delta = lmm.interim1, trace = FALSE)

    ## final
    lmm.decision1 <- analyzeData(df.sim$d, ddf = "nlme", getinfo = TRUE, trace = TRUE)
    GSD.decision1 <- update(GSD.interim1, delta = lmm.decision1, trace = FALSE)
    ## coef(GSD.decision1, type = "decision")

    ## test
    GS <- data.frame("Stage" = c(1, 2, 3), 
                     "Fbound" = c(-Inf, NA, NA), 
                     "Ebound" = c(Inf, NA, NA), 
                     "statistic.interim" = c(3.432607, NA, NA), 
                     "Cbound" = as.logical(c(NA, NA, NA)),
                     "statistic.decision" = c(NA, NA, NA))
    expect_equivalent(GS, coef(GSD.interim1, type = "boundary"), tol = 1e-5)

})


## ** At interim but not at decision
## ** Not at interim but at decision
## ** At final

## * Decreasing information
## remotes::install_github("paulowhite/DelayedGSD", ref = "a12eb76419")
debug.plannedB <- CalcBoundaries(kMax = 2,  
                                 alpha = 0.025, 
                                 beta = 0.2,  
                                 InfoR.i = c(0.58),  
                                 InfoR.d = c(0.68,1),  
                                 rho_alpha = 2,  
                                 rho_beta = 2,  
                                 method = 1,  
                                 cNotBelowFixedc = TRUE,
                                 bindingFutility= FALSE,
                                 delta = 0.6)
inflationFactor <- debug.plannedB$planned$InflationFactor

df.sim <- GenData(n = 561,  ## only 557 for analysis with method 1 (561 is for method 2)
                  N.fw = 2,
                  rand.block = c(1,1,0,0),
                  allsd = c(2.5,2.1,2.4),
                  mean0 = c(10,0,0),
                  delta = c(0,0.5,1)*0.6,
                  ar = (0.86*2)*2*5,
                  cor.01.1 = -0.15,
                  cor.ij.1 = 0.68,
                  cor.0j.1 = -0.27,
                  seed = 187033903,
                  MissProb = matrix(c(0.04807692, 0.05769231, 0.00961538, 0.88461538), 
                                    nrow = 2, 
                                    ncol = 2, 
                                    dimnames = list(c("V1 missing", "V1 not missing"),c("V2 missing", "V2 not missing")) 
                                    ),
                  DigitsOutcome = 2,
                  TimeFactor = 14,
                  DigitsTime = 0)$d

## ** interim
SelectData(df.sim, t = 247)

debug.lmmI <- analyzeData(SelectData(df.sim, t = 247),
                          ddf = "nlme", data.decision = sum(df.sim$t1 <= 247 + 1.50001*14), getinfo = TRUE, trace = TRUE)
debug.GSDI <- update(debug.plannedB, delta = debug.lmmI, trace = TRUE)

## ** decision
debug.lmmD <- analyzeData(df.sim[which(df.sim$t1 <= 247 + 1.50001*14),],
                          ddf = "nlme", getinfo = TRUE, trace = TRUE)
debug.GSDF <- update(debug.GSDI, delta = debug.lmmD, trace = FALSE)
summary(debug.GSDF)
confint(debug.GSDF)
#### hard to find a good confidence interval


## * Skipped interim analysis
GSD.3t <- CalcBoundaries(kMax=3,  ## max number of analyses (including final)
                         alpha=0.025,  ## type I error
                         beta=0.2,  ## type II error
                         InfoR.i=c(0.5,0.75),  ## planned or observed information rates
                         rho_alpha=2,  ## rho parameter for alpha error spending function
                         rho_beta=2,  ## rho parameter for beta error spending function
                         method=1,  ## use method 1 or 2 from paper H&J
                         delta=1.5,  ## effect that the study is powered for
                         InfoR.d=c(0.55,0.8,1))

set.seed(1322)
GD <- GenData(n=82*2,delta=1.3,ar=5)  ## generate data with all data for in case trial completes
GD.I <- SelectData(GD$d, t = GD$d[50,"t3"])
    
lmmI1.3t <- analyzeData(GD.I, ddf = "nlme", getinfo = TRUE, trace = TRUE)
GSDI1.3t <- update(GSD.3t, delta = lmmI1.3t)

GD.I2 <- SelectData(GD$d, t = GD$d[45,"t3"])
lmmI2.3t <- analyzeData(GD.I2, ddf = "nlme", getinfo = TRUE, trace = TRUE)
GSDI2.3t <- update(GSDI1.3t, delta = lmmI2.3t)

plot(GSDI1.3t)
plot(GSDI2.3t)
plot(GSDI2.3t, cex.estimate = c(2,0.7)) ## plot function seems to work now!

lmmF.3t <- analyzeData(GD$d, ddf = "nlme", getinfo = TRUE, trace = TRUE)
GSDF.3t <- update(GSDI2.3t, delta = lmmF.3t)

confint(GSDF.3t)
plot(GSDF.3t, cex.estimate = c(2,0.7)) 

summary(GSDF.3t)
summary(GSDF.3t, planned = "only")

##----------------------------------------------------------------------
### test-specialCases.R ends here
