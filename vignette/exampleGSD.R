### exampleGSD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Dec 16 2021 (08:56) 
## Version: 
## Last-Updated: feb 23 2022 (16:15) 
##           By: Brice Ozenne
##     Update #: 49
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(DelayedGSD)

## * Design
thear <- 5  ## accrual rate (pt per month)
theDelta.t <- 0.7500001  ## time in months to process data

## planned boundaries
theAlpha <- 0.025
theBeta <- 0.2
theDelta <- 1.5

## * Boundaries
## ** 2 interim
GSD.2t <- CalcBoundaries(kMax=2,  ## max number of analyses (including final)
                         alpha=theAlpha,  ## type I error
                         beta=theBeta,  ## type II error
                         InfoR.i=c(0.5,1),  ## planned or observed information rates
                         rho_alpha=2,  ## rho parameter for alpha error spending function
                         rho_beta=2,  ## rho parameter for beta error spending function
                         method=2,  ## use method 1 or 2 from paper H&J
                         delta=theDelta,  ## effect that the study is powered for
                         InfoR.d=0.55,
                         cNotBelowFixedc=TRUE)

GSD.2t
coef(GSD.2t, type = "boundary")
coef(GSD.2t, type = "information")
summary(GSD.2t)
plot(GSD.2t)

## ** 4 interim
GSD.4t <- CalcBoundaries(kMax=4,  ## max number of analyses (including final)
                         alpha=theAlpha,  ## type I error
                         beta=theBeta,  ## type II error
                         InfoR.i=c(0.25,0.5,0.75,1),  ## planned or observed information rates
                         rho_alpha=2,  ## rho parameter for alpha error spending function
                         rho_beta=2,  ## rho parameter for beta error spending function
                         method=1,  ## use method 1 or 2 from paper H&J
                         delta=theDelta,  ## effect that the study is powered for
                         InfoR.d=c(0.3,0.55,0.8))

GSD.4t
coef(GSD.4t, type = "information")
coef(GSD.4t, type = "boundary")
summary(GSD.4t)
plot(GSD.4t)

## * Data
## ** simulation
set.seed(1322)
GD <- GenData(n=82*2,delta=1.3,ar=thear)  #generate data with all data for in case trial completes

## ** time point at which to do the interim analysis
thet.2t <- GD$d$t3[ceiling(nrow(GD$d)/2) + ceiling(thear*theDelta.t)] ## 2 interim
## mean(x$d$t3<=thet.2t)
dataI.2t <- SelectData(GD$d,t=thet.2t)
dataD.2t <- GD$d[which(GD$d$t1 <= thet.2t + thear*theDelta.t),]
dataF.2t <- GD$d

thet.4t <- c(GD$d$t3[ceiling(nrow(GD$d)/4) + ceiling(thear*theDelta.t)],
             GD$d$t3[ceiling(nrow(GD$d)/2) + ceiling(thear*theDelta.t)], 
             GD$d$t3[ceiling(3*nrow(GD$d)/4) + ceiling(thear*theDelta.t)]) ## 4 interim,''
## c(mean(GD$d$t3<=thet.4t[1]),mean(GD$d$t3<=thet.4t[2]),mean(GD$d$t3<=thet.4t[3]))
dataI1.4t <- SelectData(GD$d,t=thet.4t[1])
dataD1.4t <- GD$d[which(GD$d$t1 <= thet.4t[1] + thear*theDelta.t),]
dataI2.4t <- SelectData(GD$d,t=thet.4t[2])
dataD2.4t <- GD$d[which(GD$d$t1 <= thet.4t[2] + thear*theDelta.t),]
dataI3.4t <- SelectData(GD$d,t=thet.4t[3])
dataD3.4t <- GD$d[which(GD$d$t1 <= thet.4t[3] + thear*theDelta.t),]
dataF.4t <- GD$d

## * First stage
## ** 2 interim
## interim
lmmI.2t <- analyzeData(dataI.2t, ddf = "nlme", getinfo = TRUE, trace = TRUE)
GSDI.2t <- update(GSD.2t, delta = lmmI.2t)

GSDI.2t
confint(GSDI.2t)
coef(GSDI.2t, type = "effect")
coef(GSDI.2t, type = "information", planned = FALSE)
coef(GSDI.2t, type = "information", planned = FALSE, predicted = FALSE)
coef(GSDI.2t, type = "boundary", planned = FALSE)
coef(GSDI.2t, type = "boundary", planned = FALSE, predicted = FALSE)
coef(GSDI.2t, type = "decision")

summary(GSDI.2t)
summary(GSDI.2t, planned = "only")
summary(GSDI.2t, planned = TRUE)
summary(GSDI.2t, planned = FALSE)
summary(GSDI.2t, planned = FALSE, predicted = FALSE)

plot(GSDI.2t, planned = FALSE)
plot(GSDI.2t, planned = TRUE)
plot(GSDI.2t, planned = "only")

## anticipate decision
analyzeData(dataI.2t, ddf = "nlme", data.decision = NROW(dataD.2t), getinfo = TRUE, trace = TRUE)
update(GSD.2t, delta = analyzeData(dataI.2t, ddf = "nlme", data.decision = NROW(dataD.2t), getinfo = TRUE, trace = TRUE))

## decision
lmmD.2t <- analyzeData(dataD.2t, ddf = "nlme", getinfo = TRUE, trace = TRUE)
GSDD.2t <- update(GSDI.2t, delta = lmmD.2t)

GSDD.2t
summary(GSDD.2t)
summary(GSDD.2t, planned = FALSE)
coef(GSDD.2t, type = "information", planned = TRUE)
coef(GSDD.2t, type = "information", planned = FALSE)
confint(GSDD.2t, k = "all")

plot(GSDD.2t, planned = FALSE)
plot(GSDD.2t, planned = TRUE)


## ** 4 interim
lmmI1.4t <- analyzeData(dataI1.4t, ddf = "nlme", getinfo = TRUE, trace = TRUE)
GSDI1.4t <- update(GSD.4t, delta = lmmI1.4t)

confint(GSDI1.4t)
coef(GSDI1.4t, type = "effect")
coef(GSDI1.4t, type = "information", planned = FALSE)
coef(GSDI1.4t, type = "information", planned = FALSE, predicted = FALSE)
coef(GSDI1.4t, type = "boundary", planned = FALSE)
coef(GSDI1.4t, type = "boundary", planned = FALSE, predicted = FALSE)
coef(GSDI1.4t, type = "decision")

GSDI1.4t
summary(GSDI1.4t, planned = "only")
summary(GSDI1.4t, planned = TRUE)
summary(GSDI1.4t, planned = FALSE)
summary(GSDI1.4t, planned = FALSE, predicted = FALSE)

plot(GSDI1.4t, planned = FALSE)
plot(GSDI1.4t, planned = TRUE)
plot(GSDI1.4t, planned = "only")

confint(GSDI1.4t, k = "all")

## add information at decision
lmmD1.4t <- analyzeData(dataD1.4t, ddf = "nlme", getinfo = TRUE, trace = TRUE)
GSDD1.4t <- update(GSDI1.4t, delta = lmmD1.4t, k = 1, type.k = "decision")


## * Second interim
lmmI2.4t <- analyzeData(dataI2.4t, ddf = "nlme", getinfo = TRUE, trace = TRUE)
GSDI2.4t <- update(GSDI1.4t, delta = lmmI2.4t)

confint(GSDI2.4t)
coef(GSDI2.4t, type = "effect")
coef(GSDI2.4t, type = "information", planned = FALSE)
coef(GSDI2.4t, type = "information", planned = FALSE, predicted = FALSE)
coef(GSDI2.4t, type = "boundary", planned = FALSE)
coef(GSDI2.4t, type = "boundary", planned = FALSE, predicted = FALSE)
coef(GSDI2.4t, type = "decision")

GSDI2.4t
summary(GSDI2.4t, planned = "only")
summary(GSDI2.4t, planned = TRUE)
summary(GSDI2.4t, planned = FALSE)
summary(GSDI2.4t, planned = FALSE, predicted = FALSE)

plot(GSDI2.4t, planned = FALSE)
plot(GSDI2.4t, planned = TRUE)
plot(GSDI2.4t, planned = "only")

confint(GSDI2.4t, k = "all")

## decision
lmmD2.4t <- analyzeData(dataD2.4t, ddf = "nlme", getinfo = TRUE, trace = TRUE)
GSDD2.4t <- update(GSDI2.4t, delta = lmmD2.4t)

plot(GSDD2.4t)

##----------------------------------------------------------------------
### exampleGSD.R ends here
