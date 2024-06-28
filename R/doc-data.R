### doc-data.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 28 2024 (13:58) 
## Version: 
## Last-Updated: jun 28 2024 (14:09) 
##           By: Brice Ozenne
##     Update #: 6
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * abeta
## ** abetaW
#' @title Simulated Data Based on a Phase 3 Trial  
#' @name simLP3trial
#' @rdname data-simLP3trial
#'
#' @description Simulated dataset with properties (e.g. correlation structure between repeated measurements and missingness patterns)
#' similar to a Phase III trial (double-blind, randomized, placebo-controlled, parallel arm) sponsored by H. Lundbeck A/S.
#' The trial evaluate the efficacy of a drug to treat a certain neurodegenerative disease. Patients are randomized 1:1 to an active and a placebo arm.
#' The primary outcome is the change from baseline to Week 12 on a continuous endpoint.
#' This endpoint is measured at Baseline, Week 6 and Week 12.
#' The endpoint of interest is a score between 0 and 20, where higher scores indicate a worse disease status.
#'
#' \itemize{
#' \item \code{id}: study participant.
#' \item \code{Z}: randomization group (0: placebo arm, 1: active arm)
#' \item \code{X1}: endpoint measured at baseline (higher is worse).
#' \item \code{X2}: endpoint measured at week 6 (higher is worse).
#' \item \code{X3}: endpoint measured at week 12 (higher is worse).
#' \item \code{t1}: time (in days) elapsed between the start of the trial and the baseline measurement, i.e., when was the patient recruited w.r.t. the start of the study.
#' \item \code{t2}: time (in days) elapsed between the start of the trial and the week 6 measurement.
#' \item \code{t3}: time (in days) elapsed between the start of the trial and the week 12 measurement.
#' }
#' 
#' @docType data
#' @usage data(simLP3trial)
#' @keywords datasets
NULL

## library(DelayedGSD)
## time_trt <- 12
## Miss11 <- 5/104 # miss both V1 and V2
## Miss12 <- 1/104 # miss V1 and but not V2
## Miss21 <- 6/104 # do not miss V1 and but miss V2
## Miss22 <- 92/104 # miss none
## MyMissProb <- matrix(c(Miss11,Miss12,Miss21,Miss22),ncol=2,nrow=2,byrow=TRUE, # to additionnally remove 1 more because some FASFL=N
##                      dimnames = list(c("V1 missing","V1 not missing"), c("V2 missing","V2 not missing")))
## simLP3trial <- GenData(n=132*2,N.fw=2,allsd=c(2.5,2.1,2.4),
##                mean0 = c(10,0,0),delta = c(0,0.5,1),ar = time_trt/2*2.5,
##                cor.01.1 = -0.15,cor.ij.1 = 0.68,cor.0j.1 = -0.27,
##                MissProb = MyMissProb,TimeFactor = time_trt/2*7,seed=78)$d
## save(simLP3trial, file = "data/simLP3trial.rda")




##----------------------------------------------------------------------
### doc-data.R ends here
