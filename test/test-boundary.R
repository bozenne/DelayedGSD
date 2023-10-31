### test-boundary.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  7 2022 (14:48) 
## Version: 
## Last-Updated: mar 22 2023 (14:48) 
##           By: Brice Ozenne
##     Update #: 29
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(testthat)
library(gsDesign)
set.seed(10)

context("Computing bounderies. \n")

NonBindingHJ <- DelayedGSD:::NonBindingHJ

## * slide 106 from CJ_DSBS-v5.pdf
test_that("Boundary calculation with non-binding futility following Jennison method",{


    bCJ <- NonBindingHJ(rho_alpha=1.345,
                        rho_beta=1.345,
                        alpha=0.025,
                        beta=0.1,
                        Kmax=3,
                        Info.max=12,
                        InfoR.i=c(3.5,6.75,12)/12,
                        InfoR.d=c(5.5,8.75)/12,
                        delta=1,  
                        abseps = 1e-06, 
                        direction="smaller",
                        sided=1
                        )

    expect_equal(round(bCJ$boundaries$l.k,3), c(-0.409,0.664,2.069))
    expect_equal(round(bCJ$boundaries$u.k,3), c(2.437,2.244,2.069))
    expect_equal(round(bCJ$boundaries$c.k,3), c(1.960,1.960,NA))
    

    bCJ2 <- Method3(rho_alpha=1.345,
                    rho_beta=1.345,
                    alpha=0.025,
                    beta=0.1,
                    Kmax=3,
                    binding=FALSE,
                    Info.max=12,
                    InfoR.i=c(3.5,6.75)/12,
                    InfoR.d=c(5.5,8.75,12)/12,
                    delta=1,  
                    abseps = 1e-06, 
                    alternative="greater"
                    )


    expect_equal(round(bCJ2$boundaries$lk,3), c(-0.409,0.664,2.069))
    expect_equal(round(bCJ2$boundaries$uk,3), c(2.437,2.244,2.069))
    expect_equal(round(bCJ2$boundaries$ck,3), c(1.960,1.960,2.069))


})

## * gsDesign
test_that("Compare boundaries to gsDesign",{

    ## binding futility bounds (1 interim)
    b1 <- gsDesign(k=2,alpha=0.025,beta=0.2,test.type=3,timing=c(0.6,1),sfu=sfPower,sfl=sfPower,sfupar=2,sflpar=2)
    b12 <- Method1(Kmax=2,delta=1.5,alpha=0.025,beta=0.2,InfoR.i=c(0.6),InfoR.d=c(0.65,1),binding=TRUE)
 
    expect_equal(b1$lower$bound, b12$boundaries[,"lk"], tol = 1e-3)
    expect_equal(b1$upper$bound, b12$boundaries[,"uk"], tol = 1e-3)

    ## binding futility bounds (2 interim)
    b2 <- gsDesign(k=3,alpha=0.025,beta=0.2,test.type=3,timing=c(0.5,0.75,1),sfu=sfPower,sfl=sfPower,sfupar=2,sflpar=2)
    b22 <- Method1(Kmax=3,delta=1.5,alpha=0.025,beta=0.2,InfoR.i=c(0.5,0.75),InfoR.d=c(0.65,0.85,1),binding=T)
 
    expect_equal(b2$lower$bound, b22$boundaries[,"lk"], tol = 1e-3)
    expect_equal(b2$upper$bound, b22$boundaries[,"uk"], tol = 1e-3)
    all.equal(b2$n.I[3],b22$coef, tol = 1e-3)

    ## non-binding futility bounds (1 interim)
    b1nb <- gsDesign(k=2,alpha=0.025,beta=0.2,test.type=4,timing=c(0.6,1),sfu=sfPower,sfl=sfPower,sfupar=2,sflpar=2)
    b12nb <- Method1(Kmax=2,delta=1.5,alpha=0.025,beta=0.2,InfoR.i=c(0.6),InfoR.d=c(0.65,1),binding=FALSE)

    expect_equal(b1nb$lower$bound, b12nb$boundaries[,"lk"], tol = 1e-3)
    expect_equal(b1nb$upper$bound, b12nb$boundaries[,"uk"], tol = 1e-3)
})

## * Stability test
test_that("Check boundaries against previous version",{

    ## binding futility bound
    test1 <- CalcBoundaries(kMax=3,
                            alpha=0.025,  
                            beta=0.1,  
                            InfoR.i=c(3.5,6.75)/12,
                            rho_alpha=1.345,
                            rho_beta=1.345,
                            method=1, 
                            cNotBelowFixedc=FALSE,
                            bindingFutility=TRUE,
                            delta=1,
                            InfoR.d=c(5.5,8.75,12)/12)

    test2 <- suppressWarnings(CalcBoundaries(kMax=3,
                                             alpha=0.025,  
                                             beta=0.1,  
                                             InfoR.i=c(3.5,6.75)/12,
                                             rho_alpha=1.345,
                                             rho_beta=1.345,
                                             method=2, 
                                             cNotBelowFixedc=FALSE,
                                             bindingFutility=TRUE,
                                             delta=1,
                                             InfoR.d=c(5.5,8.75,12)/12))

    test3 <- CalcBoundaries(kMax=3,
                            alpha=0.025,  
                            beta=0.1,  
                            InfoR.i=c(3.5,6.75)/12,
                            rho_alpha=1.345,
                            rho_beta=1.345,
                            method=3, 
                            cNotBelowFixedc=TRUE,
                            bindingFutility=TRUE,
                            delta=1,
                            InfoR.d=c(5.5,8.75,12)/12)

    expect_equal(test1$planned$lk, c(-0.22161265, 0.77399593), tol = 1e-6)
    expect_equal(test1$planned$uk, c(2.59231338, 2.39170973), tol = 1e-6)
    expect_equal(test1$planned$ck, c(1.4057564, 1.72220096, 2.05757005), tol = 1e-6)

    expect_equal(test2$planned$lk, c(-0.21343134, 0.7978885), tol = 1e-6)
    expect_equal(test2$planned$uk, c(2.59231338, 2.39169096), tol = 1e-6)
    expect_equal(test2$planned$ck, c(1.41033082, 1.73451626, 2.05495019), tol = 1e-6)

    expect_equal(test3$planned$lk, c(-0.42681035, 0.63924557), tol = 1e-6)
    expect_equal(test3$planned$uk, c(2.43747723, 2.24378001), tol = 1e-6)
    expect_equal(test3$planned$ck, c(1.95996398, 1.95996398, 2.03661657), tol = 1e-6)


    ## non-binding futility bound and ck>=1.96
    test1 <- CalcBoundaries(kMax=3,
                            alpha=0.025,  
                            beta=0.1,  
                            InfoR.i=c(3.5,6.75)/12,
                            rho_alpha=1.345,
                            rho_beta=1.345,
                            method=1, 
                            cNotBelowFixedc=TRUE,
                            bindingFutility=FALSE,
                            delta=1,
                            InfoR.d=c(5.5,8.75,12)/12)

    test2 <- CalcBoundaries(kMax=3,
                            alpha=0.025,  
                            beta=0.1,  
                            InfoR.i=c(3.5,6.75)/12,
                            rho_alpha=1.345,
                            rho_beta=1.345,
                            method=2, 
                            cNotBelowFixedc=TRUE,
                            bindingFutility=FALSE,
                            delta=1,
                            InfoR.d=c(5.5,8.75,12)/12)

    test3 <- CalcBoundaries(kMax=3,
                            alpha=0.025,  
                            beta=0.1,  
                            InfoR.i=c(3.5,6.75)/12,
                            rho_alpha=1.345,
                            rho_beta=1.345,
                            method=3, 
                            cNotBelowFixedc=TRUE,
                            bindingFutility=FALSE,
                            delta=1,
                            InfoR.d=c(5.5,8.75,12)/12)

    expect_equal(test1$planned$lk, c(-0.199562, 0.80465825), tol = 1e-6)
    expect_equal(test1$planned$uk, c(2.59231338, 2.39218946), tol = 1e-6)
    expect_equal(test1$planned$ck, c(1.95996398, 1.95996398, 2.10213671), tol = 1e-6)

    expect_equal(test2$planned$lk, c(-0.28880108, 0.76708909), tol = 1e-6)
    expect_equal(test2$planned$uk, c(2.59231338, 2.39218947), tol = 1e-6)
    expect_equal(test2$planned$ck, c(1.95996398, 1.95996398, 2.10214493), tol = 1e-6)

    expect_equal(test3$planned$lk, c(-0.40887808, 0.66365936), tol = 1e-6)
    expect_equal(test3$planned$uk, c(2.43747723, 2.24415768), tol = 1e-6)
    expect_equal(test3$planned$ck, c(1.95996398, 1.95996398, 2.0685517), tol = 1e-6)
})

##----------------------------------------------------------------------
### test-boundary.R ends here
