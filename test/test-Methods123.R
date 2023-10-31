### test-Methods --- 
##----------------------------------------------------------------------
## Author: Corine Baayen
## Created: jan  27 2022 (14:10) 
## Version: 
## Last-Updated: mar 22 2023 (14:45) 
##           By: Brice Ozenne
##     Update #: 3
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

test_that("Boundary calculation with non-binding futility following Jennison Method 3",{
  ## slide 106 from CJ_DSBS-v5.pdf
  
  bCJ <- Method3(rho_alpha=1.345,
                      rho_beta=1.345,
                      alpha=0.025,
                      beta=0.1,
                      Kmax=3,
                      Info.max=12,
                      InfoR.i=c(3.5,6.75)/12,
                      InfoR.d=c(5.5,8.75,12)/12,
                      delta=1,  
                      abseps = 1e-06, 
                      alternative="greater"
  )
  expect_equal(round(bCJ$boundaries$lk,3), c(-0.409,0.664,2.069))
  expect_equal(round(bCJ$boundaries$uk,3), c(2.437,2.244,2.069))
  expect_equal(round(bCJ$boundaries$ck,3), c(1.960,1.960,2.069))
  
})


test_that("Boundary calculation with Method 1",{
  
    ## binding futility bound
    b1 <- gsDesign(k=2,alpha=0.025,beta=0.2,test.type=3,timing=c(0.6,1),sfu=sfPower,sfl=sfPower,sfupar=2,sflpar=2)
    b12 <- Method1(Kmax=2,delta=1.5,alpha=0.025,beta=0.2,InfoR.i=c(0.6),InfoR.d=c(0.65,1),binding=T)
   
  expect_equal(b1$lower$bound, b12$boundaries[,"lk"], tol = 1e-4)
  expect_equal(b1$upper$bound, b12$boundaries[,"uk"], tol = 1e-4)
  
                                        #non-binding futility bound
  b1nb <- gsDesign(k=2,alpha=0.025,beta=0.2,test.type=4,timing=c(0.6,1),sfu=sfPower,sfl=sfPower,sfupar=2,sflpar=2)
  b12nb <- Method1(Kmax=2,delta=1.5,alpha=0.025,beta=0.2,InfoR.i=c(0.6),InfoR.d=c(0.65,1),binding=F)
  
  expect_equal(b1nb$lower$bound, b12nb$boundaries[,"lk"], tol = 1e-4)
  expect_equal(b1nb$upper$bound, b12nb$boundaries[,"uk"], tol = 1e-4)
  
  #3 analyses
  b2 <- gsDesign(k=3,alpha=0.025,beta=0.2,test.type=3,timing=c(0.5,0.75,1),sfu=sfPower,sfl=sfPower,sfupar=2,sflpar=2)
  b22 <- Method1(Kmax=3,delta=1.5,alpha=0.025,beta=0.2,InfoR.i=c(0.5,0.75),InfoR.d=c(0.65,0.85,1),binding=T)
   
  expect_equal(b2$lower$bound, b22$boundaries[,"lk"], tol = 1e-4)
  expect_equal(b2$upper$bound, b22$boundaries[,"uk"], tol = 1e-4)
  expect_equal(b2$n.I[3],b22$coef, tol = 1e-4)
  
})

##----------------------------------------------------------------------
### test-Methods ends here
