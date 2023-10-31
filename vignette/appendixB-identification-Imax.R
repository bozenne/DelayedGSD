### appendixB-identification-Imax.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 12 2022 (09:12) 
## Version: 
## Last-Updated: jul 20 2022 (10:34) 
##           By: Brice Ozenne
##     Update #: 41
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(truncnorm)
library(numDeriv)
library(DelayedGSD)
library(mvtnorm)

## * Setting
## CalcBoundaries(kMax = 2)

l1 <- 0.3788015
u1 <- 2.4977055
u2 <- 2.00025
l2 <- 2.00025

n1 <- 50
n2 <- 50
w1 <- sqrt(n1/(n1+n2))
w2 <- sqrt(n2/(n1+n2))

Sigma <- matrix(c(1,w1,w1,1),2,2)
Imax <- 3.64016

delta <- 1.5
mu <- c(delta*sqrt(Imax*n1/(n1+n2)), delta*sqrt(Imax))
  

alpha <- 0.02500
alpha1 <- 0.00625
alpha2 <- alpha - alpha1
beta <- 0.20
beta1 <- 0.05
beta2 <- beta - beta1

Imax <- 3.64016

## * Lemma B.1

## ** Value
GS.B1 <- pmvnorm(lower = c(l1,u2),
                 upper = c(u1,Inf),
                 mean = rep(0,2),
                 sigma= Sigma)
GS.B1
## [1] 0.01875021

integrate(f = function(x){
    pnorm((u2 - w1*x)/w2,lower.tail = FALSE)*dnorm(x)
}, lower = l1, upper = u1, rel.tol = 1e-8)$value
## [1] 0.01875021

(pnorm(u1)-pnorm(l1))-integrate(f = function(x){pnorm((u2 - w1*x)/w2)*dnorm(x)}, lower = l1, upper = u1, rel.tol = 1e-8)$value
## [1] 0.01875021

(pnorm(u1)-pnorm(l1))*(1-integrate(f = function(x){pnorm((u2 - w1*x)/w2)*dtruncnorm(x, a = l1, b = u1)}, lower = -Inf, upper = Inf, rel.tol = 1e-8)$value)
## [1] 0.01875021

## FUN.I <- function(x,y){
##     (w2*x+w1*y<=u2)*dnorm(x)*dtruncnorm(y, a = l1, b = u1)
## }
## (pnorm(u1)-pnorm(l1))*(1-integrate(f = function(x){sapply(x, FUN = function(iX){integrate(f = FUN.I, x = iX, lower = l1, upper = u1, rel.tol = 1e-7)$value})}, lower = -10, upper = 10, rel.tol = 1e-7)$value)

## library(pracma)
## (pnorm(u1)-pnorm(l1))*(1-integral2(fun = FUN.I, xmin = -10, xmax = +10, ymin = l1, ymax = u1, reltol = 1e-8, maxlist = 1e4)$Q)

library(cubature)
FUN.I2 <- function(x){
    (w2*x[1]+w1*x[2]<=u2)*dnorm(x[1])*dtruncnorm(x[2], a = l1, b = u1)
}
res <- (pnorm(u1)-pnorm(l1))*(1-hcubature(f = FUN.I2, lowerLimit = c(-10,l1), upperLimit = c(10,u1), tol = 1e-6)$integral)
res
## [1] 0.01875063

## ** Derivative (du2/dl1)
findU <- function(L){uniroot(f = function(X){
    pmvnorm(lower = c(L,X),
            upper = c(u1,Inf),
            mean = rep(0,2),
            sigma= Sigma)-alpha2
}, lower = 0, upper = 3)$root}
u2.sol <- findU(l1)
u2.sol
## [1] 2.000239
du2.sol <- jacobian(findU, l1)[1,1]
du2.sol
## [1] -0.05680709

pnorm(u1) - pnorm(l1) - integrate(f = function(z){pnorm((u2.sol - w1*z)/w2)*dnorm(z)}, lower = l1, upper = u1, rel.tol = 1e-8)$value - alpha2
- dnorm(l1) + pnorm((u2.sol - w1*l1)/w2)*dnorm(l1)  - integrate(f = function(z){(du2.sol/w2)*dnorm((u2.sol - w1*z)/w2)*dnorm(z)}, lower = l1, upper = u1, rel.tol = 1e-8)$value
- w2 * dnorm(l1) * (1 - pnorm((u2.sol - w1*l1)/w2)) / integrate(f = function(z){dnorm((u2.sol - w1*z)/w2)*dnorm(z)}, lower = l1, upper = u1, rel.tol = 1e-8)$value - du2.sol

- w2 * dnorm(l1) * (1 - pnorm((u2.sol - w1*l1)/w2)) / integrate(f = function(z){dnorm((u2.sol - w1*z)/w2)*dnorm(z)}, lower = l1, upper = u1, rel.tol = 1e-8)$value
## [1] -0.05677705

## * Lemma B.2

## ** Value
GS.B2 <- pmvnorm(lower = c(l1,-Inf),
                 upper = c(u1,l2),
                 mean = mu,
                 sigma = Sigma)
GS.B2
## [1] 0.1499911

FUN.Ibis <- function(x){
    (l1<=x[1])*(x[1]<=u1)*(x[2]<=l2)*dmvnorm(cbind(x[1],x[2]),
                                             mean = mu,
                                             sigma = Sigma)
}
res <- cubature::hcubature(f = FUN.Ibis, lowerLimit = c(l1,-10), upperLimit = c(u1,10))$integral
res
## [1] 0.1499911

FUN.Ibis <- function(x){
    (x[2]<=l2)*dmvnorm(cbind(x[1],x[2]),
                       mean = mu,
                       sigma = Sigma)
}
res <- cubature::hcubature(f = FUN.Ibis, lowerLimit = c(l1,-10), upperLimit = c(u1,10))$integral
res
## [1] 0.1499911

integrate(f = function(x){
    pnorm((l2-w1*x)/w2, mean = (mu[2]-w1*mu[1])/w2)*dnorm(x, mean = mu[1])
}, lower = l1, upper = u1, rel.tol = 1e-8)$value
## [1] 0.1499911

integrate(f = function(x){
    pnorm((l2-w1*x)/w2, mean = delta*sqrt(Imax*n2/(n1+n2)))*dnorm(x, mean = delta*sqrt(Imax*n1/(n1+n2)))
}, lower = l1, upper = u1, rel.tol = 1e-8)$value
## [1] 0.1499911

integrate(f = function(x){
    pnorm((l2-w1*x)/w2 - delta*sqrt(Imax*n2/(n1+n2)))*dnorm(x - delta*sqrt(Imax*n1/(n1+n2)))
}, lower = l1, upper = u1, rel.tol = 1e-8)$value
## [1] 0.1499911

integrate(f = function(x){
    pnorm((l2-w1*x - w1*delta*sqrt(Imax*n1/(n1+n2)))/w2 - delta*sqrt(Imax*n2/(n1+n2)))*dnorm(x)
}, lower = l1 - delta*sqrt(Imax*n1/(n1+n2)), upper = u1 - delta*sqrt(Imax*n1/(n1+n2)), rel.tol = 1e-8)$value
## [1] 0.1499911

deltastar1 <- delta*sqrt(Imax*n1/(n1+n2))
deltastar2 <- delta*sqrt(Imax*n1/(n1+n2))
l1star <- l1 - deltastar1
u1star <- u1 - deltastar1
integrate(f = function(x){
    pnorm((l2 - w1*(x + deltastar1))/w2 - deltastar2)*dnorm(x)
}, lower = l1star, upper = u1star, rel.tol = 1e-8)$value
## [1] 0.1499911

(pnorm(u1star) - pnorm(l1star)) * integrate(f = function(x){
    pnorm((l2 - w1*(x + deltastar1))/w2 - deltastar2)*dtruncnorm(x, a = l1star, b = u1star)
}, lower = l1star, upper = u1star, rel.tol = 1e-8)$value
## [1] 0.1499911


## ** Derivative (du2/dl1)
delta <- 1.5

findL <- function(Info){uniroot(f = function(X){ ## Info <- Imax
    muInfo <- c(delta*sqrt(Info*n1/(n1+n2)), delta*sqrt(Info))
    L1 <- qnorm(beta1, mean = muInfo[1])
    pmvnorm(lower = c(L1,-Inf),
            upper = c(u1,X),
            mean = muInfo,
            sigma= Sigma)-beta2
}, lower = 0, upper = 3)$root}

l2.loc <- findL(Imax)
l2.loc
## [1] 2.000316
dl2.loc <- jacobian(findL, Imax)[1,1]
dl2.loc
## [1] 0.4134894

## dl1/dImax
jacobian(function(x){qnorm(beta1, mean = delta*sqrt(x*n1/(n1+n2)))}, Imax)
delta*sqrt(n1/(4*Imax*(n1+n2)))



beta2 - integrate(f = function(x){
    pnorm((l2 - w1*(x + delta*sqrt(n1*Imax/(n1+n2)))) / w2 - delta*sqrt(n2*Imax/(n1+n2)))*dnorm(x)
}, lower = l1 - delta*sqrt(n1*Imax/(n1+n2)), upper = u1 - delta*sqrt(n1*Imax/(n1+n2)), rel.tol = 1e-8)$value
## [1] 8.935686e-06

integrandL <- function(x){
    pnorm((l2 - w1*(x + delta*sqrt(n1*Imax/(n1+n2)))) / w2 - delta*sqrt(n2*Imax/(n1+n2)))*dnorm(x)
}
dintegrandL <- function(x){
    dnorm((l2 - w1*(x + delta*sqrt(n1*Imax/(n1+n2)))) / w2 - delta*sqrt(n2*Imax/(n1+n2)))*dnorm(x)
}
beta2 - integrate(f = integrandL, lower = l1 - delta*sqrt(n1*Imax/(n1+n2)), upper = u1 - delta*sqrt(n1*Imax/(n1+n2)), rel.tol = 1e-8)$value
## [1] 8.935686e-06

## jacobian(function(x){integrate(f = integrandL, lower = qnorm(beta1, mean = delta*sqrt(n1*x/(n1+n2))) - delta*sqrt(n1*x/(n1+n2)), upper = u1 - delta*sqrt(n1*x/(n1+n2)), rel.tol = 1e-8)$value}, Imax)
## -0.004486283
term1 <- - delta*sqrt(n1/(4*Imax*(n1+n2))) * integrandL(u1 - delta*sqrt(n1*Imax/(n1+n2)))
term2 <- - (delta*sqrt(n1/(4*Imax*(n1+n2))) - delta*sqrt(n1/(4*Imax*(n1+n2)))) * integrandL(l1 - delta*sqrt(n1*Imax/(n1+n2)))

jacobian(function(y){integrate(f = function(x){
    pnorm((l2 - w1*(x + delta*sqrt(n1*y/(n1+n2)))) / w2 - delta*sqrt(n2*y/(n1+n2)))*dnorm(x)
},
lower = l1 - delta*sqrt(n1*Imax/(n1+n2)), upper = u1 - delta*sqrt(n1*Imax/(n1+n2)), rel.tol = 1e-8)$value}, Imax)
##           [,1]
## [1,] -0.093665

termI <- integrate(f = dintegrandL, lower = l1 - delta*sqrt(n1*Imax/(n1+n2)), upper = u1 - delta*sqrt(n1*Imax/(n1+n2)), rel.tol = 1e-8)$value
term3 <- (- w1*delta*sqrt(n1/(4*Imax*(n1+n2)))/w2 - delta*sqrt(n2/(4*Imax*(n1+n2)))) * termI
term4 <- dl2.loc/w2 * termI

term1 + term2 + term3 + term4
## [1] 0.0003724049
term1 + term3 + term4
## [1] 0.0003724049

## w2*(term1+term3)/termI

w2*delta*sqrt(n1)/(2*termI*sqrt((n1+n2)*Imax)) * pnorm((l2 - w1*u1) / w2 - delta*sqrt(n2*Imax/(n1+n2)))*dnorm(u1 - delta*sqrt(n1*Imax/(n1+n2))) + delta * (w1*sqrt(n1)+w2*sqrt(n2)) / (2*sqrt((n1+n2)*Imax))
## [1] 0.4134894

##----------------------------------------------------------------------
### appendixB-identification-Imax.R ends here
