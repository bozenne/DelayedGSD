% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Method2.R
\name{Method2}
\alias{Method2}
\title{Calculate boundaries for a group sequential design using Method 2}
\usage{
Method2(
  rho_alpha = 2,
  rho_beta = 2,
  alpha = 0.025,
  beta = 0.2,
  Kmax,
  Info.max = NULL,
  InfoR.i = NULL,
  InfoR.d = NULL,
  delta = 0,
  abseps = 1e-06,
  alternative = "greater",
  binding = TRUE,
  Trace = FALSE,
  nWhileMax = 30,
  toldiff = 1e-05,
  tolcoef = 1e-04,
  mycoefMax = 1.2,
  mycoefL = 1,
  myseed = 2902,
  cMin = -Inf
)
}
\arguments{
\item{rho_alpha}{rho parameter of the rho-family spending functions (Kim-DeMets) for alpha}

\item{rho_beta}{rho parameter of the rho-family spending functions (Kim-DeMets) for beta}

\item{alpha}{type I error}

\item{beta}{type II error}

\item{Info.max}{maximum information needed for given beta (type II-error), delta (expected difference), alpha (type I-error) and Kmax. It can be given if it is known. Otherwise it is computed from the  values given for alpha, beta, delta and Kmax.}

\item{InfoR.i}{Expected or observed (wherever possible) information rates at the interim analyses 1:(Kmax-1)}

\item{InfoR.d}{Expected or observed information rates at all potential decision analyses and final analyses 1:Kmax}

\item{delta}{expected effect under the alternative (should be on the scale of the test statistc for which If and Info.max relate to one over the variance, e.g. delta=expected log(Hazard ratio))}

\item{abseps}{tolerance for precision when finding roots or computing integrals}

\item{alternative}{a character string specifying the alternative hypothesis, \code{"greater"} or \code{"less"}.
H0 \eqn{\theta=0} vs H1 \eqn{theta<0} (\code{"less"}) or theta > 0 (\code{"greater"}).
Note that in Jennison and Turnbull's book chapter (2013) they consider \code{"greater"}.}

\item{binding}{whether we assume binding futility boundaries}

\item{Trace}{Used only if Info.max=NULL. Whether to print informations to follow the progression of the (root finding) algorithm to compute Info.max (from  alpha, beta, delta and Kmax).}

\item{nWhileMax}{Used only if Info.max=NULL. Maximum number of steps in the (root finding) algorithm to compute Info.max (from  alpha, beta, delta and Kmax)}

\item{toldiff}{Used only if Info.max=NULL. Maximum tolerated difference between lower and upper bounds at anaylis Kmax (which souhld be zero), in the root finding algorithm, to find the value of Info.max}

\item{tolcoef}{Used only if Info.max=NULL. Maximum tolerated difference before stopping the search (in the root finding algorithm), between two successive values for the multiplier coeficient 'coef' such that Info.max=coef*If  (some values for coef are given in Table 7.6 page 164 Jennison's book. The value of "If" (which stands for Information for Fixed design) corresponds to the information we need if Kmax=1)}

\item{mycoefMax}{Used only if Info.max=NULL. Upper limit of the interval of values in which we search for the multiplier coeficient 'coef' such that Info.max=coef*If (in the root finding algorithm).}

\item{mycoefL}{Used only if Info.max=NULL. Lower limit of the interval (see mycoefMax)}

\item{myseed}{seed for producing reproducible results. Because we call functions which are based on Monte-Carlo compuation (pmvnorm)}

\item{cMin}{minimun possible value c for the decision analysis, typically that for a fixed sample test (H & J page 10)}

\item{kMax}{max number of analyses (including final)}
}
\description{
Calculate boundaries for a group sequential design with delayed endpoints based planned information using an error spending approach.
}
\examples{

Example to check that code matches
b1 <- CalcBoundaries(kMax=2,  #max number of analyses (including final)
                    alpha=0.025,  #type I error
                    beta=0.2,  #type II error
                    InfoR.i=c(0.6,1),  #planned information rates
                    rho_alpha=2,  #rho parameter for alpha error spending function
                    rho_beta=2,  #rho parameter for beta error spending function
                    method=1,  #use method 1 or 2 from paper H&J
                    delta=1.5,  #effect that the study is powered for
                    InfoR.d=0.65,
                    bindingFutility=TRUE)

b12 <- Method2(Kmax=2,Info.max=b1$Info.max,delta=1.5,alpha=0.025,InfoR.i=c(0.6),InfoR.d=c(0.65,1))


b1FT <- CalcBoundaries(kMax=2,  #max number of analyses (including final)
                    alpha=0.025,  #type I error
                    beta=0.2,  #type II error
                    InfoR.i=c(0.6,1),  #planned information rates
                    rho_alpha=2,  #rho parameter for alpha error spending function
                    rho_beta=2,  #rho parameter for beta error spending function
                    method=1,  #use method 1 or 2 from paper H&J
                    delta=1.5,  #effect that the study is powered for
                    InfoR.d=0.65,
                    bindingFutility=TRUE)

b12FT <- Method2(Kmax=2,Info.max=b1FT$Info.max,delta=1.5,binding=TRUE,alpha=0.025,InfoR.i=0.6,InfoR.d=c(0.65,1))

b12FTNoImax <- Method2(Kmax=2,Info.max=NULL,delta=1.5,binding=TRUE,alpha=0.025,InfoR.i=0.6,InfoR.d=c(0.65,1))



all.equal(b1$uk, b12$boundaries[,"b.k"])
all.equal(b1$lk, b12$boundaries[,"a.k"])

b14nb <- Method2(Kmax=2,delta=1.5,alpha=0.025,Trace=T,InfoR.i=0.6,InfoR.d=c(0.65,1),binding=F)

b14b <- Method2(Kmax=2,delta=1.5,alpha=0.025,Trace=T,InfoR.i=0.6,InfoR.d=c(0.65,1),binding=T)
}
