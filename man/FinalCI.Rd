% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FinalCI.R
\name{FinalCI}
\alias{FinalCI}
\title{calculate confidence intervals at the end of the study}
\usage{
FinalCI(
  Info.d,
  Info.i,
  ck,
  ck.unrestricted,
  lk,
  uk,
  reason.interim,
  kMax,
  conf.level,
  estimate,
  statistic,
  method,
  bindingFutility,
  cNotBelowFixedc,
  continuity.correction,
  tolerance,
  conclusion,
  FCT.p_value
)
}
\arguments{
\item{Info.d}{Information at all decision analyses up to stage where study was stopped}

\item{Info.i}{Information at all interim analyses up to stage where study was stopped}

\item{ck, ck.unrestricted}{decision boundaries for all decision analyses up to stage where study was stopped (should include final boundary if stopped at final analysis).
ck is possibly with restriction (when cNotBelowFixedc=TRUE) and ck.unrestricted always without.}

\item{lk}{lower bounds up to stage where study was stopped}

\item{uk}{upper bounds up to stage where study was stopped}

\item{reason.interim}{motivation for stopping or continuing at interim. Use to handle special cases (skipped interim because reach Imax, ...)}

\item{kMax}{maximum number of analyses}

\item{conf.level}{confidence level (to get a 100*(1-alpha)\% CI)}

\item{estimate}{naive estimate (e.g. using  ML or REML).}

\item{statistic}{naive test statistic (e.g. using  ML or REML).}

\item{method}{method 1, 2 or 3}

\item{bindingFutility}{[logical]  whether the futility stopping rule is binding.}

\item{cNotBelowFixedc}{[logical] whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)}

\item{continuity.correction}{[logical] whether to add the probability of stopping between ck and ck.uncorrected to ensure continuity of the p-value across stages.
When used the p-value will always be greater than this probability of stopping bettwen ck and ck.uncorrected.}

\item{tolerance}{[numeric] acceptable discrepancy to the objective level when evaluating the confidence intervals and median unbiased estimate.}

\item{FCT.p_value}{[function] function used to compute the p-value.}
}
\description{
calculate confidence intervals at the end of the study
}
