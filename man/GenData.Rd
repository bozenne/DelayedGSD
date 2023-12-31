% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenData.R
\name{GenData}
\alias{GenData}
\title{Generate longitudnal (full) data}
\usage{
GenData(
  n = 52,
  N.fw = 2,
  rand.block = c(1, 1, 0, 0),
  allsd = rep(3, N.fw + 1),
  mean0 = rep(0, N.fw + 1),
  delta = rep(0, N.fw + 1),
  ar = 0.86 * 2,
  cor.01.1 = -0.15,
  cor.ij.1 = 0.68,
  cor.0j.1 = -0.27,
  seed = 24082020,
  MissProb = NULL,
  DigitsOutcome = NULL,
  TimeFactor = 1,
  DigitsTime = NULL
)
}
\arguments{
\item{n}{sample size}

\item{N.fw}{number of follow-up measurements (equally spaced, after baseline)}

\item{rand.block}{for block randomization}

\item{allsd}{vector of sd of primary outcome at baseline and end of follow-up (main outcome)}

\item{mean0}{mean outcome at each visit in control group}

\item{delta}{treatment effect (i.e. difference in mean) on primary outcome at each visit}

\item{ar}{accrual rate (average, unit is per time between the equally spaced visits)}

\item{cor.01.1}{correlation between outcome at baseline and at first visit  (main outcome)}

\item{cor.ij.1}{correlation between outcome at two consecutive follow-up measurements  (main outcome)}

\item{cor.0j.1}{correlation between outcome at baseline and at any visit after the first visit  (main outcome)}

\item{seed}{integer for the random seed generator state.}

\item{MissProb}{Missingness probability. Should be an N.fw-dimensional array with in each dimension the proportion of missing and non missing (in that order).
For instance for two follow-up, should be a matrix containing the probability of missing both (1,1),
missing only the first timepoint (1,2)
missing only the second timepoint (2,1)
or have complete data (2,2).}

\item{DigitsOutcome}{Number of digits to round the outcome values (NULL means no rounding)}

\item{TimeFactor}{Multiply the times by a factor (e.g. 14 if time between two follow-up visit should be approx 14 days)}

\item{DigitsTime}{Number of digits to round the times (NULL means no rounding)}
}
\value{
ff
}
\description{
yy
}
\details{
zz
}
\examples{
x <- GenData()
head(x$d,n=20)

}
\author{
Paul Blanche
}
