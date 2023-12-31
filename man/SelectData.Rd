% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SelectData.R
\name{SelectData}
\alias{SelectData}
\title{Select available longitudnal data for analysis at a specific follow-up time.}
\usage{
SelectData(d, t)
}
\arguments{
\item{d}{data set generated with GenData function.}

\item{t}{time at which we want the available data.}
}
\value{
input dataset  with less rows and NA where appropriate.
}
\description{
yy
}
\details{
xxx
}
\examples{
x <- GenData(n=50)
head(x$d,n=20)
tail(x$d)
y <- SelectData(x$d,t=4)
head(y)
tail(y)
PlotProgress(x$d,at=4)

#----- select when half of the subjects have one follow-up measuement---
# which depends on accrual rate (ar) 

x <- GenData(n=35)
thear <- 10
thet <- x$d$t2[ceiling(nrow(x$d)/2) + ceiling(thear)]
PlotProgress(x$d,at=thet)
y <- SelectData(x$d,t=thet)
y

}
\author{
Paul Blanche
}
