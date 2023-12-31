% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.delayedGSD.R
\name{coef.delayedGSD}
\alias{coef.delayedGSD}
\title{Treatment Effect Estimates for a Group Sequential Design with Delayed Endpoints}
\usage{
\method{coef}{delayedGSD}(
  object,
  type = "effect",
  planned = NULL,
  predicted = TRUE,
  k = NULL,
  type.k = NULL,
  method = NULL,
  ...
)
}
\arguments{
\item{object}{object of class \code{delayedGSD}, typically output from \code{\link{CalcBoundaries}}.}

\item{type}{[character] Shoudl the estimate effect (\code{"effect"}), boundaries (\code{"boundaries"}), information  (\code{"information"}), or decision  (\code{"deicsion"})
be output. The estimate is only displayed for the current stage while the other are displayed for all stages.}

\item{planned}{[logical] Should the planned or estimated boundaries/information be output.}

\item{predicted}{[logical] Should the predicted information/boundaries at decision based on interim data be output (when relevant).}

\item{k}{[integer] stage relative to which the estimates should be output.}

\item{type.k}{[character] type of stage relative to which the estimates should be output: \code{"interim"}, \code{"decision"}, or \code{"final"}.}

\item{method}{[character] type of estimate to output: can be \code{"ML"} or  \code{"corrected ML"}, the latter accounting for the group sequential design.}

\item{...}{not used, for compatibility with the generic method.}
}
\description{
Extract estimate relative to the treatment effect at a specific stage of a group sequential design with delayed endpoints.
By default extract value for the latest stage that has been performed.
}
