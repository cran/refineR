% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{print.RWDRI}
\alias{print.RWDRI}
\title{Standard print method for objects of class 'RWDRI'}
\usage{
\method{print}{RWDRI}(
  x,
  RIperc = c(0.025, 0.975),
  CIprop = 0.95,
  pointEst = c("fullDataEst", "medianBS"),
  ...
)
}
\arguments{
\item{x}{(object) of class 'RWDRI'}

\item{RIperc}{(numeric) value specifying the percentiles, which define the reference interval}

\item{CIprop}{(numeric) value specifying the central region for estimation of confidence intervals}

\item{pointEst}{(character) specifying the point estimate determination: (1) using the full dataset ("fullDataEst"),
(2) calculating the median model from all bootstrap samples ("medianBS"), (2) works only if NBootstrap > 0}

\item{...}{additional arguments passed forward to other functions.}
}
\value{
No return value. Instead, a summary is printed.
}
\description{
Standard print method for objects of class 'RWDRI'
}
\author{
Christopher Rank \email{christopher.rank@roche.com}
}
