% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.RWDRI}
\alias{plot.RWDRI}
\title{Standard plot method for objects of class 'RWDRI'}
\usage{
\method{plot}{RWDRI}(
  x,
  Scale = c("original", "transformed", "zScore"),
  RIperc = c(0.025, 0.975),
  Nhist = 60,
  showCI = TRUE,
  showPathol = FALSE,
  scalePathol = TRUE,
  showBSModels = FALSE,
  showValue = TRUE,
  CIprop = 0.95,
  pointEst = c("fullDataEst", "medianBS"),
  xlim = NULL,
  ylim = NULL,
  xlab = NULL,
  ylab = NULL,
  title = NULL,
  ...
)
}
\arguments{
\item{x}{(object) of class 'RWDRI'}

\item{Scale}{(character) specifying if percentiles are calculated on the original scale ("or") or the transformed scale ("tr") or the z-Score scale ("z")}

\item{RIperc}{(numeric) value specifying the percentiles, which define the reference interval (default c(0.025, 0.975))}

\item{Nhist}{(integer) number of bins in the histogram (derived automatically if not set)}

\item{showCI}{(logical) specifying if the confidence intervals are shown}

\item{showPathol}{(logical) specifying if the estimated pathological distribution shall be shown}

\item{scalePathol}{(logical) specifying if the estimated pathological distribution shall be weighted with the ration of pathol/non-pathol}

\item{showBSModels}{(logical) specifying if the estimated bootstrapping models shall be shown}

\item{showValue}{(logical) specifying if the exact value of the estimated reference intervals shall be shown above the plot}

\item{CIprop}{(numeric) value specifying the central region for estimation of confidence intervals}

\item{pointEst}{(character) specifying the point estimate determination: (1) using the full dataset ("fullDataEst"),
(2) calculating the median model from the bootstrap samples ("medianBS"), (2) works only if NBootstrap > 0}

\item{xlim}{(numeric) vector specifying the limits in x-direction}

\item{ylim}{(numeric) vector specifying the limits in y-direction}

\item{xlab}{(character) specifying the x-axis label}

\item{ylab}{(character) specifying the y-axis label}

\item{title}{(character) specifying plot title}

\item{...}{additional arguments passed forward to other functions}
}
\value{
The applied plot limits in x-direction (xlim) are returned.
}
\description{
Standard plot method for objects of class 'RWDRI'
}
\author{
Christopher Rank \email{christopher.rank@roche.com}, Tatjana Ammer \email{tatjana.ammer@roche.com}
}
