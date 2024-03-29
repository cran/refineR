% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/algoInvHist.R
\name{testParam}
\alias{testParam}
\title{Helper function to find optimal parameters lambda, mu and sigma.}
\usage{
testParam(
  lambdaVec,
  bestParam,
  Data,
  HistData,
  startValues,
  NIter,
  alpha = 0.01,
  alphaMcb = 0.1
)
}
\arguments{
\item{lambdaVec}{(numeric) transformation parameter for inverse Box-Cox transformation}

\item{bestParam}{(numeric) vector containing best guess for lambda, mu, sigma, P, cost}

\item{Data}{(numeric) values specifying percentiles or data points comprising pathological
and non-pathological values}

\item{HistData}{(list) with histogram data}

\item{startValues}{(list) with start search regions for mu and sigma}

\item{NIter}{(integer) specifying the number of iterations for optimized grid-search}

\item{alpha}{(numeric) specifying the confidence region used for selection of histogram bins in cost calculation}

\item{alphaMcb}{(numeric) specifying the confidence level defining the maximal allowed counts below the asymmetric confidence region}
}
\value{
(numeric) vector with best parameters for lambda, mu, sigma, P, cost.
}
\description{
Helper function to find optimal parameters lambda, mu and sigma.
}
\author{
Tatjana Ammer \email{tatjana.ammer@roche.com}
}
