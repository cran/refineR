\docType{package} 
\name{refineR-package} 
\alias{refineR} 
\alias{refineR-package}
\title{refineR: Reference Interval Estimation using Real-World Data (RWD)} 
\description{ This package includes the implementation of the refineR algorithm (Ammer et al., 2021) which is an indirect method for the estimation of 
reference intervals using Real-World Data (RWD). It takes routine measurements of diagnostic tests, containing pathological and non-pathological samples 
as input and uses sophisticated statistical methods to derive a model describing the distribution of the non-pathological samples. 
This distribution can then be used to derive reference intervals. Main function of this package is \code{\link{findRI}} that takes an input data set 
and tries to find a model that best explains the non-pathological distribution. Furthermore, the package offers functions for printing \code{\link{print.RWDRI}} and plotting
\code{\link{plot.RWDRI}} the results of the algorithm operating on S3-objects of class 'RWDRI'.
}
\details{   
 \tabular{ll}{     
 Package: \tab refineR\cr
 Type: \tab Package \cr     
 Version: \tab 1.5.1 \cr     
 Date: \tab 2022-07-19 \cr
 License: \tab GPL (>=3) \cr     
 LazyLoad: \tab yes \cr   
 } 
} 
\author{
	Tatjana Ammer \email{tatjana.ammer@roche.com}, 
    Christopher M Rank \email{christopher.rank@roche.com},
    Andre Schuetzenmeister \email{andre.schuetzenmeister@roche.com}
} 
\references{
	Ammer, T., Schuetzenmeister, A., Prokosch, HU., Rauh, M., Rank, C.M., Zierk, J. refineR: A Novel Algorithm for Reference Interval Estimation from Real-World Data. Sci Rep 11, 16023 (2021).
} 
\keyword{package}
