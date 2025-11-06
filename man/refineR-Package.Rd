\docType{package} 
\name{refineR-package} 
\alias{refineR} 
\alias{refineR-package}
\title{refineR: Reference Interval Estimation using Real-World Data (RWD)} 
\description{ This package includes the implementation of the refineR algorithm (Ammer et al., 2021) which is an indirect method for the estimation of 
reference intervals using Real-World Data (RWD). It takes routine measurements of diagnostic tests, containing pathological and non-pathological samples 
as input and uses sophisticated statistical methods to derive a model describing the distribution of the non-pathological samples. 
This distribution can then be used to derive reference intervals. Version 2.0.0 added functions to check whether a reference interval is valid for the population from which the RWD originates.
Main function of this package is \code{\link{findRI}} that takes an input data set and tries to find a model that best explains the non-pathological distribution.
Furthermore, the package offers functions for printing \code{\link{print.RWDRI}} and plotting \code{\link{plot.RWDRI}} the results of the algorithm operating on S3-objects of class 'RWDRI'.
Verification is performed with the \code{\link{verifyRI}} function.
}
\details{   
 \tabular{ll}{     
 Package: \tab refineR\cr
 Type: \tab Package \cr     
 Version: \tab 2.0.0 \cr     
 Date: \tab 2025-11-05 \cr
 License: \tab GPL (>=3) \cr     
 LazyLoad: \tab yes \cr   
 } 
} 
\author{
	Matthias Beck \email{matthias.beck.mb1@roche.com},
	Tatjana Ammer \email{tatjana.ammer@roche.com},
	Christopher M Rank \email{christopher.rank@roche.com},
	Andre Schuetzenmeister \email{andre.schuetzenmeister@roche.com}
}
\references{
Ammer, T., Schuetzenmeister, A., Prokosch, H.U., Rauh, M., Rank, C.M., Zierk, J. (2021). refineR: A Novel Algorithm for Reference Interval Estimation from Real-World Data. Sci Rep 11, 16023. \doi{10.1038/s41598-021-95301-2}
Ammer, T., Schützenmeister, A., Rank, C. M., Doyle, K. (2023). Estimation of Reference Intervals from Routine Data Using the refineR Algorithm-A Practical Guide. (2023). J Appl Lab Med 8 (1), 84-91. \doi{10.1093/jalm/jfac101}
Beck, M., Dufey, F., Ammer, T., Schützenmeister, A., Zierk, J., Rank, C.M., Rauh, M. (2025). VeRUS: verification of reference intervals based on the uncertainty of sampling.  Clin Chem Lab Med. \doi{10.1515/cclm-2025-0728}.
} 
\keyword{package}
