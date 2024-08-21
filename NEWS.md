# refineR version 1.6.2

The new version includes the following updates and improvements: 

* Restriction of the search region of lambda to the range 0.0 to 1.0 (before it was 0.0 to 1.5)

* Slight modification of model estimation for small sample sizes (N < 1000)

* Performance improvements, especially for data sets with small sample sizes (N < 1000)

* RIbench benchmark score: 0.300, Failure rate: 0.0, Implausible Results: 0.069% 


# refineR version 1.6.1

The new version includes the following updates and improvements: 

* Revised vignette and literature information

* Slight modification of search region for the fraction of non-pathological samples
 
* RIbench benchmark score: 0.307, Failure rate: 0.0, Implausible Results: 0.122%  


# refineR version 1.6.0

The new version includes the following updates and improvements: 

* Using the truncated normal distribution now as basic model assumption 

* Additional regularization term in cost function to control the fit at the tails of the distribution

* The median bootstrapping model (pointEst="medianBS") now corresponds to one parametric model selected from all bootstrapping models

* Some modifications in plot function, e.g. adjusted default `xlim`
 
* Slight modification of print function
 
* RIbench benchmark score: 0.307, Failure rate: 0.0, Implausible Results: 0.122%  


# refineR version 1.5.1

The new version includes the following updates and improvements: 

* Adjusted computation of the `roundingBase` to now also handle data not rounded to  power of 10 (e.g. data that was rounded prior to unit conversion)

* Vignette ('refineR_package') demonstrating the main functions of the package

* RIbench benchmark score: 0.307, Failure rate: 0.0, Implausible Results: 0.625%  


# refineR version 1.5.0 

The new version includes the following updates and improvements: 

* More fine-grained search region for lambda (`lambdaVec`)

* Option to use the two-parameter (modified) Box-Cox transformation
  (`findRI(Data = Data, model = "modBoxCox")`)
  
* Update of calculation of costs: new factor to account for small deviations from the assumption of a unimodal distribution of non-pathological samples
	   
* Adapted definition of region of test results that characterizes the
  non-pathological distribution 
  
* Improved performance for skewed distributions 

* RIbench benchmark score: 0.307, Failure rate: 0.0, Implausible Results: 0.625%  
Detailed description can be found in   Ammer, T., Ammer, T.,Schuetzenmeister, A., Prokosch, HU.,  Zierk, J., Rank, C.M., Rauh, M. RIbench: A Proposed Benchmark for the Standardized Evaluation of Indirect Methods for Reference Interval Estimation. Clinical Chemistry (2022) 
  

# refineR version 1.0.0 

* Initial version put on CRAN 

* Detailed description of the method can be found in Ammer, T.,
  Schuetzenmeister, A., Prokosch, HU., Rauh, M., Rank, C.M., Zierk, J. refineR:
  A Novel Algorithm  for Reference Interval Estimation from Real-World Data. Sci
  Rep 11, 16023 (2021).