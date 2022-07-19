# refineR version 1.5.1

The new version includes the following updates and improvements: 

* Adjusted detection of rounding to now also handle data not rounded to the power of 10 (e.g. data that was rounded prior to unit conversion)

* Added a vignette ('refineR_package') demonstrating the main functions of the package


# refineR version 1.5.0 

The new version includes the following updates and improvements: 

* More fine-grained search region for lambda (`lambdaVec`)

* Option to use the two-parameter (modified) Box-Cox transformation
  (`findRI(Data = Data, model = "modBoxCox")`)
  
* Update of calculation of costs: new factor to account for small deviations from the assumption of a unimodal distribution of non-pathological samples
	
* Adapted definition of region of test results that characterizes the
  non-pathological distribution 
  
* Improved performance for skewed distributions 



# refineR version 1.0.0 

* Initial version put on CRAN 

* Detailed description of the method can be found in Ammer, T.,
  Schuetzenmeister, A., Prokosch, HU., Rauh, M., Rank, C.M., Zierk, J. refineR:
  A Novel Algorithm  for Reference Interval Estimation from Real-World Data. Sci
  Rep 11, 16023 (2021).