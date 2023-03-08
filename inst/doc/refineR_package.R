## ----global_options, echo=FALSE, eval=TRUE------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(fig.width=7, fig.height=5, fig.align='center', echo=TRUE, eval=TRUE, warning=FALSE, message=FALSE)
			  
# increasing the width of the stdout-stream
options(width=200)

## ----load_testcase4, echo=TRUE------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# load refineR package and load data
library(refineR)
head(testcase4)

## ----run_refineR_default, echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# run refineR estimation and print resulting RWDRI object
fit <- findRI(Data = testcase4)
print(fit)

## ----getRI_default, echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# compute reference intervals using the estimated model parameters
getRI(fit)

## ----plot_default, echo=TRUE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# plot the estimated model 
plot(fit)

## ----run_refineR_bootstrap, echo=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------------------
# run refineR estimation with 20 bootstrap iterations
fit.bs <- findRI(Data = testcase4, NBootstrap = 20)
print(fit.bs)

## ----run_refineR_modBoxCox, echo=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------------------
# run refineR estimation with alternative model (two-parameter (modified) Box-Cox transformation)
fit.mbc <- findRI(Data = testcase4, model = "modBoxCox")
print(fit.mbc)

## ----print_refineR_param, echo=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# compute 2.5%, 50% (median), 97.5% percentiles for the estimated model 
getRI(fit, RIperc = c(0.025, 0.5, 0.975))

# print 2.5%, 50% (median), 97.5% percentiles and estimated model parameters
print(fit, RIperc = c(0.025, 0.5, 0.975))

## ----print_refineR_param_bs, echo=TRUE----------------------------------------------------------------------------------------------------------------------------------------------------------------
# compute percentiles for estimated model with bootstrapping using the median as point estimate  
getRI(fit.bs, RIperc = c(0.025, 0.975), pointEst = "medianBS")

# print percentiles for estimated model with bootstrapping using the median as point estimate and estimated model parameters
print(fit.bs, RIperc = c(0.025, 0.975), pointEst = "medianBS")

## ----plot_param, echo=TRUE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# plot estimated model with bootstrapping with adjusted function arguments
plot(fit.bs, RIperc = c(0.025, 0.5, 0.975), pointEst = "medianBS", xlim = c(0,100), xlab = "Concentration [U/L]", 
		title = "Testcase 4")

## ----plot_showPathol, echo=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# plot estimated model with bootstrapping showing the difference between raw input data and estimated model 
# 		(i.e. 'pathological distribution'), wihtout showing the estimated reference limits
plot(fit.bs, showPathol = TRUE, showValue = FALSE, pointEst = "medianBS", 
		title = "Testcase 4 with pathological distribution")


