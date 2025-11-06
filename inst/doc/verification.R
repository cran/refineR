## ----global_options, echo=FALSE, eval=TRUE------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(fig.width = 7, fig.height = 5, fig.align = "center", echo = TRUE, eval = TRUE, warning = FALSE, message = FALSE)
# increasing the width of the stdout-stream
options(width = 200)

## ----loading_example, echo=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(refineR)
fit <- findRI(Data = testcase1)
print(fit, uncertaintyRegion = "uncertaintyMargin")

## ----plot_example, echo=TRUE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# The uncertainty margins are shown in the plot
plot(fit, uncertaintyRegion = "uncertaintyMargin")

## ----verification_1, echo=TRUE------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Suppose the currently used reference interval is [9.1, 31.5]
verifyRI(RIdata = fit, RIcand = c(9.1, 31.5))

## ----one_sided, echo=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
verifyRI(RIdata = fit, RIcand = 67, RIperc = 0.99)

## ----show_possible_inputs_of_verifyRI, fig.show='hold'------------------------------------------------------------------------------------------------------------------------------------------------
# compare RWDRI object with a numeric reference interval
verifyRI(RIdata = fit, RIcand = c(9.1, 33.5),
 title = "Comparison of RWDRI with Numeric RI",
 printResults = FALSE)
# compare two numeric RIs
verifyRI(RIdata = c(4, 26), RIcand = c(1, 29),
 title = "Comparison of two Numeric RIs",
 printResults = FALSE)

# compare two RWDRI objects
# custom RWDRI object describing the reference distribution
# Mu and Sigma define the mean and standard deviation of a normal distribution.
# Then the inverse Box-Cox transformation is applied to the data with the power parameter Lambda.
# The Shift parameter is used to shift the distribution to the desired location.
custom_RWDRI <- list(Mu = 20, Sigma = 5, Lambda = 0.9, Shift = 0)
class(custom_RWDRI) <- "RWDRI"
verifyRI(RIdata = fit, RIcand = custom_RWDRI,
 title = "Comparison of two RWDRI objects",
 printResults = FALSE)

## ----similarity, echo=TRUE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Suppose the currently used reference interval is [9.1, 55]
getRISimilarity(RIdata = fit, RIcand = c(9.1, 55))

## ----output_verifyRI, echo=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# We deactivate the printing of the results and plotting
verification_result <- verifyRI(RIdata = fit, RIcand = c(9.1, 33), printResults = FALSE, generatePlot = FALSE)
verification_tab <- verification_result$RIVerificationTab
# The verification table also contains the values of the power parameter (Lambda) of each estimated model.
# The Lambda values of the models of the candidate distribution and the local distribution are the same when and RWDRI object is compared to a numeric candidate RI.
print(paste0("All Lambdas are equal: ", identical(verification_tab$RICandLambda, verification_tab$RIdataLambda)))
# For this example, the Lambda values are irrelevant, as they are the same.
cols_of_interest <- colnames(verification_tab)
cols_of_interest <- cols_of_interest[!cols_of_interest %in% c("RICandLambda", "RIdataLambda")]

knitr::kable(verification_tab[, cols_of_interest])

## ----adjust_uncertainty_margins, fig.show='hold'------------------------------------------------------------------------------------------------------------------------------------------------------

verifyRI(RIdata = custom_RWDRI, RIcand = c(9.1, 31.5), title = "n=120, UMprop=0.9")
# increasing n may be desired to get a more strict verification
verifyRI(RIdata = custom_RWDRI, RIcand = c(9.1, 31.5), n = 1000, UMprop = 0.95, title = "n=1000, UMprop=0.95")

## ----adjust_uncertainty_margins_2, echo=TRUE----------------------------------------------------------------------------------------------------------------------------------------------------------
plot(fit, uncertaintyRegion = "uncertaintyMargin", n = 1000, UMprop = 0.95, asymmetryCorr = TRUE)
getRI(fit, n = 1000, UMprop = 0.95, asymmetryCorr = TRUE)
print(fit, uncertaintyRegion = "uncertaintyMargin", n = 1000, UMprop = 0.95, asymmetryCorr = TRUE)

## ----adjust_plot_1, echo=TRUE, fig.show='hold'--------------------------------------------------------------------------------------------------------------------------------------------------------
scaleOptions <- c("original", "splitXAxis", "transformed")
labels <- c("default", "Laboratory", "Timepoint")
for(i in 1:3) {

if(i == 1){
  candLabel <- dataLabel <- NULL
}
if(i == 2){
  candLabel <- "Lab A"
  dataLabel <- "Lab B"
}
if(i == 3){
  candLabel <- "Timepoint A"
  dataLabel <- "Timepoint B"
}

  verifyRI(
    RIdata = fit, RIcand = c(9.1, 33), printResults = FALSE,
    xlab = paste("x-axis scale:", scaleOptions[i]),
    title = paste("Verification plot with", scaleOptions[i], "scale"),
    Scale = scaleOptions[i],
    candLabel = candLabel,
    dataLabel = dataLabel
  )
}

## ----similarity_2, echo=TRUE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Suppose the currently used reference interval is [9.1, 55]
getRISimilarity(RIdata = fit, RIcand = c(9.1, 55), UMprop = 0.95, asymmetryCorr = FALSE, Overlap = "OverlapPointEst")

## ----equivalence_limits, echo=TRUE--------------------------------------------------------------------------------------------------------------------------------------------------------------------
verifyRI(
  RIdata = fit, RIcand = c(9.1, 33), printResults = TRUE,
  marginType = "EL",
  title = "Verification with Equivalence Limits"
)

