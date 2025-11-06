#' Verify Reference Intervals
#'
#' @description This function verifies the reference intervals based on the provided data.
#' 
#' @param RIdata		specifying the RI of the local population: (1) (object) of class \code{RWDRI} or (2) (numeric) representation of reference limits
#' @param RIcand		specifying the RI that needs to be verified: (1) (object) of class \code{RWDRI} or (2) (numeric) representation of reference limits
#' @param RIperc		(numeric) value specifying the percentiles, which define the reference interval (default c(0.025, 0.975))
#' @param UMprop		(numeric) specifying the width of the confidence interval approximated by the uncertainty margins (default 0.9)
#' @param pointEst		(character) specifying the point estimate determination: (1) using the full dataset ("fullDataEst"),
#' 						(2) calculating the median model from the bootstrap samples ("medianBS"), (2) works only if NBootstrap > 0
#' @param printResults	(logical) if TRUE, results are printed in the console.
#' @param generatePlot	(logical) if TRUE, a plot is generated.
#' @param Scale			(character) specifying the scale of the plot: (1) "original", (2) "transformed", (3) "splitXAxis"
#' @param xlab			(character) optional x-axis label for the plot
#' @param title			(character) optional title for the plot
#' @param verbose		(logical) specifying if additional warning messages are printed
#' @param ...			(list) specifying non-default parameters for the calculation of the margins
#'
#' @returns 			(list) containing the verification results
#'
#' @author Matthias Beck \email{matthias.beck.mb1@@roche.com}
#'
#' @examples
#'	\dontrun{
#'
#'	# standard usecase: 
#'	# comparison of the RI estimated from the local population with a numeric candidate
#'	# RI from, e.g., literature
#'	RIdata <- findRI(testcase1)
#'	RIcand <- c(10, 20)
#'	verifyRI(RIdata, RIcand, RIperc = c(0.025, 0.975))
#'
#' # Test if two refineR models are equivalent with stricter criteria, i.e., larger n
#'	model1 <- list(Mu = 20, Sigma = 11, Shift = 8, Lambda = 1,
#'								Method = "manual", roundingBase = NA)
#'	model2 <- list(Mu = 3.41, Sigma = 0.504, Shift = 1, Lambda = 0.06,
#'								Method = "manual", roundingBase = NA)
#'	class(model1) <- class(model2) <- "RWDRI"
#'	verifyRI(RIdata = model1, RIcand = model2, UMprop = 0.95, n = 1e4)
#'
#' # verify two numeric RIs
#' verifyRI(RIdata = c(2.4, 28), RIcand = c(3.1, 29.75))
#' }
#' 
#' @export
#'
verifyRI <- function(RIdata, RIcand, RIperc = c(0.025, 0.975), UMprop = 0.9, pointEst = c("fullDataEst", "medianBS"),
						printResults = TRUE, generatePlot = TRUE, Scale = c("original", "transformed", "splitXAxis"),
						xlab = NULL, title = NULL, verbose = TRUE, ...) {

	stopifnot(inherits(RIdata, "RWDRI") || (is.numeric(RIdata) && length(RIdata) >= 2))

	if((is.numeric(RIdata)) && (is.numeric(RIcand) && length(RIcand) == 1)){
		stop("RIdata must be the result of refineR (RWDRI class) if RIcand is a single number")
	}

	if(is.numeric(RIcand)){
		stopifnot(length(RIcand) == length(RIperc))
		stopifnot(length(RIcand) >= 1)
	}

	if(is.numeric(RIdata)){
		stopifnot(length(RIdata) == length(RIperc))
		stopifnot(length(RIdata) >= 1)
	}

	pointEst <- match.arg(pointEst[1], choices = c("fullDataEst", "medianBS"))
	Scale <- match.arg(Scale[1], choices = c("original", "transformed", "splitXAxis"))

	stopifnot(is.logical(verbose))
	stopifnot(is.logical(generatePlot))
	stopifnot(is.null(title) || is.character(title))
	stopifnot(is.numeric(RIperc), all(RIperc >= 0), all(RIperc <= 1))

	if(verbose && length(RIperc) == 1 && is.numeric(RIcand)){
		message("Verification involves interpolation as only one limit of RIcand is provided. Interpret results with caution.")
	}

	# sort the input vectors
	RIpercSorting <- sort(RIperc, index.return = TRUE)
	RIperc <- RIpercSorting$x
	if(is.numeric(RIdata)){
		RIdata <- RIdata[RIpercSorting$ix]
		if(is.unsorted(RIdata)){
			stop("Order of numbers in RIdata does not match the order of percentiles in RIperc")
		}
	}
	if(is.numeric(RIcand)){
		RIcand <- RIcand[RIpercSorting$ix]
		if(is.unsorted(RIcand)){
			stop("Order of numbers in RIcand does not match the order of percentiles in RIperc")
		}
	}

	rm(RIpercSorting)

	# select the marginType
	param <- list(...)
	if("marginType" %in% names(param)){
		marginType <- param$marginType
		marginType <- match.arg(marginType[1], choices = c("VeRUS", "asymptotic", "EL"))
		param$marginType <- NULL
	} else {
		marginType <- "VeRUS"
	}
	param$UMprop <- UMprop

	if(marginType == "asymptotic"){
		# to ensure compatibility with older versions of the function
		marginType <- "VeRUS"
	}

	# get labels for the plot
	candLabel <- param$candLabel
	dataLabel <- param$dataLabel
	candLabel <- if(is.null(candLabel)) "RI Cand" else as.character(candLabel)
	dataLabel <- if(is.null(dataLabel)) "RI Data" else as.character(dataLabel)
	# remove from param so that no warning is thrown in checkInvalidArgs
	param$candLabel <- param$dataLabel <- NULL

	# get the verification arguments
	tmp <- getVerificationArgs(RIdata, RIcand = RIcand, RIperc = RIperc, pointEst = pointEst, verbose = verbose)
	dataTestArgs <- tmp$dataTestArgs
	candTestArgs <- tmp$candTestArgs
	rm(tmp)

	# Calculate margins
	if(marginType == "VeRUS"){
		# Check for and inform about wrong arguments
		checkInvalidArgs(param, c("n", "UMprop", "asymmetryCorr"), marginType, verbose)

		asymmetryCorr <- TRUE
		if(!is.null(param$asymmetryCorr) && is.logical(param$asymmetryCorr) && !param$asymmetryCorr){
			asymmetryCorr <- FALSE
		}

		# Calculate margins
		marginsRIcand <- getRIMargins(RI = candTestArgs$RI, RIperc = candTestArgs$RIperc,
										lambda = candTestArgs$Lambda, shift = candTestArgs$Shift,
										n = ifelse(is.null(param$n), 120, param$n),
										UMprop = ifelse(is.null(param$UMprop), 0.9, param$UMprop),
										asymmetryCorr = asymmetryCorr 
										)

		marginsRIdata <- getRIMargins(RI = dataTestArgs$RI, RIperc = dataTestArgs$RIperc,
										lambda = dataTestArgs$Lambda, shift = dataTestArgs$Shift,
										n = ifelse(is.null(param$n), 120, param$n),
										UMprop = ifelse(is.null(param$UMprop), 0.9, param$UMprop),
										asymmetryCorr = asymmetryCorr
										)

	} else if(marginType == "EL"){

		# Check for and inform about wrong arguments
		checkInvalidArgs(param, c("pCVA.exp", "n", "with.bias", "CIprop"), marginType, verbose)

		# Calculate margins
		marginsRIcand <- getEquivalenceLimits(RI = candTestArgs$RI, RIperc = candTestArgs$RIperc,
		                                      CIprop = ifelse(is.null(param$CIprop), 0.9, param$CIprop),
		                                      pCVA.exp = ifelse(is.null(param$pCVA.exp), 0.5, param$pCVA.exp),
		                                      n = ifelse(is.null(param$n), NULL, param$n),
		                                      with.bias = ifelse(is.logical(param$with.bias), param$with.bias, FALSE)
		                                    )

		marginsRIdata <- getEquivalenceLimits(RI = dataTestArgs$RI, RIperc = dataTestArgs$RIperc,
		                                      CIprop = ifelse(is.null(param$CIprop), 0.9, param$CIprop),
		                                      pCVA.exp = ifelse(is.null(param$pCVA.exp), 0.5, param$pCVA.exp),
		                                      n = ifelse(is.null(param$n), NULL, param$n),
		                                      with.bias = ifelse(is.logical(param$with.bias), param$with.bias, FALSE)
		                                     )
	} else {
		stop("marginType must be one of 'asymptotic', 'PU', 'EL'")
	}

	# Create verification table
	result <- createVerificationTab(RIperc = RIperc, marginsRIdata = marginsRIdata, marginsRIcand = marginsRIcand)

	# generate plot
	if (generatePlot) {

		verificationTab <- result$RIVerificationTab
		marginOverlap <- c()
		for (i in 1:nrow(verificationTab)) {
			if(!is.na(verificationTab$OverlapPointEst[i]) && verificationTab$OverlapPointEst[i]){
				marginOverlap <- c(marginOverlap, "PEOverlap")
			} else if(!is.na(verificationTab$OverlapMargins[i]) && verificationTab$OverlapMargins[i]){
				marginOverlap <- c(marginOverlap, "MarOverlap")
			} else {
				marginOverlap <- c(marginOverlap, "noOverlap")
			}
		}

		marginsRIdataPlot <- marginsRIdata[marginsRIdata$Percentile %in% RIperc, ]
		marginsRIcandPlot <- marginsRIcand[marginsRIcand$Percentile %in% RIperc, ]


		plotRIVerification(mar1 = marginsRIcandPlot, mar2 = marginsRIdataPlot, marginOverlap=marginOverlap, Scale=Scale, title = title, xlab = xlab, candLabel=candLabel, dataLabel=dataLabel)
		rm(marginsRIdataPlot, marginsRIcandPlot, marginOverlap)
	}

	if(marginType == "VeRUS"){
		result$RIVerificationTab$RICandLambda <- candTestArgs$Lambda
		result$RIVerificationTab$RIDataLambda <- dataTestArgs$Lambda
	} else {
		result$RIVerificationTab$RICandLambda <- 0 #Haeckel et al. assume log-normal
		result$RIVerificationTab$RIDataLambda <- 0 #Haeckel et al. assume log-normal
	}

	#print result
	if(printResults) {
		# Print header
		marginTypeStr <- ifelse(marginType == "VeRUS", "VeRUS", "EL")
		headerText <- paste0("Verification of Reference Interval", " [", marginTypeStr, "]")
		cat("\n", headerText, "\n")
		cat(paste(rep("-", nchar(headerText)), collapse = ""), "\n\n")

		printVerificationTable(result$RIVerificationTab, RIperc)
		if(inherits(RIdata, "RWDRI") && "Data" %in% names(RIdata) && is.numeric(RIdata$Data) && length(RIdata$Data) > 0){
			cat("\n")
			cat("Fraction of data within/outside of reference interval:")
			printDataFractionWithinOutsideRI(result$RIVerificationTab, RIperc, RIdata$Data)
		}
	}

	invisible(result)
}

#' Check Invalid Arguments
#'
#' @description This function checks if the arguments in param are valid
#' 
#' @param param			(list)  providing non-default parameters for the calculation of the margins
#' @param validArgs		(character) vector of valid arguments for the calculation of the margins
#' @param marginType	(character) specifying calculation of the margins: using "VeRUS", or Equivalence Limits ("EL")
#' @param verbose		(logical) specifying if additional warning messages are printed
#'
#' @returns warning message if invalid arguments are found
#'
#' @author Matthias Beck \email{matthias.beck.mb1@@roche.com}
#'
checkInvalidArgs <- function(param, validArgs, marginType, verbose) {
	invalidArgs <- setdiff(names(param), validArgs)
		if(verbose && length(invalidArgs) > 0) {
			message(paste0("'", paste(invalidArgs, collapse = "', '"), "' are not valid arguments in param if marginType is '", marginType, "'"))
		}
}


#' Get the correct values for RI, RIperc, Lambda, and Shift for the verification
#'
#' @description This function adapts the list of test arguments based on the margin type and additional test parameters
#'
#' @param RIdata	specifying the RI of the local population: (1) (object) of class RWDRI or (2) (numeric) representation of reference limits
#' @param RIcand	specifying the RI that needs to be verified: (1) (object) of class RWDRI or (2) (numeric) representation of reference limits
#' @param RIperc	(numeric) value specifying the percentiles, which define the reference interval (default c(0.025, 0.975))
#' @param pointEst	(character) specifying the point estimate determination: (1) using the full dataset ("fullDataEst"),
#' 							(2) calculating the median model from the bootstrap samples ("medianBS"), (2) works only if NBootstrap > 0
#' @param verbose  (logical) specifying if additional warning messages are printed
#' 
#' @returns (list) containing a list with the "RI", "Lambda", "Shift" parameter for  RIdata and RIcand each
#' 
#' @author Matthias Beck \email{matthias.beck.mb1@@roche.com}
#'
getVerificationArgs <- function(RIdata, RIcand, RIperc, pointEst = c("fullDataEst", "medianBS"), verbose = TRUE) {
	pointEst <- match.arg(pointEst[1], choices = c("fullDataEst", "medianBS"))

	RIdataIsValidRWDRI <- inherits(RIdata, "RWDRI") && is.numeric(RIdata$Lambda) && is.numeric(RIdata$Shift) && is.numeric(RIdata$Mu) && is.numeric(RIdata$Sigma)
	RIcandIsValidRWDRI <- inherits(RIcand, "RWDRI") && is.numeric(RIcand$Lambda) && is.numeric(RIcand$Shift) && is.numeric(RIcand$Mu) && is.numeric(RIcand$Sigma)

	stopifnot(RIdataIsValidRWDRI || is.numeric(RIdata))
	stopifnot(RIcandIsValidRWDRI || is.numeric(RIcand))
	stopifnot(is.numeric(RIperc))

	# determine if Bootstrap or full data estimate should be used
	BSPerformed <- FALSE
	if(pointEst == "medianBS" && RIdataIsValidRWDRI) {
		BSPerformed  <- (length(RIdata$MuBS) > 0 && length(RIdata$SigmaBS) > 0 && length(RIdata$LambdaBS) > 0 && length(RIdata$ShiftBS) > 0) &&
						(!is.na(RIdata$LambdaMed) && !is.na(RIdata$MuMed) && !is.na(RIdata$SigmaMed) && !is.na(RIdata$ShiftMed))


		if(BSPerformed && RIcandIsValidRWDRI) {
		BSPerformed  <- (length(RIcand$MuBS) > 0 && length(RIcand$SigmaBS) > 0 && length(RIcand$LambdaBS) > 0 && length(RIcand$ShiftBS) > 0) &&
		                (!is.na(RIcand$LambdaMed) && !is.na(RIcand$MuMed) && !is.na(RIcand$SigmaMed) && !is.na(RIcand$ShiftMed))
		}
	}
	usedEst <- ifelse(BSPerformed, "medianBS", "fullDataEst")
	if(verbose && pointEst == "medianBS" && !BSPerformed) {
		warning("Bootstrap estimates are not available. Full data estimates will be used.")
	}

		# select appropriate arguments for the verification
	if(is.numeric(RIcand)) {
		if(length(RIcand) == 1 && !RIdataIsValidRWDRI) {
			stop("If candidate RI is a single number, the data RI must be a valid RWDRI object")
		}

		if(length(RIcand) == 1 && RIdataIsValidRWDRI){
		# special case: RIcand is a single number
			if (RIperc >= 0.5 - 1e-3 && RIperc <= 0.5 + 1e-3) {
				RIperc <- c(0.1, RIperc)
				RIcand <- c(NA, RIcand)
				RIcand[1] <- getRI(RIdata, RIperc = RIperc[1], pointEst = usedEst, UCMargins = FALSE)$PointEst
			} else if(RIperc < 0.5){
				RIperc <- c(RIperc, 1 - RIperc)
				RIcand <- c(RIcand, NA)
				RIcand[2] <- getRI(RIdata, RIperc = RIperc[2], pointEst = usedEst, UCMargins = FALSE)$PointEst
			} else if(RIperc > 0.5){
      	RIperc <- c(1 - RIperc, RIperc)
				RIcand <- c(NA, RIcand)
				RIcand[1] <- getRI(RIdata, RIperc = RIperc[1], pointEst = usedEst, UCMargins = FALSE)$PointEst
			}
			if(RIcand[2] < RIcand[1]){
				stop("Interpolation failed: Model does not fit the data")
			}
		}

		if(RIdataIsValidRWDRI) {
			#Data RWDRI and candidate is numeric
			candTestArgs <- list(RI = RIcand, RIperc = RIperc,
				Lambda = ifelse(usedEst == "medianBS", RIdata$LambdaMed, RIdata$Lambda),
				Shift  = ifelse(usedEst == "medianBS", RIdata$ShiftMed, RIdata$Shift))

  		dataTestArgs <- list(RI = getRI(RIdata, RIperc = RIperc, pointEst = usedEst, UCMargins = FALSE)$PointEst,
				RIperc = RIperc,
  			Lambda = ifelse(usedEst == "medianBS", RIdata$LambdaMed, RIdata$Lambda),
  			Shift  = ifelse(usedEst == "medianBS", RIdata$ShiftMed,  RIdata$Shift))

		} else if (is.numeric(RIdata)) {
			# Both numeric
			candTestArgs <- list(RI = RIcand, RIperc = RIperc, Lambda = 0, Shift  = 0)
			dataTestArgs <- list(RI = RIdata, RIperc = RIperc, Lambda = 0, Shift  = 0)
		}

	} else if (RIcandIsValidRWDRI) {
		# candidate RWDRI
		tmp <- getRIfromRWDRI(RWDRI = RIcand, RIperc = RIperc, pointEst = usedEst)
		RIperc <- tmp$RIperc
		candTestArgs <- list(RI = tmp$RI, RIperc = tmp$RIperc,
			Lambda = ifelse(usedEst == "medianBS", RIcand$LambdaMed, RIcand$Lambda),
			Shift  = ifelse(usedEst == "medianBS", RIcand$ShiftMed,  RIcand$Shift))

		if(RIdataIsValidRWDRI) {
			# Both RWDRI
			tmp <- getRIfromRWDRI(RWDRI = RIdata, RIperc = RIperc, pointEst = usedEst)
			RIperc <- tmp$RIperc
			dataTestArgs <- list(RI = tmp$RI, RIperc = tmp$RIperc,
				Lambda = ifelse(usedEst == "medianBS", RIdata$LambdaMed, RIdata$Lambda),
				Shift  = ifelse(usedEst == "medianBS", RIdata$ShiftMed,  RIdata$Shift))

		} else if (is.numeric(RIdata)) {
			# Data numeric and candidate RWDRI
			stop("Consider switching the inputs of RIdata and RIcand; RIdata should be the result of refineR")
		}
	}

	stopifnot(length(dataTestArgs$RI) == length(candTestArgs$RI) && length(candTestArgs$RI) == length(RIperc))
	#dataTestArgs$RIperc <- candTestArgs$RIperc <- RIperc

	return(list(dataTestArgs = dataTestArgs, candTestArgs = candTestArgs))
}

#' Get Reference Interval from RWDRI
#' @description Helper function to get the reference interval from a RWDRI object. For internal use only.
#' @param RWDRI (RWDRI) The RWDRI object containing the reference interval information.
#' @param RIperc (numeric) The percentiles to interpolate.
#' @param pointEst (character) The point estimate to use ("fullDataEst" or "medianBS").
#' @returns (list) A list containing the interpolated reference interval and percentiles.
getRIfromRWDRI <- function(RWDRI, RIperc, pointEst){
	# The function makes sure that the RI is calculated for at least two percentiles
	# This is necessary so that the UMs can be calculated in the next step with getRIMargins
  stopifnot(inherits(RWDRI, "RWDRI"))
  stopifnot(is.numeric(RIperc), length(RIperc) >= 1, all(RIperc > 0 & RIperc < 1))
  stopifnot(pointEst %in% c("fullDataEst", "medianBS"))

  if(length(RIperc) == 1){
    if(RIperc >= 0.5 - 1e-3 && RIperc <= 0.5 + 1e-3){
      RIperc <- c(0.1, RIperc)

    } 
    else if(RIperc < 0.5){
      RIperc <- c(RIperc, 1 - RIperc)

    } else if(RIperc > 0.5){
      RIperc <- c(1 - RIperc, RIperc)

  	}
	}
  tmp <- getRI(RWDRI, RIperc = RIperc, pointEst = pointEst, UCMargins = FALSE)
  return(list(RI = tmp$PointEst, RIperc = RIperc))
}



#' Create Verification Table
#'
#' @description	Method to summarize the UMS of RIdata and RIcand and checks of the UMs overlap
#'
#' @param RIperc		(numeric) value specifying the percentiles, which define the reference interval (default c(0.025, 0.975))
#' @param marginsRIdata	(data.frame) Margins for the data-derived reference intervals
#' @param marginsRIcand	(data.frame) Margins for the candidate reference intervals
#' 
#' @returns				(list) containing whether the point estimate and margin overlap tests passed, and the verification table
#'
#' @author Matthias Beck \email{matthias.beck.mb1@@roche.com}
#'
createVerificationTab <- function(RIperc, marginsRIdata, marginsRIcand) {
		# create verification table
		RIVerificationTab <- data.frame(
		Percentile       = marginsRIcand$Percentile,
		RICandPointEst   = marginsRIcand$PointEst,
		RICandMarginLow  = marginsRIcand$MarginLow,
		RICandMarginHigh = marginsRIcand$MarginHigh,
		RIDataPointEst   = marginsRIdata$PointEst,
		RIDataMarginLow  = marginsRIdata$MarginLow,
		RIDataMarginHigh = marginsRIdata$MarginHigh,
		OverlapPointEst  = NA,
		OverlapMargins   = NA
	)

	# Adpating the table to the case of a single percentile
	if(length(RIperc) == 1) {
		relevantRows <- which(RIVerificationTab$Percentile == RIperc)
		RIVerificationTab <- RIVerificationTab[relevantRows,]
	}
	# Check if the point estimates and margins overlap
	for(row in seq_len(nrow(RIVerificationTab))) {

		marginCandLim <- c(RIVerificationTab$RICandMarginLow[row], RIVerificationTab$RICandMarginHigh[row])
		marginDataLim <- c(RIVerificationTab$RIDataMarginLow[row], RIVerificationTab$RIDataMarginHigh[row])

		RIcand <- RIVerificationTab$RICandPointEst[row]
		RIdata <- RIVerificationTab$RIDataPointEst[row]

		RIVerificationTab[row, 'OverlapPointEst'] <- all(
			c(RIcand >= marginDataLim[1]  && RIcand <= marginDataLim[2],
				RIdata >= marginCandLim[1]  && RIdata <= marginCandLim[2])
		)
		if(is.na(RIVerificationTab[row, 'OverlapPointEst'])) RIVerificationTab[row, 'OverlapPointEst'] <- FALSE

		RIVerificationTab[row, 'OverlapMargins'] <- any(
			c(marginCandLim[1:2] >= marginDataLim[1] & marginCandLim[1:2] <= marginDataLim[2],
				marginDataLim[1:2] >= marginCandLim[1] & marginDataLim[1:2] <= marginCandLim[2])
		)
		if(is.na(RIVerificationTab[row, 'OverlapMargins'])) RIVerificationTab[row, 'OverlapMargins'] <- FALSE
	}

	testPassedPointEst <- all(RIVerificationTab$OverlapPointEst)
	testPassedMargins  <- all(RIVerificationTab$OverlapMargins)

	return(list(testPassedPointEst = testPassedPointEst, testPassedMargins = testPassedMargins, RIVerificationTab = RIVerificationTab))
}


#' Calculate uncertainty margins for a reference interval using the asymptotic method
#'
#' @param RI 			(numeric) vector of length >= 2 representing the lower and upper limits of the reference interval.
#' @param RIperc		(numeric) value specifying the percentiles, which define the reference interval (default \code{c(0.025, 0.975)})
#' @param UMprop		(numeric) value between 0 and 1 representing the confidence level for the uncertainty margins.
#' @param lambda		(numeric) value representing the power parameter for the Box-Cox transformation.
#' @param shift			(numeric) value representing the shift parameter for the Box-Cox transformation.
#' @param asymmetryCorr	(logical) value if the asymmetry correction shall be applied for extremely skewed distributions
#' @param n				(numeric) value representing the sample size for which the sampling uncertainty shall be taken into account
#' @param verbose		(logical) specifying if additional warning messages are printed
#'
#' @returns				(data.frame) containing the calculated uncertainty margins for the reference interval.
#' 
#' @examples
#' \dontrun{
#'	getRIMargins(RI = c(12, 65), lambda = 0, shift = 0)
#'	getRIMargins(RI = c(12, 65), lambda = 1, shift = 0)
#'  the following examples should return NAs with and without asymmetry correction
#'  # Next examples should return NAs
#'	getRIMargins(RI = c(78.6, 116.7), lambda = 0.30673426, shift = 120.898)
#'	getRIMargins(RI = c(78.6, 116.7), lambda = 0.30673426, shift = 120.898, asymmetryCorr = FALSE)
#' }
#' 
#' @author Christopher Rank \email{christopher.rank@@roche.com}
#'
getRIMargins <- function(RI, RIperc = c(0.025, 0.975), UMprop = 0.9, lambda = 0, shift = 0, asymmetryCorr = FALSE, n = 120, verbose = TRUE) {
	# check input
	stopifnot(is.finite(RI))
	stopifnot(is.finite(RIperc), length(RIperc) > 1, all(RIperc>0 & RIperc<1))
	stopifnot(length(RI) == length(RIperc))
	stopifnot(is.numeric(UMprop), length(UMprop)==1, UMprop>=0, UMprop<=1)
	stopifnot(is.numeric(lambda), length(lambda)==1)
	stopifnot(is.numeric(shift))
	stopifnot(is.logical(asymmetryCorr))
	stopifnot(is.numeric(n), length(n)==1)
	stopifnot(is.logical(verbose))

	# sort RI input vectors
	RIindex <- sort(RIperc, index.return=TRUE)$ix
	RIperc  <- RIperc_filterd <- RIperc[RIindex]
	RI      <- RI[RIindex]

	# restrict values close to zero
	if(any(RI < 1e-6) && verbose) warning("RI values are coerced to be greater then or equal to 1e-6")
	RI		<- pmax(1e-6, RI)

	RITra <- BoxCox(x=RI-shift, lambda=lambda)

	#values that are not finite after Box-Cox transformation are not used for the calculation of the uncertainty margins
	if(any(!is.finite(RITra))){
		if(verbose){
		  warning("Some RI values are not finite after Box-Cox transformation and are excluded from the calculation of the uncertainty margins.")
		}

		RI_tmp <- RI[which(is.finite(RITra))]
		RIperc_filterd <- RIperc[which(is.finite(RITra))]
		if(length(RI_tmp) < 2){
			if(verbose) warning("Uncertainty margins cannot be calculated.")

			tabRIMargins <- data.frame(Percentile=RIperc, PointEst=RI, MarginLow=NA, MarginHigh=NA)
			return(tabRIMargins)
		} else {
			RITra  <- BoxCox(x=RI_tmp-shift, lambda=lambda)
		}
	}

	fit   <- lm(y ~ x, data=data.frame(x=qnorm(p=RIperc_filterd), y=RITra), na.action = na.fail)
	mu    <- as.numeric(coef(fit)[1])
	sigma <- as.numeric(coef(fit)[2])

  bestAsymmetryCorr <- 0
		
	if(asymmetryCorr)
	{
		# calculate central 95% RI
		RI95 <- invBoxCox(x=qnorm(p=c(0.025, 0.975), mean=mu, sd=sigma), lambda=lambda)+shift
		RI95[RI95<1e-3]   <- 1e-3
		RI95[is.na(RI95)] <- 1e-3

		diffRI95 <- RI95[2]-RI95[1]
		if(diffRI95 < 1e-5){
			if(verbose) warning("Uncertainty margins cannot be calculated.")
			tabRIMargins <- data.frame(Percentile=RIperc, PointEst=RI, MarginLow=NA, MarginHigh=NA)
			return(tabRIMargins)
		}

		# calculate initial ratio between width of upper margin and width of RI for a central 95% RI with n=120 and UMprop=0.9
		UcLimit975     <- as.numeric(perc_ci_asymptotic(mu=mu, sigma=sigma, n=120, prob=0.975, conf.level=0.9))
		RIMargin975    <- invBoxCox(x=UcLimit975, lambda=lambda)+shift


		# check if NAs were introduced by inverse Box-Cox transformation
		if(all(is.finite(RIMargin975))) {
			ratioMargin975 <- (RIMargin975[2]-RIMargin975[1])/diffRI95
		} else if(sum(is.finite(RIMargin975)) == 1) {
			RIMargin975[!is.finite(RIMargin975)] <- 1e-6
			ratioMargin975 <- (RIMargin975[2]-RIMargin975[1])/diffRI95
			if(ratioMargin975 <= 0) ratioMargin975 <- Inf
		} else {
			ratioMargin975 <- Inf
		}

		thresholdRatio    <- 0.75

		# Correct asymmetry of uncertainty margins for extremely skewed distributions
		if(ratioMargin975>=thresholdRatio)
		{
			# vector with potential values for asymmetry correction
			asymmetryCorrTemp <- 0.04*diff(range(RI95))*(0:4095)^4/4095^4

			indexLower  <- 0
			indexUpper  <- 4096
			indexSearch <- 1

			for(i in 1:12)
			{
				# set next search index
				if(indexSearch == floor((indexUpper+indexLower)/2)) break

				indexSearch <- (indexUpper+indexLower)/2

				RITra95 <- BoxCox(x=RI95-shift+asymmetryCorrTemp[indexSearch], lambda=lambda)
				mu95    <- 0.5*(RITra95[1] + RITra95[2])
				sigma95 <- (RITra95[2]-mu95)/qnorm(p=0.975)

				# update ratio between width of upper margin and width of RI for a central 95% RI with n=120 and UMprop=0.9
				UcLimit975     <- as.numeric(perc_ci_asymptotic(mu=mu95, sigma=sigma95, n=120, prob=0.975, conf.level=0.9))
				RIMargin975    <- invBoxCox(x=UcLimit975, lambda=lambda)+shift-asymmetryCorrTemp[indexSearch]

				if(all(is.finite(RIMargin975))){
				ratioMargin975 <- (RIMargin975[2]-RIMargin975[1])/diffRI95
				} else {
					# RIMargin975 is NA if the Box-Cox is done with negative values
					# the shift must be increased to avoid negative values
					indexLower 		 <- indexSearch
					next
				} 

				# adapt lower or upper index and make search region smaller
				if(ratioMargin975<thresholdRatio)
				{
					indexUpper 		    <- indexSearch
					bestAsymmetryCorr <- asymmetryCorrTemp[indexSearch]
				} else
				{
					indexLower <- indexSearch
				}
			}
		}

		if(bestAsymmetryCorr == 0){
			#if(verbose) warning("Asymmetry correction did not converge")
		} else {
			# update mu and sigma with asymmetry correction
			RITra <- BoxCox(x=RI-shift+bestAsymmetryCorr, lambda=lambda)	
			fit   <- lm(y ~ x, data=data.frame(x=qnorm(p=RIperc_filterd), y=RITra), na.action = na.fail)
			mu    <- as.numeric(coef(fit)[1])
			sigma <- as.numeric(coef(fit)[2])
		}
	}

	UcLimits  <- as.numeric(perc_ci_asymptotic(mu=mu, sigma=sigma, n=n, prob=RIperc, conf.level=UMprop))
	RIMargins <- invBoxCox(x=UcLimits, lambda=lambda)+shift-bestAsymmetryCorr

	RIMargins[RIMargins<0] <- 0
	RIMargins[!is.finite(RIMargins)] <- NA
	tabRIMargins <- data.frame(Percentile=RIperc, PointEst=RI, MarginLow=RIMargins[(1:(length(RIMargins)/2))*2-1], MarginHigh=RIMargins[(1:(length(RIMargins)/2))*2])

	tabRIMargins <- tabRIMargins[order(tabRIMargins$Percentile), ]
	return(tabRIMargins)
}


#' Function to approximate the sampling uncertainty of quantiles using the asymptotic method
#'
#' @param mu			(numeric) mean of the distribution
#' @param sigma			(numeric) standard deviation of the distribution
#' @param n				(numeric) integer value indicating the sample size
#' @param prob			(numeric) quantile value(s) for which to calculate the confidence interval(s)
#' @param conf.level	(numeric) confidence level for the interval(s)
#' 
#' @returns (matrix) of confidence interval(s) of quantiles of the normal distribution
#' 
#' @author Florian Dufey \email{florian.dufey@@roche.com}
#' 
#' @references Serfling RJ. Approximation theorems of mathematical statistics. NY: John Wiley & Sons; 1980:121 p.
#'
#' @examples 
#' \dontrun{
#'	perc_ci_asymptotic(mu = 0, sigma = 1, n = 120, prob = c(0.025, 0.975), conf.level = 0.9)
#' }
#' 
perc_ci_asymptotic <- function(mu, sigma, n, prob, conf.level) {
	stopifnot(is.numeric(prob) && all(prob >= 0 & prob <= 1))
	stopifnot(is.numeric(conf.level) && conf.level > 0 && conf.level < 1)
	if(length(prob) > 1) {
		tmp <- sapply(prob, perc_ci_asymptotic, mu=mu, sigma=sigma, n=n, conf.level=conf.level)
		colnames(tmp) <- paste0("Q", 100*prob, "%")
		return(tmp)
	}
	alpha<- 1-conf.level
	perc <- qnorm(p=prob, mean=mu, sd=sigma)
	se     <- sigma*sqrt(prob*(1-prob)/n)/dnorm(qnorm(prob))
	o      <- perc+qnorm(p=c(alpha/2, conf.level+alpha/2))*se
	names(o) <- paste0(100*conf.level, "%", c("LCL", "UCL"))
	o
}


#' Calculate equivalence limits
#'
#' @param RI		(numeric) vector of length 2 representing the lower and upper limits of the reference interval
#' @param RIperc	(numeric) value specifying the percentiles, which define the reference interval (default c(0.025, 0.975))
#' @param CIprop	(numeric) specifying the width of the confidence interval used to determine Equivalence Limits (default 0.9)
#' @param pCVA.exp	(numeric) value greater than 0 representing the exponent for the calculation of pCV.A
#' @param with.bias	(logical) value indicating whether to consider bias according to Haekel 2015
#' @param n			(numeric) value representing the sample size for consideration of bias
#' 
#' @returns 		(data.frame) containing the calculated equivalence limits
#' 
#' @author Matthias Beck \email{matthias.beck.mb1@@roche.com}
#' 
#' 
getEquivalenceLimits <- function(RI, RIperc = c(0.025, 0.975), CIprop = 0.9, pCVA.exp = 0.5, with.bias = FALSE, n = NULL) {
# references
# Haeckel, R., Wosniok, W. & Arzideh, F. (2016). Equivalence limits of reference intervals for
# partitioning of population data. Relevant differences of reference limits.
# LaboratoriumsMedizin, 40(3), 199-205. https://doi.org/10.1515/labmed-2016-0002
# 
# Haeckel, R., Wosniok, W., Gurr, E. & Peil, B. (2015). Permissible limits for uncertainty of measurement
# in laboratory medicine. Clinical Chemistry and Laboratory Medicine (CCLM), 53(8), 1161-1171.
# https://doi.org/10.1515/cclm-2014-0874

	if (RI[1] < 0 || RI[2] < 0 || RI[1] >= RI[2]) {
		stop("(permissible_uncertainty) please check limits, only positive values allowed.")
	}

	if (!is.numeric(RIperc) || length(RIperc) != 2 || abs((0.5 - RIperc[1]) - (RIperc[2] - 0.5)) > 1e-10) {
		stop("RIperc must be a numeric vector of length 2 with values that have to be symmetrical with 0.5 in the middel.")
	}

	if (!is.numeric(CIprop) || CIprop < 0 || CIprop > 1) {
		stop("CIprop must be a numeric value between 0 and 1.")
	}

	if (!is.numeric(pCVA.exp) || pCVA.exp < 0) {
		stop("pCVA.exp must be a numeric value greater than 0.")
	}

	if (with.bias && !(is.numeric(n) && (n %% 1 == 0))) {
		stop("Please provide a sample size for consideration of random variation")
	}
	
	s.E.ln <- (log(RI[2]) - log(RI[1])) / (2 * qnorm(RIperc[2])) # EQ 5
	med.ln <- (log(RI[2]) + log(RI[1])) / 2 # EQ 5

	CV.E <- 100 * (sqrt(exp(s.E.ln^2) - 1)) # EQ 8

	pCV.A <- (CV.E - 0.25)^pCVA.exp # EQ 9
	psA.med <- pCV.A * 0.01 * exp(med.ln) # EQ 10

	psA1 <- 0.2 * psA.med + RI[1] * (0.8 * psA.med) / exp(med.ln)
	psA2 <- 0.2 * psA.med + RI[2] * (0.8 * psA.med) / exp(med.ln)

	alpha <- 1 - CIprop

	if (with.bias) { # Haekel 2015

		factUB <- qt((1 - alpha) / 2, n - 1) / sqrt(n) # Haekel 2015  appendix
		factUB <- sqrt(factUB^2 + factUB^2) # Haekel 2015  EQ 9

		psA1 <- psA1 * factUB # Haekel 2015  EQ 10
		psA2 <- psA2 * factUB # Haekel 2015  EQ 10
	}

	factoPsA <- qnorm(CIprop) # Haekel 2015  EQ 13 for with.randVar

	EL <- data.frame(
			Percentile = RIperc, PointEst = RI,
			# pD = c(psA1 * factoPsA, psA2 * factoPsA), # for testing of the fucntion only
			MarginLow  = c(RI[1] - psA1 * factoPsA, RI[2] - psA2 * factoPsA),
			MarginHigh = c(RI[1] + psA1 * factoPsA, RI[2] + psA2 * factoPsA)
		)
	return(EL)
}


#' Plot method for RI verification
#' 
#' @param mar1			(data.frame) with output of function getRIMargins()
#' @param mar2			(data.frame) with output of function getRIMargins()
#' @param marginOverlap	(character) vector length of mar1; specifying if and how the margins overlap; options are: "PEOverlap", "MarsOverlap", "noOverlap"
#' @param Scale			(character) specifying if percentiles are shown on the original scale ("original") or the transformed scale ("transformed") or a split view of the x-axis ("splitXAxis")
#' @param lambda		(numeric) specifying the power parameter (skewness) of the assumed distribution
#' @param xlab			(character) specifying the x-axis label
#' @param title			(character) specifying plot title
#' @param ...				(additional arguments) candLabel (character) specifying the label for the candidate RI; dataLabel (character) specifying the label for the data-derived RI
#' 
#' @return					(NULL) Instead, a plot is generated.
#' 
#' @author Christopher Rank \email{christopher.rank@@roche.com}
#' 
plotRIVerification <- function(mar1, mar2, marginOverlap, Scale=c("original", "transformed", "splitXAxis"), lambda=0, xlab=NULL, title=NULL, ...)
{
	stopifnot(nrow(mar1) == nrow(mar2))
	stopifnot(nrow(mar1) == length(marginOverlap))
	stopifnot(all(mar1$Percentile == mar2$Percentile))

	Scale <- match.arg(Scale[1], choices = c("original", "transformed", "splitXAxis"))

	args <- list(...)
	candLabel <- args$candLabel
	dataLabel <- args$dataLabel
	candLabel <- if(is.null(candLabel)) "RI Cand" else as.character(candLabel)
	dataLabel <- if(is.null(dataLabel)) "RI Data" else as.character(dataLabel)

	pointEstOrig1 <- signif(mar1$PointEst, 3)
	pointEstOrig2 <- signif(mar2$PointEst, 3)

	if(Scale=="transformed")
	{
		mar1$PointEst   <- BoxCox(pmax(1e-6, mar1$PointEst),   lambda=lambda)
		mar1$MarginLow  <- BoxCox(pmax(1e-6, mar1$MarginLow),  lambda=lambda)
		mar1$MarginHigh <- BoxCox(pmax(1e-6, mar1$MarginHigh), lambda=lambda)
		mar2$PointEst   <- BoxCox(pmax(1e-6, mar2$PointEst),   lambda=lambda)
		mar2$MarginLow  <- BoxCox(pmax(1e-6, mar2$MarginLow),  lambda=lambda)
		mar2$MarginHigh <- BoxCox(pmax(1e-6, mar2$MarginHigh), lambda=lambda)

	} else if(Scale=="splitXAxis")
	{
		xlabLabels <- NULL
		xlabPos    <- NULL

		for(i in 1:nrow(mar1))
		{		
			xlabPretty <- pretty(range(c(as.numeric(unlist(mar1[i,-1 ])), as.numeric(unlist(mar2[i,-1 ]))), na.rm = TRUE), n=floor(10/nrow(mar1)))

			xlimSplit <- range(xlabPretty)
			xlimSplit <- xlimSplit + c(-0.05, 0.05)*diff(xlimSplit)

			fit <- lm(y~x, data=data.frame(x=xlimSplit, y=c((i-1), i)*100/nrow(mar1)))

			mar1$PointEst[i]   <- predict(fit, newdata=data.frame(x=mar1$PointEst[i]))
			mar1$MarginLow[i]  <- predict(fit, newdata=data.frame(x=mar1$MarginLow[i]))
			mar1$MarginHigh[i] <- predict(fit, newdata=data.frame(x=mar1$MarginHigh[i]))
			mar2$PointEst[i]   <- predict(fit, newdata=data.frame(x=mar2$PointEst[i]))
			mar2$MarginLow[i]  <- predict(fit, newdata=data.frame(x=mar2$MarginLow[i]))
			mar2$MarginHigh[i] <- predict(fit, newdata=data.frame(x=mar2$MarginHigh[i]))

			xlabLabels <- c(xlabLabels, xlabPretty)
			xlabPos	   <- c(xlabPos, predict(fit, newdata=data.frame(x=xlabPretty)))
		}
	}

	if (is.null(xlab) & Scale=="original" | Scale=="splitXAxis")
		xlab <- "Concentration [Units]"
	else if (is.null(xlab) & Scale=="transformed")
		xlab <- "Transformed Scale"

	if (is.null(title))	
		title <- paste0("Verification of Reference Interval: ", ifelse(all(marginOverlap!="noOverlap"), "TEST PASSED", "TEST NOT PASSED"))

	#bottom, left, top, and right
	par(mar=c(5.1, 3.1, 4.1, 2.1))

	if(Scale == "splitXAxis")
		xlim <- c(0, 100)
	else
		xlim <- range(c(as.numeric(unlist(mar1[,-1 ])), as.numeric(unlist(mar2[,-1 ]))), na.rm = TRUE)

	plot(x=0, y=-1000, xlim=xlim, ylim=c(0.72, 2.4), xlab=xlab, ylab="", xaxt="n", yaxt="n", main=title, cex.main=1.3)

	if(Scale != "splitXAxis")
	{
		xlabPos <- pretty(xlim, n = 10)
		axis(1, at = xlabPos, cex.axis = 0.85)
		addGrid(x=xlabPos, y=rep(NA, times=length(xlabPos)), lty=2)

		rect(xleft=min(mar1$PointEst), ybottom=1.7, xright=max(mar1$PointEst), ytop=1.9, border=NA, col=as.rgb("grey30", alpha=0.65))
		rect(xleft=min(mar2$PointEst), ybottom=1.0, xright=max(mar2$PointEst), ytop=1.2, border=NA, col=as.rgb("royalblue3", alpha=0.75))

	} else
	{
		abline(v=(1:(nrow(mar1)-1))*100/nrow(mar1))
		axis(1, at=xlabPos, labels=xlabLabels, cex.axis = 0.85)
		addGrid(x=xlabPos, y=rep(NA, times=length(xlabPos)), lty=2)

		for(i in (1:(nrow(mar1)-1))*100/nrow(mar1))
		{
			lines(x=i+c(-1.5, 1.5), y=c(0.72, 0.79), lwd=2)
			lines(x=i+c(-1.5, 1.5), y=c(0.75, 0.82), lwd=2)
		}
	}

	for(i in 1:nrow(mar1))
	{
		rect(xleft=mar1$MarginLow[i], ybottom=1.5, xright=mar1$MarginHigh[i], ytop=1.7, border=NA, col=as.rgb("grey30", alpha=0.65/2.5))	
		rect(xleft=mar2$MarginLow[i], ybottom=1.2, xright=mar2$MarginHigh[i], ytop=1.4, border=NA, col=as.rgb("royalblue3", alpha=0.75/2.5))	

		xleft  <- max(mar1$MarginLow[i], mar2$MarginLow[i])
		xright <- min(mar1$MarginHigh[i], mar2$MarginHigh[i]) 

		if(marginOverlap[i] == "PEOverlap")
		{
			rect(xleft=xleft, ybottom=1.4, xright=xright, ytop=1.5, border=NA, col=as.rgb("green2", alpha=0.36))

		} else if(marginOverlap[i] == "MarOverlap")
		{
			rect(xleft=xleft, ybottom=1.4, xright=xright, ytop=1.5, border=NA, col=as.rgb("gold2", alpha=0.40))

		} else if(marginOverlap[i] == "noOverlap")
		{
			rect(xleft=xleft, ybottom=1.4, xright=xright, ytop=1.5, border=as.rgb("red2", alpha=0.5), col=NA, lwd=3)
		}

		lines(x=rep(mar1$PointEst[i], 2), y=c(1.5, 2.1), col="grey30", lty=1, lwd=3)
		lines(x=rep(mar2$PointEst[i], 2), y=c(0.8, 1.4), col="blue2", lty=1, lwd=3)
	}

	abline(h=1.45, lty=2, lwd=2, col="grey50")
	text(x=mar1$PointEst, y=2.15, labels=pointEstOrig1, col="grey30")
	text(x=mar2$PointEst, y=0.75, labels=pointEstOrig2, col="blue2")

	mtext(text=dataLabel, side=2, line=1, at=1.1, cex=1.3, col="blue2")
	mtext(text=candLabel, side=2, line=1, at=1.8, cex=1.3, col="grey30")
	mtext(text="Reference Intervals with Uncertainty Margins", col="black", cex=1.0)

	legendNames <- c("No overlap", "Overlap of margins", "Overlap with PE")  

	if(Scale == "splitXAxis")
		legend("top", border=c(as.rgb("red2", alpha=0.5), "black", "black"), fill=c(as.rgb("white", alpha=0.80), as.rgb("gold2", alpha=0.40), col=as.rgb("green2", alpha=0.36)), legend=legendNames, ncol=3, cex=1, bg="white", x.intersp=0.5)
	else
		legend("top", border=c(as.rgb("red2", alpha=0.5), "black", "black"), fill=c(as.rgb("white", alpha=0.80), as.rgb("gold2", alpha=0.40), col=as.rgb("green2", alpha=0.36)), legend=legendNames, ncol=3, cex=1, bty="n", x.intersp=0.5)

	#bottom, left, top, and right
	par(mar=c(5.1, 4.1, 4.1, 2.1))
}


#' Print Verification Table
#'
#' @description This function prints a formatted verification table for reference intervals
#'
#' @param verificationTab	(data.frame) containing the verification table. Must include the following columns:
#'	- Percentile: The percentile values.
#'	- RIDataPointEst: Point estimates for RI Data.
#'	- RIDataMarginLow: Lower margins for RI Data.
#'	- RIDataMarginHigh: Upper margins for RI Data.
#'	- RITestPointEst: Point estimates for RI Cand.
#'	- RITestMarginLow: Lower margins for RI Cand.
#'	- RITestMarginHigh: Upper margins for RI Cand.
#'	- OverlapPointEst: Logical indicating if point estimates overlap.
#'	- OverlapMargins: Logical indicating if margins overlap.
#'
#' @param RIperc			(numeric) value specifying the percentiles, which define the reference interval (default c(0.025, 0.975))
#'
#' @returns None. Prints the verification table to the console.
#'
#' @examples
#' \dontrun{
#'	df <- data.frame(
#'		Percentile = c(0.025, 0.975),
#'		RICandPointEst = c(12, 65),
#'		RICandMarginLow = c(10.9013, 60.2780),
#'		RICandMarginHigh = c(13.86800, 72.64152),
#'		RIDataPointEst = c(13, 69),
#'		RIDataMarginLow = c(11.60130, 63.77848),
#'		RIDataMarginHigh = c(14.39861, 74.22152),
#'		OverlapPointEst = c(TRUE, TRUE),
#'		OverlapMargins = c(TRUE, TRUE)
#'	)
#'  
#'	RIperc <- c(0.025, 0.975)
#'	printVerificationTable(df, RIperc)
#' }
printVerificationTable <- function(verificationTab, RIperc) {
	# Input validation
	stopifnot(is.data.frame(verificationTab))

	# Filter verification table
	i_row <- sapply(verificationTab$Percentile, function(p) any(abs(p - RIperc) < 1e-10))
	if(sum(i_row) == 0) {
		stop("No matching percentiles found in the verification table.")
	}

	verificationTab <- verificationTab[i_row, ]

	required_cols <- c("Percentile", "RIDataPointEst", "RICandPointEst", 
					"RIDataMarginLow", "RIDataMarginHigh",
					"OverlapPointEst", "OverlapMargins")
	stopifnot(all(required_cols %in% colnames(verificationTab)))

	Ndecimal <- 3

	# Align numeric vectors
	aligned_nums <- list(
		data_points = alignVec(verificationTab$RIDataPointEst),
		data_low    = alignVec(verificationTab$RIDataMarginLow),
		data_high   = alignVec(verificationTab$RIDataMarginHigh),
		cand_points = alignVec(verificationTab$RICandPointEst),
		cand_low    = alignVec(verificationTab$RICandMarginLow),
		cand_high   = alignVec(verificationTab$RICandMarginHigh)
	)

	# Generate formatted rows
	rows <- lapply(seq_len(nrow(verificationTab)), function(i) {
		overlap <- if (is.finite(verificationTab$OverlapPointEst[i]) && verificationTab$OverlapPointEst[i]) {
			"Overlap with PE"
		} else if (is.finite(verificationTab$OverlapMargins[i]) && verificationTab$OverlapMargins[i]) {
			"Overlap of margins"
		} else {
			"No overlap"
		}

		list(
			Perc = sprintf("%4.1f%%", verificationTab$Percentile[i] * 100),
			`RI Data` = sprintf(
				"%s (%s; %s)", 
				aligned_nums$data_points[i],
				aligned_nums$data_low[i],
				aligned_nums$data_high[i]
			),
			`RI Cand` = sprintf(
				"%s (%s; %s)", 
				aligned_nums$cand_points[i],
				aligned_nums$cand_low[i],
				aligned_nums$cand_high[i]
			),
			Overlap = overlap
		)
	})

	# Calculate column widths
	headers <- c("Perc", "RI Data", "RI Cand", "Overlap")
		# Calculate column widths
	stopifnot(length(headers) == length(rows[[1]]))
	col_widths <- numeric()
	for (i in seq_len(length(headers))) {
		col_widths <- c(col_widths, max(nchar(headers[i]), sapply(rows, function(row) nchar(row[[i]]))))
	}

	# Determine test result
	test_result <- if (all(!is.finite(verificationTab$OverlapPointEst))){
		"TEST NOT PASSED"
	} else if (all(verificationTab$OverlapPointEst, na.rm = TRUE)) {
		"TEST PASSED (all point estimates within uncertainty margins)"
	} else if (all(verificationTab$OverlapMargins, na.rm = TRUE)) {
		"TEST PASSED (all uncertainty margins overlap)"
	} else {
		"TEST NOT PASSED"
	}

	# Print Verification Result
	cat("Verification result:", test_result, "\n\n")

	# Print table header
	printTableHeader(headers, col_widths)

	# Print table rows
	printTableRows(lapply(rows, function(r) unlist(r)), col_widths)

	cat("\n")
}

#' Helper function to print Table Header
#'
#' @param headers 		(character) vector of column headers
#' @param col_widths 	(numeric)	vector of column widths
#' 
#' @returns 			(NULL) Prints the header to the console.
#' 
#' @author Matthias Beck \email{matthias.beck.mb1@@roche.com}
#' 
printTableHeader <- function(headers, col_widths) {
	# Create header format string
	header_format <- paste(sprintf("%%-%ds", col_widths), collapse = " | ")
	separator_format <- paste(sprintf("%%%ds", col_widths), collapse = " | ")

	# Print headers
	line <- do.call(sprintf, c(header_format, as.list(headers)))
	cat(line, "\n")

	# Print separator
	seperator <- lapply(col_widths, function(w) paste(rep("-", w), collapse = ""))
	line <- do.call(sprintf, c(separator_format, seperator))
	cat(line, "\n")
}


#' Helper function to print Table Rows
#'
#' @param rows 			(list) of lists containing row data.
#' @param col_widths 	(numeric) vector of column widths.
#'
#' @returns 			(NULL) Prints the rows to the console.
#' 
#' @author Matthias Beck \email{matthias.beck.mb1@@roche.com}
#' 
#' 
printTableRows <- function(rows, col_widths) {
	# Create row format string
	row_format <- paste(sprintf("%%%ds", col_widths), collapse = " | ")
	row_format <- paste(row_format, "\n", sep = "")

	# Print each row
	for (row in rows) {
		line <- do.call(sprintf, 
		as.list(c(row_format, unname(unlist(row))))
		)
		cat(line)
	}
}


#' Print Data Fractions Within and Outside Reference Intervals
#'
#' This function prints a table showing the fractions of data points that are below, within, and above specified reference intervals.
#' The function handles both lower and upper limits, or just one of them, based on the input.
#'
#' @param tab		(dataframe) must contain columns `Percentile`, `RICandPointEst`, and `RIDataPointEst`.
#' @param RIperc	(numeric) value specifying the percentiles, which define the reference interval (default c(0.025, 0.975))
#' @param data		(numeric) vector of unknown length, possibly containing NAs
#'
#' @returns 		(NULL) Prints the verification table to the console
#' 
#' @author Matthias Beck \email{matthias.beck.mb1@@roche.com}
#' 
printDataFractionWithinOutsideRI <- function(tab, RIperc, data){
	# Input validation
	stopifnot(is.data.frame(tab))
	stopifnot(is.numeric(RIperc), all(RIperc > 0 & RIperc < 1))

	# Filter verification table
	i_row <- sapply(tab$Percentile, function(p) any(abs(p - RIperc) < 1e-10))
	if(sum(i_row) == 0) {
		stop("No matching percentiles found in the verification table.")
	}

	required_cols <- c("Percentile", "RICandPointEst", "RIDataPointEst")
	stopifnot(all(required_cols %in% colnames(tab)))
	stopifnot(is.numeric(data))
	data <- data[is.finite(data)]

	tab <- tab[i_row, ]
	if(length(RIperc) > 2) {
		contains_025 <- any(abs(RIperc - 0.025) < 1e-10) && any(abs(tab$Percentile - 0.025) < 1e-10)
		contains_975 <- any(abs(RIperc - 0.975) < 1e-10) && any(abs(tab$Percentile - 0.975) < 1e-10)
		if(contains_025 && contains_975) {
			# keep the standard percentiles
			rows <- c(which.min(abs(tab$Percentile - 0.025)),
				which.min(abs(tab$Percentile - 0.975)))
		} else {
			# keep smallest and largest percentile
			rows <- c(which.min(tab$Percentile),
				which.max(tab$Percentile))
		}
		tab <- tab[rows, ]
		rm(rows)
	}

	# Helper function to calculate fractions
	calcFractions <- function(data, lower = -Inf, upper = Inf) {
		data_length <- length(data)
		below <- sum(data < lower, na.rm = TRUE) / data_length * 100
		within <- sum(data >= lower & data <= upper, na.rm = TRUE) / data_length * 100
		above <- sum(data > upper, na.rm = TRUE) / data_length * 100
		return(c(below, within, above))
	}
		

	if(length(tab$Percentile) == 1 && tab$Percentile < 0.5){

		RI   <- alignVec(c(tab$RIDataPointEst, tab$RICandPointEst))
		fracData <- calcFractions(data, lower = tab$RIDataPointEst)
		fracCand <- calcFractions(data, lower = tab$RICandPointEst)
		frac <- format(c(fracData, fracCand), nsmall = 1, digits = 2)

		rows <- list(
		list("RI Data", sprintf("[%s,    ]", RI[1]), paste0(frac[1], "%"), paste0(frac[2], "%")),
		list("RI Cand", sprintf("[%s,    ]", RI[2]), paste0(frac[4], "%"), paste0(frac[5], "%"))
		)

		headers <- c("Interval", "Reference Limits", "Below", "Within")

	} else if (length(tab$Percentile) == 1 && tab$Percentile > 0.5) {

		fracData <- calcFractions(data, upper = tab$RIDataPointEst)
		fracCand <- calcFractions(data, upper = tab$RICandPointEst)
		RI   <- alignVec(c(tab$RIDataPointEst, tab$RICandPointEst))
		frac <- format(c(fracData, fracCand), nsmall = 1, digits = 2)

		rows <- list(
		list("RI Data", sprintf("[   , %s]", RI[1]), paste0(frac[2], "%"), paste0(frac[3], "%")),
		list("RI Cand", sprintf("[   , %s]", RI[2]), paste0(frac[5], "%"), paste0(frac[6], "%"))
		)

		headers <- c("Interval", "Reference Limits", "Within", "Above")

	} else if ((length(tab$Percentile) == 1 && tab$Percentile == 0.5)) {
		fracData <- c(sum(data < tab$RIDataPointEst[1]), sum(data > tab$RIDataPointEst[1])) * 100 / length(data)
		fracCand <- c(sum(data < tab$RICandPointEst[1]), sum(data > tab$RICandPointEst[1])) * 100 / length(data)
		RI   <- alignVec(c(tab$RIDataPointEst, tab$RICandPointEst))
		frac <- format(c(fracData, fracCand), nsmall = 1, digits = 2)

		rows <- list(
		list("RI Data", sprintf("%s", RI[1]), paste0(frac[1], "%"), paste0(frac[2], "%")),
		list("RI Cand", sprintf("%s", RI[2]), paste0(frac[3], "%"), paste0(frac[4], "%"))
		)

		headers <- c("Interval", "Median", "Below", "Above")

	} else {
		fracData <- calcFractions(data, lower = tab$RIDataPointEst[1], upper = tab$RIDataPointEst[2])
		fracCand <- calcFractions(data, lower = tab$RICandPointEst[1], upper = tab$RICandPointEst[2])
		RI   <- alignVec(c(tab$RIDataPointEst, tab$RICandPointEst))
		frac <- format(c(fracData, fracCand), nsmall = 1, digits = 2)

		rows <- list(
		list("RI Data", sprintf("[%s , %s]", RI[1], RI[2]), paste0(frac[1], "%"), paste0(frac[2], "%"), paste0(frac[3], "%")),
		list("RI Cand", sprintf("[%s , %s]", RI[3], RI[4]), paste0(frac[4], "%"), paste0(frac[5], "%"), paste0(frac[6], "%"))
	)

		headers <- c("Interval", "Reference Limits", "Below", "Within", "Above")
	}

	# Calculate column widths
	stopifnot(length(headers) == length(rows[[1]]))
	col_widths <- numeric()
	for (i in seq_len(length(headers))) {
		col_widths <- c(col_widths, max(nchar(headers[i]), sapply(rows, function(row) nchar(row[[i]]))))
	}

	cat("\n")
	printTableHeader(headers, col_widths)
	printTableRows(rows, col_widths)
	cat("\n")
}


#' Calculate similarity of two reference intervals
#' 
#' @param RIdata		specifying the RI of the local population: (1) (object) of class \code{RWDRI} or (2) (numeric) representation of reference limits
#' @param RIcand		specifying the RI that needs to be verified: (1) (object) of class \code{RWDRI} or (2) (numeric) representation of reference limits
#' @param RIperc		(numeric) value specifying the percentiles, which define the reference interval (default c(0.025, 0.975))
#' @param pointEst		(character) specifying the point estimate determination: (1) using the full dataset ("fullDataEst"),
#' 						(2) calculating the median model from the bootstrap samples ("medianBS"), (2) works only if NBootstrap > 0
#' @param UMprop		(numeric) specifying the confidence level for the uncertainty margins
#' @param Overlap		(character) specifying the overlap criteria for the verification process. Options are: 
#' 							(1) uncertainty margins overlap "OverlapMargins" and 
#' 							(2) point estimates must be within the uncertainty margins "OverlapPointEst"
#' @param verbose		(logical) specifying if additional warning messages are printed
#' @param printResults	(logical) specifying if the results are printed to the console
#' @param ...			arguments to overwrite the default values of the Uncertainty Margin calculation
#' 
#' @returns				(data.frame) containing the similarity of the two reference intervals
#' 
#' @author Matthias Beck \email{matthias.beck.mb1@@roche.com}
#' 
#' @examples 
#' \dontrun{
#' 	example <- list(Mu = 3.41, Sigma = 0.504, Shift = 1, Lambda = 0.06,
#' 							  Method = "manual", roundingBase = NA)
#' 	class(example) <- "RWDRI"
#' 	test <- getRISimilarity(example, c(4, 55.5))
#' 	getRISimilarity(c(4, 55.5), c(6, 58))
#' 	getRISimilarity(c(4, 55.5), c(6, 58), UMprop = 0.95)
#' }
#' @export 
#' 
getRISimilarity <- function(RIdata, RIcand, RIperc = c(0.025, 0.975), pointEst = c("fullDataEst", "medianBS"), 
							UMprop = 0.9, Overlap = c("OverlapMargins", "OverlapPointEst"),
							printResults = TRUE, verbose = TRUE, ...){

	# check input
	stopifnot(inherits(RIdata, "RWDRI") || all(is.finite(RIdata)))
	stopifnot(inherits(RIcand, "RWDRI") || all(is.finite(RIcand)))
	stopifnot(is.numeric(RIperc))
	stopifnot(all(RIperc > 0 & RIperc < 1))
	stopifnot(is.logical(verbose))
	stopifnot(is.logical(printResults))

	Overlap <- match.arg(Overlap[1], c("OverlapMargins", "OverlapPointEst"))
	pointEst <- match.arg(pointEst[1], choices = c("fullDataEst", "medianBS"))

	 RIpercSorting <- sort(RIperc, index.return = TRUE)
	 RIperc <- RIpercSorting$x
	if(is.numeric(RIcand)){
		stopifnot(length(RIcand) == length(RIperc))
		RIcand <- RIcand[RIpercSorting$ix]
	}

	if(is.numeric(RIdata)){
		stopifnot(length(RIdata) == length(RIperc))
		RIdata <- RIdata[RIpercSorting$ix]
	}

	if(is.numeric(RIdata)){
		stopifnot(length(RIdata) == length(RIperc))
		stopifnot(length(RIdata) >= 1)
	}

	if(verbose && length(RIperc) == 1 && is.numeric(RIcand)){
		message("getRISimilarity involves interpolation as only one limit of RIcand is provided. Interpret results with caution.")
	}

	param <- list(...)
	checkInvalidArgs(param = param, validArgs = c("asymmetryCorr"), marginType = "VeRUS", verbose)


	asymmetryCorr <- TRUE
	if(!is.null(param$asymmetryCorr) && is.logical(param$asymmetryCorr) && !param$asymmetryCorr){
		asymmetryCorr <- FALSE
	}

	# get test arguments in the right format for the getRIMargins function
	tmp <- getVerificationArgs(RIdata = RIdata, RIcand = RIcand, RIperc = RIperc, pointEst = pointEst, verbose = verbose)
	dataTestArgs <- tmp$dataTestArgs
	candTestArgs   <- tmp$candTestArgs
	rm(tmp)

	stopifnot(length(dataTestArgs$RIperc) == length(candTestArgs$RIperc))
	stopifnot(all(abs(dataTestArgs$RIperc - candTestArgs$RIperc) <= 1e-10))

	ParamMarginsRICand <- list(RI = candTestArgs$RI, RIperc = candTestArgs$RIperc,
								lambda = candTestArgs$Lambda, shift = candTestArgs$Shift,
								n = NA,
								UMprop = UMprop,
								asymmetryCorr = asymmetryCorr
								)

  ParamMarginsRIData <- list(RI = dataTestArgs$RI, RIperc = dataTestArgs$RIperc,
								lambda = dataTestArgs$Lambda, shift = dataTestArgs$Shift,
								n = NA,
								UMprop = UMprop,
								asymmetryCorr =asymmetryCorr
								)

  # initialize similarity table
	SimilarityTab <- data.frame(
	    Percentile = ParamMarginsRICand$RIperc,
	    RIdata = ParamMarginsRIData$RI,
	    RIcand = ParamMarginsRICand$RI,
	    max_sample_size = NA,
	    s_value = NA)

	SimilarityTab <- SimilarityTab[SimilarityTab$Percentile %in% RIperc, ]

	# calculate similarity for each percentile
	for(i in seq_along(RIperc)){

		# calculate similarity
		upperLimit <- 1e7
		lowerLimit <- 0
		counter <- 0
		n <- NA
		perc <- RIperc[i]

		while(TRUE){

			n <- floor((upperLimit + lowerLimit)/2)
			ParamMarginsRICand$n <- ParamMarginsRIData$n <- n

			marginsRICand <- do.call(getRIMargins, ParamMarginsRICand)
			marginsRIData <- do.call(getRIMargins, ParamMarginsRIData)
			marN <- createVerificationTab(RIperc = RIperc, marginsRIdata = marginsRIData, marginsRIcand = marginsRICand)$RIVerificationTab
			overlapAtN <- marN[marN$Percentile == perc, Overlap]

			ParamMarginsRICand$n <- ParamMarginsRIData$n <- n + 1
			marginsRICand <- do.call(getRIMargins, ParamMarginsRICand)
			marginsRIData <- do.call(getRIMargins, ParamMarginsRIData)
			marN <- createVerificationTab(RIperc = RIperc, marginsRIdata = marginsRIData, marginsRIcand = marginsRICand)$RIVerificationTab
			overlapAtNPlus1 <- marN[marN$Percentile == perc, Overlap]

			if (overlapAtN && overlapAtNPlus1) {
				lowerLimit <- n
			} else if (overlapAtN && !overlapAtNPlus1) {
				break
			} else {
				upperLimit <- n
			}

			if(n < 25){
				n <- "n < 25" 
				break
			}
			if(n > 1e6 && counter > 30){
				n <- "n > 1e6"
				break
			}

			if(counter > 50){
				warning(paste0("Binary search did not converge [", perc, "]"))
				n <- "Error"
				break
			}
			counter <- counter + 1
		}

		SimilarityTab[SimilarityTab$Percentile == perc, "max_sample_size"] <- n
		if(is.numeric(n)){
			SimilarityTab[SimilarityTab$Percentile == perc, "s_value"] <- 6 / n # 0.05 * 120/n
		} else if( n == "n < 25"){
			SimilarityTab[SimilarityTab$Percentile == perc, "s_value"] <- paste0("> ", 6/25)
		} else if( n == "n > 1e6"){
			SimilarityTab[SimilarityTab$Percentile == perc, "s_value"] <- paste0("< ", 0.05 * 120/1e6)
		}

	}

	# find the minimum sample size for all percentiles
	sample_sizes <- SimilarityTab$max_sample_size
	sample_sizes[sample_sizes == "Error"] <- NA
	sample_sizes[sample_sizes == "n > 1e6"] <- Inf

	min_sample_size <- NA
	if(any(sample_sizes == "n < 25")) {
		min_sample_size <- "n < 25"
	} else {
		sample_sizes <- as.numeric(sample_sizes)
		min_sample_size <- min(sample_sizes[is.numeric(sample_sizes)], na.rm = TRUE)
	}

	if(printResults){
		printSimilarityTable(SimilarityTab, RIperc)
		cat("\n")
		cat("*Higher sample sizes indicate greater similarity between the reference limits.\n")
	}

	invisible(list(n = min_sample_size, SimilarityTab = SimilarityTab))
}



#' Helper function to print the results of getRISimilarity
#'
#' @description This function prints a formatted similarity table showing reference intervals and maximum sample sizes.
#'
#' @param SimilarityTab	(data.frame) containing the similarity table. Must include:
#'	- Percentile: The percentile values.
#'	- RIdata: Reference interval data points.
#'	- RIcand: Reference interval candidate points.
#'	- max_sample_size: Maximum sample sizes.
#'	- s_value: S-values indicating the level of similarity.
#'
#' @param RIperc		(numeric) value specifying the percentiles, which define the reference interval
#'
#' @returns				(NULL) Prints the similarity table to the console.
#' 
printSimilarityTable <- function(SimilarityTab, RIperc) {
	# Input validation
	stopifnot(is.data.frame(SimilarityTab))

	# Filter verification table
	i_row <- sapply(SimilarityTab$Percentile, function(p) any(abs(p - RIperc) < 1e-10))
	if(sum(i_row) == 0) {
		stop("No matching percentiles found in the verification table.")
	}

	SimilarityTab <- SimilarityTab[i_row, ]

	required_cols <- c("Percentile", "RIdata", "RIcand", "max_sample_size", "s_value")
	stopifnot(all(required_cols %in% colnames(SimilarityTab)))

	# Align numeric vectors
	aligned_nums <- list(
		RIdata = alignVec(SimilarityTab$RIdata),
		RIcand = alignVec(SimilarityTab$RIcand),
		max_sample_size = alignVec(c(SimilarityTab$max_sample_size, "Max Sample Size*"))[1:length(SimilarityTab$max_sample_size)],
		s_value = alignVec(c(SimilarityTab$s_value, "S-Value"))[1:length(SimilarityTab$s_value)]
	)

	# Generate formatted rows
	rows <- lapply(seq_len(nrow(SimilarityTab)), function(i) {
		list(
			Perc = sprintf("%4.1f%%", SimilarityTab$Percentile[i] * 100),
			`RI Data` = sprintf("%s", aligned_nums$RIdata[i]),
			`RI Cand` = sprintf("%s", aligned_nums$RIcand[i]),
			`Max Sample Size` = sprintf("%s", aligned_nums$max_sample_size[i]),
			`S-Value` = sprintf("%s", aligned_nums$s_value[i])
		)
	})

	# Calculate column widths
	headers <- c("Perc", "RI Data", "RI Cand", "Max Sample Size*", "S-Value")

	stopifnot(length(headers) == length(rows[[1]]))
	col_widths <- numeric()
	for (i in seq_len(length(headers))) {
		col_widths <- c(col_widths, max(nchar(headers[i]), sapply(rows, function(row) nchar(row[[i]]))))
	}

	# Print table header
	cat("\nSimilarity Table\n")
	cat(paste(rep("-", 16), collapse = ""), "\n\n")
	printTableHeader(headers, col_widths)

	# Print table rows
	printTableRows(lapply(rows, function(r) unlist(r)), col_widths)

	cat("\n")
}

#' Helper function to align the content of printed tables
#' 
#' @description This function aligns vectors by padding with spaces to the right
#' 
#' @param x			(numeric or character) vector of values to be aligned
#' @param digits	(integer) number of digits to be displayed for numeric values
#' 
#' @returns			(character) vector of aligned values
#' 
#' @author Christopher Rank \email{christopher.rank@@roche.com}
#' 
alignVec <- function(x, digits = 3) {
	stopifnot(is.numeric(digits), digits > 0, as.integer(digits) == digits)

	# Initialize variables
	lengthChar <- digits
	isNumeric <- sapply(x, function(val) !is.na(suppressWarnings(as.numeric(val))))

	# Get maximum width needed
	for(i in seq_along(x)) {
		if(isNumeric[i]){
			tmp <- nchar(as.character(signif(as.numeric(x[i]), digits)))
		} else {
			tmp <- nchar(as.character(trimws(x[i])))
		}
		if(is.na(tmp)) tmp <- 2
		lengthChar <- max(lengthChar, tmp)
	}

	# Format each element
	xChar <- character(length(x))

	for (i in seq_along(x)) {
		if(is.na(x[i])) {
			xChar[i] <- format("", width = lengthChar, justify = "right")
		} else if(isNumeric[i]) {
				xChar[i] <- format(signif(as.numeric(x[i]), digits), 
						width = lengthChar, 
							justify = "right")
		} else {
			xChar[i] <- format(as.character(x[i]), 
						width = lengthChar, 
						justify = "right")
		}
	}
	return(xChar)
}