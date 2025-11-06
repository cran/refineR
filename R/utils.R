#' One-parameter Box-Cox transformation.
#' 
#' @param x				(numeric) data to be transformed
#' @param lambda		(numeric) Box-Cox transformation parameter
#' 
#' @return (numeric) vector with Box-Cox transformation of x
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}

BoxCox <- function(x, lambda) {

	if(abs(lambda) < 1e-20)
		x <- log(x)
	else
		x <- (x^lambda-1) / lambda
	x
}


#' Inverse of the one-parameter Box-Cox transformation.
#' 
#' @param x				(numeric) data to be transformed
#' @param lambda		(numeric) Box-Cox transformation parameter
#' 
#' @return (numeric) vector with inverse Box-Cox transformation of x
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}

invBoxCox <- function(x, lambda) {

	if(abs(lambda) < 1e-20)
		x <- exp(x)
	else
		x <- (lambda*x+1)^(1/lambda) # if(lambda*x+1) is negative & lambda < 1-> result is NaN
	x
}


#' Approximate calculation of CDF of normal distribution.
#' 
#' @param q			(numeric) vector of quantiles of data points
#' @param pNormVal	(numeric) vector of lookup table for pNorm
#' @param mean		(numeric) vector of mean values
#' @param oneOverSd	(numeric) reciprocal vector of sd values
#' @param oneOverH	(numeric) defining the precision of the approximation
#' 
#' @return (numeric) vector of approximate CDFs of normal distribution. 
#'
#' @author Christopher Rank \email{christopher.rank@@roche.com}

pnormApprox <- function(q, pNormVal, mean = 0, oneOverSd = 1, oneOverH = 10) {	

	index <- ((q-mean)*oneOverSd + 5.2)*oneOverH
	index[index < 0] <- 0
	index[index >= (length(pNormVal)-2)] <- length(pNormVal)-1.5

	indexInt <- as.integer(index)

	w <- index-indexInt

	res <- (1-w)*pNormVal[indexInt+1] + w*pNormVal[indexInt+2]

	return(res)
}


#' Estimate density of distribution employing the R package "ash" using R-wrapper function.
#' 
#' @param x			(numeric) vector of data points
#' @param ab		(numeric) vector of lower and higher truncation limit of density estimation
#' @param nbin		(integer) specifying the number of bins used for density estimation
#' @param m			(integer) specifying the width of the smoothing kernel(s) used for density estimation
#' @param kopt	    (integer) vector specifying the smoothing kernel
#' @param normToAB	(logical) specifying if the density is normed to the interval ab or to all data points in x
#' 
#' @return (list) with density estimation (x values, y values, m and ab). 
#'
#' @author Christopher Rank \email{christopher.rank@@roche.com}, Tatjana Ammer \email{tatjana.ammer@@roche.com}

ashDensity <- function(x, ab, nbin, m, kopt = c(2, 1), normToAB = FALSE) {	

	invisible(capture.output(d <- ash1(bins = bin1(x = x, ab = ab, nbin = nbin), m = m, kopt = kopt)))

	if(normToAB)
		PCorr <- 1.0
	else
		PCorr <- sum(x >= ab[1] & x < ab[2])/length(x)

	list(x = d$x, y = d$y*PCorr, m = m, ab = ab, nbin = nbin, PCorr = PCorr)
}



#' Estimate rounding base of the input data. 
#' 
#' @param x		(numeric) vector of data points
#' 
#' @return (numeric) with estimated rounding base (e.g. 0.001 when rounded to 3 digits)
#' 
#' @author Christopher Rank \email{christopher.rank@@roche.com}, Tatjana Ammer \email{tatjana.ammer@@roche.com}

findRoundingBase <- function(x) {

	ab <- as.numeric(quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))

	# generate table with difference between neighboring values 
	diffVal  <- table(round(diff(sort(unique(x[x>=ab[1] & x<=ab[2]]))), digits=10))
	# remove number of those that have a difference smaller than e-10
	diffVal <- diffVal[names(diffVal) != "0"]

	# select differences that occur in >= 1% of cases
	diffVal  <- diffVal[diffVal>=0.01*sum(diffVal)]

	# find multiples of rounding base
	if(length(diffVal) > 1)
	{
		for(i in 1:length(diffVal))
		{
			if(diffVal[i] > 0)
			{
				ratio <- as.numeric(names(diffVal))/as.numeric(names(diffVal[i]))

				selection <- which(ratio>1.5 & abs(round(ratio, digits=9)-round(ratio))<1e-20)

				diffVal[i] <- diffVal[i] + sum(diffVal[selection])

				diffVal[selection] <- 0
			}
		}
	}

	# select differences that occur in >= 10% of cases
	diffVal  <- diffVal[diffVal>=0.1*sum(diffVal)]

	roundingBase <- NA

	# if dataset has only finite number of unique values (e.g. rounded data)
	if (length(diffVal) > 0) {

		# determine rounding base (maximum step size of discrete values)
		roundingBase <- max(as.numeric(names(diffVal)))
	}

	return(roundingBase)
}


#' Function to calculate the RI from a set of parameters
#' 
#' @param mu				(numeric) mean of the distribution
#' @param sigma			(numeric) standard deviation of the distribution
#' @param lambda			(numeric) Box-Cox transformation parameter
#' @param shift			(numeric) shift of the distribution
#' @param RIperc		(numeric) vector of percentiles for which the reference interval should be calculated
#' 
#' @return (numeric) vector of reference interval values for the given percentiles
#' 
#' @references Freeman, R. Modarres / Statistics & Probability Letters 76 (2006) P 767
#' @author Matthias Beck \email{matthias.beck.mb1@@roche.com}, Christopher Rank \email{christopher.rank@@roche.com}, Tatjana Ammer \email{tatjana.ammer@@roche.com}
#'

cdfTruncatedBoxCox <- function(mu, sigma, lambda, shift, RIperc){

	stopifnot(is.numeric(RIperc), min(RIperc)>=0, max(RIperc)<=1)
	stopifnot(is.numeric(mu), is.numeric(sigma), is.numeric(lambda), is.numeric(shift))
	if(any(!is.finite(c(mu, sigma, lambda, shift))) || sigma <= 0)
	{
		return(rep(NA, length(RIperc)))
	}

	RIperc <- sort(RIperc)

	# formula for truncated normal distribution
	if(lambda >= 0) {
		RI <- pnorm(-1/lambda, mean=mu, sd=sigma) + RIperc*(1 - pnorm(-1/lambda, mean=mu, sd=sigma))
		RI <- qnorm(RI, mean=mu, sd=sigma)
	} else {
		# for lambda < 0
		RI <- qnorm(RIperc, mean=mu, sd=sigma)
	}

	RI <- invBoxCox(RI, lambda = lambda)
	RI <- RI + shift
	RI[RI<0 | is.na(RI)] <- 0
	return(RI) # remove negative values
}


#' Method to calculate reference intervals (percentiles) for objects of class 'RWDRI'
#' 
#' @param x				(object) of class 'RWDRI'
#' @param RIperc		(numeric) value specifying the percentiles, which define the reference interval
#' @param CIprop		(numeric) value specifying the central region for estimation of confidence intervals
#' @param pointEst		(character) specifying the point estimate determination: (1) using the full dataset ("fullDataEst"),
#' 						(2) calculating the median model from all bootstrap samples ("medianBS"), (2) works only if NBootstrap > 0
#' @param Scale			(character) specifying if percentiles are calculated on the original scale ("or") or the transformed scale ("tr") or the z-Score scale ("z")
#' @param UMprop		(numeric) value specifying the central region for estimation of uncertainty margins
#' @param ...       calcUCMargins argument (logical) disabling the calculation of uncertainty margins when set to FALSE,
#' 									n argument (integer) specifying the theoretical sample size used for uncertainty margin calculation default (n = 120),	
#' 									asymmetryCorr argument (logical) disabling the asymmetry correction when set to FALSE
#'
#' @return				(data.frame) with columns for percentile, point estimate, bootstrap-based confidence intervals and uncertainty margins.
#'
#' @author Christopher Rank \email{christopher.rank@@roche.com}, Tatjana Ammer \email{tatjana.ammer@@roche.com}, Matthias Beck \email{matthias.beck.mb1@@roche.com}
#'
getRI <- function(x, RIperc = c(0.025, 0.975), CIprop = 0.95, pointEst = c("fullDataEst", "medianBS"), Scale = c("original", "transformed", "zScore"), UMprop = 0.90, ...) {

	stopifnot(class(x) == "RWDRI")
	stopifnot(is.numeric(RIperc), min(RIperc)>=0, max(RIperc)<=1)
	stopifnot(length(RIperc) > 0)
	stopifnot(is.numeric(CIprop), length(CIprop)==1, CIprop>=0, CIprop<=1)
	pointEst <- match.arg(pointEst[1], choices = c("fullDataEst", "medianBS"))
	Scale    <- match.arg(Scale[1], choices = c("original", "transformed", "zScore"))

	calcUCMargins <- TRUE
	asymmetryCorr <- TRUE
	n <- 120

	args <- list(...)
	if("UCMargins" %in% names(args) && !(args$UCMargins))
	{

		calcUCMargins <- args$UCMargins
	}
	if("asymmetryCorr" %in% names(args) && !(args$asymmetryCorr))
	{
		asymmetryCorr <- args$asymmetryCorr
	}
	if("n" %in% names(args) && is.numeric(args$n) && length(args$n) == 1 && args$n > 0)
	{
		n <- args$n
	}

	RIperc	 <- sort(RIperc)
	RIResult <- data.frame(Percentile = RIperc, PointEst = NA, CILow = NA, CIHigh = NA, UMLow = NA, UMHigh = NA)

	if (!is.na(x$Mu) && !is.na(x$Sigma) && !is.na(x$Lambda) && !is.na(x$Shift)) {

		BSPerformed  <- (length(x$MuBS)>0 && length(x$SigmaBS)>0 && length(x$LambdaBS)>0 && length(x$ShiftBS)>0)

		# extract model parameters
		lambda <- ifelse(pointEst=="medianBS" && BSPerformed, x$LambdaMed, x$Lambda)
		mu 	   <- ifelse(pointEst=="medianBS" && BSPerformed, x$MuMed, 	  x$Mu)
		sigma  <- ifelse(pointEst=="medianBS" && BSPerformed, x$SigmaMed,  x$Sigma)
		shift  <- ifelse(pointEst=="medianBS" && BSPerformed, x$ShiftMed,  x$Shift)

		RI <- cdfTruncatedBoxCox(mu = mu, sigma = sigma, lambda = lambda, shift = shift, RIperc = RIperc)

		RIforMargins <- RI[is.finite(RI)]
		if(length(RIforMargins) == 0)
		{
			calcUCMargins <- BSPerformed <- FALSE
		}

		if(calcUCMargins) {
			stopifnot(is.numeric(UMprop), length(UMprop)==1, UMprop>=0, UMprop<=1)
			# impute RI in case of a single value
			if(length(RIforMargins) == 1) {
				RIpercTmp <- RIperc[RI %in% RIforMargins] # should be length 1
				# calculate symmetrical RI
				if(RIpercTmp >= 0.5 - 1e-3 && RIpercTmp <= 0.5 + 1e-3)
				{
					RIpercMargin <- c(0.1, RIpercTmp)
				} else {
					RIpercMargin <- c(RIpercTmp, 1 - RIpercTmp)
				}
			} else {
				RIpercMargin <- RIperc
			}

			RIforMargins <- cdfTruncatedBoxCox(mu = mu, sigma = sigma, lambda = lambda, shift = shift, RIperc = RIpercMargin)
			RIpercMargin <- sort(RIpercMargin) # cdfTruncatedBoxCox sorts the percentiles
			RIMargins <- getRIMargins(RI = RIforMargins, RIperc=RIpercMargin, UMprop=UMprop, lambda=lambda, shift=shift, asymmetryCorr = asymmetryCorr, n =n)
			RIMargins <- RIMargins[abs(RIMargins$Percentile - RIperc) <= 1e-5, ]

			RIResult$UMLow  <- RIMargins$MarginLow
			RIResult$UMHigh <- RIMargins$MarginHigh

			if (Scale == "transformed" | Scale == "zScore") {

				RIResult$UMLow  <- suppressWarnings(BoxCox(RIResult$UMLow-shift,  lambda = lambda))
				RIResult$UMHigh <- suppressWarnings(BoxCox(RIResult$UMHigh-shift, lambda = lambda))

				if(Scale == "zScore")
				{
					RIResult$UMLow  <- (RIResult$UMLow  - mu) / sigma
					RIResult$UMHigh <- (RIResult$UMHigh - mu) / sigma
				}
			}

		} else {
			RIResult$UMLow  <- RIResult$UMHigh <- NA
		}

		if (Scale == "transformed" | Scale == "zScore") {
			
			RI <- suppressWarnings(BoxCox(RI-shift, lambda = lambda))			
			
			if(Scale == "zScore") {
				RI <- (RI - mu) / sigma
			}
		}
		
		RIResult$PointEst <- RI

		# reference intervals for Bootstrap samples
		if (BSPerformed) {
			NBootstrap <- min(x$NBootstrap, length(x$MuBS), length(x$SigmaBS), length(x$LambdaBS), length(x$ShiftBS))

			for (p in seq_along(RIperc)) {
				# Calculate the reference interval for each bootstrap sample
				RIBS <- sapply(seq_len(NBootstrap), function(i) {
					cdfTruncatedBoxCox(
						mu = x$MuBS[i],
						sigma = x$SigmaBS[i],
						lambda = x$LambdaBS[i],
						shift = x$ShiftBS[i],
						RIperc = RIperc[p]
					)
				})

				RIBS <- RIBS[!is.na(RIBS)]

				if (Scale == "transformed" | Scale == "zScore") {

					RIBS <- suppressWarnings(BoxCox(RIBS-shift, lambda = lambda))

					if(Scale == "zScore")
						RIBS <- (RIBS - mu) / sigma
				} 

				RIResult$CILow[p]  <- as.numeric(quantile(x = RIBS, probs = (1-CIprop)/2, na.rm = TRUE))
				RIResult$CIHigh[p] <- as.numeric(quantile(x = RIBS, probs = 1-(1-CIprop)/2, na.rm = TRUE))
			}
		}
	}

	return(RIResult)
}


#' Standard print method for objects of class 'RWDRI'
#' 
#' @param x					(object) of class 'RWDRI'
#' @param RIperc			(numeric) value specifying the percentiles, which define the reference interval
#' @param uncertaintyRegion	(character) specifying the type of the uncertainty region around point estimates 
#' @param CIprop			(numeric) value specifying the central region for estimation of confidence intervals
#' @param UMprop			(numeric) value specifying the central region for estimation of uncertainty margins
#' @param pointEst			(character) specifying the point estimate determination: (1) using the full dataset ("fullDataEst"),
#' 							(2) calculating the median model from all bootstrap samples ("medianBS"), (2) works only if NBootstrap > 0
#' @param ...			additional arguments passed forward to other functions.
#'
#' @return				No return value. Instead, a summary is printed.
#' @examples
#' \dontrun{
#' x <- refineR::findRI(refineR::testcase1)
#' print(x, uncertaintyRegion = "bootstrapCI")
#' print(x, uncertaintyRegion = "uncertaintyMargin")
#' }
#'
#' @author Christopher Rank \email{christopher.rank@@roche.com}
#'
#' @method print RWDRI

print.RWDRI <- function(x, RIperc = c(0.025, 0.975), uncertaintyRegion = c("bootstrapCI", "uncertaintyMargin"), CIprop = 0.95, UMprop = 0.9, pointEst = c("fullDataEst", "medianBS"), ...) {

	stopifnot(class(x) == "RWDRI")
	stopifnot(is.numeric(RIperc), min(RIperc)>=0, max(RIperc)<=1)
	uncertaintyRegion <- match.arg(uncertaintyRegion[1], choices = c("bootstrapCI", "uncertaintyMargin"))
	pointEst        <- match.arg(pointEst[1], choices = c("fullDataEst", "medianBS"))

	args <- list(...)
	if("n" %in% names(args) && is.numeric(args$n) && length(args$n) == 1 && args$n > 0)
	{
		n <- as.integer(args$n)
	} else {
		n <- as.integer(120)
	}

	if(uncertaintyRegion == "bootstrapCI")
	{
		stopifnot(is.numeric(CIprop), length(CIprop)==1, CIprop>=0, CIprop<=1)
	} else if(is.null(CIprop) && uncertaintyRegion == "uncertaintyMargin") {
		stopifnot(is.numeric(UMprop), length(UMprop)==1, UMprop>=0, UMprop<=1)
	}

	# calculate reference intervals
	if(uncertaintyRegion == "bootstrapCI") {
		RI <- getRI(x = x, RIperc = RIperc, CIprop = CIprop, UMprop = UMprop, pointEst = pointEst, UCMargins = FALSE)
	} else if(uncertaintyRegion == "uncertaintyMargin") {
		RI <- getRI(x = x, RIperc = RIperc, CIprop = CIprop, UMprop = UMprop, pointEst = pointEst, UCMargins = TRUE, n = n)
	}


	uncertaintyRegionExists <- (uncertaintyRegion=="bootstrapCI" && any(!is.na(RI$CILow))     && any(!is.na(RI$CIHigh))) |
											 (uncertaintyRegion=="uncertaintyMargin" && any(!is.na(RI$UMLow)) && any(!is.na(RI$UMHigh)))

	cat("\nEstimated Reference Interval\n")
	cat("------------------------------------------------\n")

	#check if reference intervals are na
	for (i in seq_along(RIperc)) {
		limit     <- "     median ["
		if(RIperc[i] < 0.5)
			limit <- "lower limit ["
		if(RIperc[i] > 0.5)
			limit <- "upper limit ["

		cat(paste0(limit, ifelse(RIperc[i]*100<10, " ", ""), format(round(RIperc[i]*100, 1), nsmall = 1), "%]: ", signif(RI$PointEst[i], 3)))
		if(uncertaintyRegion=="bootstrapCI" & !is.na(RI$CILow[i]) & !is.na(RI$CIHigh[i]))
			cat(paste0(" (", signif(RI$CILow[i], 3), "; ", signif(RI$CIHigh[i], 3), ")\n"))
		else if(uncertaintyRegion=="uncertaintyMargin" & !is.na(RI$UMLow[i]) & !is.na(RI$UMHigh[i]))
			cat(paste0(" (", signif(RI$UMLow[i], 3), "; ", signif(RI$UMHigh[i], 3), ")\n"))	
		else
			cat("\n")
	}

	if(uncertaintyRegionExists){
		ucText <- ifelse(uncertaintyRegion == "bootstrapCI", "bootstrapped CI",
		 paste0("uncertainty margin [n = ", n, "]"))
		if(uncertaintyRegion == "uncertaintyMargin") prop <- signif(UMprop*100, 2)
		if(uncertaintyRegion == "bootstrapCI")       prop <- signif(CIprop*100, 2)
		cat(paste0("\nWidth of the ", ucText, ": ", prop, "% \n"))
	}
	cat("\nModel Parameters\n")
	cat("------------------------------------------------\n")	

	cat(paste0("     method: ", x$Method, " (v", x$PkgVersion, ")\n"))
	cat(paste0("      model: ", x$Model, "\n"))
	cat(paste0("     N data: ", length(x$Data), "\n"))

	if(!is.null(x$NBootstrap) && x$NBootstrap > 0)
		cat(paste0("N bootstrap: ", x$NBootstrap, "\n"))

	if(!is.na(x$roundingBase))
		cat(paste0("    rounded: yes (base: ", x$roundingBase, ")\n"))
	else
		cat(paste0("    rounded: no\n"))

	if (!is.null(x$AgeMin) && !is.null(x$AgeMax))
		cat(paste0("  Age range: ", x$AgeMin, " to ", x$AgeMax, " years\n"))

	if (!is.null(x$Group))
		cat(paste0("     Gender: ", paste(x$Group, collapse=", "), "\n"))

	BSPerformed  <- (length(x$MuBS)>0 && length(x$SigmaBS)>0 && length(x$LambdaBS)>0 && length(x$ShiftBS)>0)

	cat(paste0("  point est: ", ifelse(pointEst=="medianBS" && BSPerformed, "medianBS", "fullDataEst"), "\n"))

	if(uncertaintyRegionExists)
		cat(paste0("uncertainty: ", uncertaintyRegion, "\n"))

	cat(paste0("     lambda: ", ifelse(pointEst=="medianBS" && BSPerformed, signif(x$LambdaMed, 3), signif(x$Lambda, 3)), "\n"))
	cat(paste0("         mu: ", ifelse(pointEst=="medianBS" && BSPerformed, signif(x$MuMed, 3), 	   signif(x$Mu, 3)), 	 "\n"))
	cat(paste0("      sigma: ", ifelse(pointEst=="medianBS" && BSPerformed, signif(x$SigmaMed, 3),  signif(x$Sigma, 3)),  "\n"))
	cat(paste0("      shift: ", ifelse(pointEst=="medianBS" && BSPerformed, signif(x$ShiftMed, 3),  signif(x$Shift, 3)),  "\n"))
	if(!is.null(x$CostMed) || !is.null(x$Cost))
		cat(paste0("       cost: ", ifelse(pointEst=="medianBS" && BSPerformed, signif(x$CostMed, 3),   signif(x$Cost, 3)),   "\n"))
	if(!is.null(x$PMed) || !is.null(x$P))
	cat(paste0("NP fraction: ", ifelse(pointEst=="medianBS" && BSPerformed, signif(x$PMed, 3), 	   signif(x$P, 3)), 	 "\n"))
}
