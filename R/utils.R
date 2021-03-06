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
		x <- (lambda*x+1)^(1/lambda) # if(lambda*x+1) is negative -> result is NaN
		
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


#' Method to calculate reference intervals (percentiles) for objects of class 'RWDRI'
#' 
#' @param x				(object) of class 'RWDRI'
#' @param RIperc		(numeric) value specifying the percentiles, which define the reference interval
#' @param CIprop		(numeric) value specifying the central region for estimation of confidence intervals
#' @param pointEst		(character) specifying the point estimate determination: (1) using the full dataset ("fullDataEst"),
#' 						(2) calculating the median from all bootstrap samples ("medianBS"), (2) works only if NBootstrap > 0
#' 						(3) calculating the mean from all bootstrap samples ("meanBS"), (3) works only if NBootstrap > 0
#' @param Scale			(character) specifying if percentiles are calculated on the original scale ("Or") or the transformed scale ("Tr")
#' 
#' @return				(data.frame) with columns for percentile, point estimate and confidence intervals. 
#' 
#' @author Christopher Rank \email{christopher.rank@@roche.com}, Tatjana Ammer \email{tatjana.ammer@@roche.com}

getRI <- function(x, RIperc = c(0.025, 0.975), CIprop = 0.95, pointEst = c("fullDataEst", "medianBS", "meanBS"), Scale = c("original", "transformed")) {
	
	stopifnot(class(x) == "RWDRI")
	stopifnot(is.numeric(RIperc))
	stopifnot(is.numeric(CIprop))
	pointEst <- match.arg(pointEst[1], choices = c("fullDataEst", "medianBS", "meanBS"))
	Scale    <- match.arg(Scale[1], choices = c("original", "transformed"))
	
	RIperc	 <- sort(RIperc)
	RIResult <- data.frame(Percentile = RIperc, PointEst = NA, CILow = NA, CIHigh = NA)
	
	if (!is.na(x$Mu) & !is.na(x$Sigma) & !is.na(x$Lambda) & !is.na(x$Shift)) {
			
		RI <- qnorm(p = RIperc, mean = x$Mu, sd = x$Sigma)
	
		if (Scale == "original") {
			RI <- invBoxCox(RI, lambda = x$Lambda)
			RI <- RI + x$Shift
						
			RI[RI<0] <- 0		
			RI[is.na(RI)] <- 0 
		}
		
		RIResult$PointEst <- RI
	
		# reference intervals for Bootstrap samples
		if (Scale == "original" & length(x$MuBS) > 0 & length(x$SigmaBS) > 0 & length(x$LambdaBS) > 0 & length(x$ShiftBS) > 0 & length(x$CostBS) > 0) {
			
			for (i in 1:length(RIperc)) {
				
				RIBS <- qnorm(p = RIperc[i], mean = x$MuBS, sd = x$SigmaBS)
								
				for (l in 1:length(RIBS)) {
					if (!is.na(x$MuBS[l]) & !is.na(x$SigmaBS[l]) & !is.na(x$LambdaBS[l]) & !is.na(x$ShiftBS[l])) {
						RIBS[l]	<- max(0, invBoxCox(RIBS[l], x$LambdaBS[l]) + x$ShiftBS[l], na.rm = TRUE)
						
					} else {
						RIBS[l] <-NA
					}
				}				
				
				RIResult$CILow[i]  <- as.numeric(quantile(x = RIBS, probs = (1-CIprop)/2, na.rm = TRUE))
				RIResult$CIHigh[i] <- as.numeric(quantile(x = RIBS, probs = 1-(1-CIprop)/2, na.rm = TRUE))
				
				if(pointEst == "medianBS") {
					RIResult$PointEst[i] <- median(RIBS, na.rm = TRUE)					
				} 		
				if(pointEst == "meanBS") {
					RIResult$PointEst[i] <- mean(RIBS, na.rm = TRUE)					
				}
			}
		}
	}
	
	return(RIResult)	
}


#' Standard print method for objects of class 'RWDRI'
#' 
#' @param x				(object) of class 'RWDRI'
#' @param RIperc		(numeric) value specifying the percentiles, which define the reference interval
#' @param CIprop		(numeric) value specifying the central region for estimation of confidence intervals
#' @param pointEst		(character) specifying the point estimate determination: (1) using the full dataset ("fullDataEst"),
#' 						(2) calculating the median from all bootstrap samples ("medianBS"), (2) works only if NBootstrap > 0
#' 						(3) calculating the mean from all bootstrap samples ("meanBS"), (3) works only if NBootstrap > 0
#' @param ...			additional arguments passed forward to other functions.
#' 
#' @return				No return value. Instead, a summary is printed.
#' 
#' @author Christopher Rank \email{christopher.rank@@roche.com}
#' 
#' @method print RWDRI

print.RWDRI <- function(x, RIperc = c(0.025, 0.975), CIprop = 0.95, pointEst = c("fullDataEst", "medianBS", "meanBS"), ...) {
	
	stopifnot(class(x) == "RWDRI")
	stopifnot(is.numeric(RIperc))
	stopifnot(is.numeric(CIprop))
	
	pointEst <- match.arg(pointEst[1], choices = c("fullDataEst", "medianBS", "meanBS"))
	
	# calculate reference intervals
	RI <- getRI(x = x, RIperc = RIperc, CIprop = CIprop, pointEst = pointEst)
	
	cat("\nReference Intervals\n")
	cat("------------------------------------------------\n")
	
	#check if reference intervals are na
	for (i in 1:length(RIperc)) {
		limit     <- "     median ["
		if(RIperc[i] < 0.5)
			limit <- "lower limit ["
		if(RIperc[i] > 0.5)
			limit <- "upper limit ["
		
		cat(paste0(limit, ifelse(RIperc[i]*100<10, " ", ""), format(round(RIperc[i]*100, 1), nsmall = 1), "% perc]: ", signif(RI$PointEst[i], 3))) 		
		if(!is.na(RI$CILow[i]) & !is.na(RI$CIHigh[i]))
			cat(paste0(" (", signif(RI$CILow[i], 3), "; ", signif(RI$CIHigh[i], 3), ")\n"))	
		else
			cat("\n")			
	}		
	
	cat("\nModel Parameters\n")
	cat("------------------------------------------------\n")	
	
	cat(paste0("     method: ", x$Method, " (v", x$PkgVersion, ")\n"))
	cat(paste0("      model: ", x$Model, "\n"))
	cat(paste0("     N data: ", length(x$Data), " (data points)\n"))

	if(!is.na(x$roundingBase))
		cat(paste0("    rounded: yes (base: ", x$roundingBase, ")\n"))
	else
		cat(paste0("    rounded: no\n"))
	
	if (!is.null(x$AgeMin) & !is.null(x$AgeMax))
		cat(paste0("  Age range: ", x$AgeMin, " to ", x$AgeMax, " years\n"))
	
	if (!is.null(x$Group))
		cat(paste0("     Gender: ", paste(x$Group, collapse=", "), "\n"))		
	
	cat(paste0("     lambda: ", signif(x$Lambda, 3), "\n"))
	cat(paste0("         mu: ", signif(x$Mu, 3), "\n"))
	cat(paste0("      sigma: ", signif(x$Sigma, 3), "\n"))
	cat(paste0("      shift: ", signif(x$Shift, 3), "\n"))
	cat(paste0("       cost: ", signif(x$Cost, 3), "\n"))
	cat(paste0("NP fraction: ", signif(x$P, 3), "\n"))	
}

