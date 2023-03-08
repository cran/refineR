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
#' 						(2) calculating the median model from all bootstrap samples ("medianBS"), (2) works only if NBootstrap > 0 						
#' @param Scale			(character) specifying if percentiles are calculated on the original scale ("or") or the transformed scale ("tr") or the z-Score scale ("z")
#' 
#' @return				(data.frame) with columns for percentile, point estimate and confidence intervals. 
#' 
#' @author Christopher Rank \email{christopher.rank@@roche.com}, Tatjana Ammer \email{tatjana.ammer@@roche.com}

getRI <- function(x, RIperc = c(0.025, 0.975), CIprop = 0.95, pointEst = c("fullDataEst", "medianBS"), Scale = c("original", "transformed", "zScore")) {
	
	stopifnot(class(x) == "RWDRI")
	stopifnot(is.numeric(RIperc) & min(RIperc)>=0 & max(RIperc)<=1)
	stopifnot(is.numeric(CIprop) & length(CIprop)==1 & CIprop>=0 & CIprop<=1)
	pointEst <- match.arg(pointEst[1], choices = c("fullDataEst", "medianBS"))
	Scale    <- match.arg(Scale[1], choices = c("original", "transformed", "zScore"))
	
	RIperc	 <- sort(RIperc)
	RIResult <- data.frame(Percentile = RIperc, PointEst = NA, CILow = NA, CIHigh = NA)
	
	if (!is.na(x$Mu) & !is.na(x$Sigma) & !is.na(x$Lambda) & !is.na(x$Shift)) {
			
		BSPerformed  <- (length(x$MuBS)>0 & length(x$SigmaBS)>0 & length(x$LambdaBS)>0 & length(x$ShiftBS)>0)
		
		# extract model parameters
		lambda <- ifelse(pointEst=="medianBS" & BSPerformed, x$LambdaMed, x$Lambda)
		mu 	   <- ifelse(pointEst=="medianBS" & BSPerformed, x$MuMed, 	  x$Mu)
		sigma  <- ifelse(pointEst=="medianBS" & BSPerformed, x$SigmaMed,  x$Sigma)
		shift  <- ifelse(pointEst=="medianBS" & BSPerformed, x$ShiftMed,  x$Shift)
		
		# formula for truncated normal distribution
		RI <- pnorm(-1/lambda, mean=mu, sd=sigma) + RIperc*(1 - pnorm(-1/lambda, mean=mu, sd=sigma))
		RI <- qnorm(RI, mean=mu, sd=sigma)
			
		if (Scale == "original") {
			RI <- invBoxCox(RI, lambda = lambda)
			RI <- RI + shift
						
			RI[RI<0] <- 0		
			RI[is.na(RI)] <- 0 
			
		} else if(Scale == "zScore")
		{
			RI <- qnorm(RIperc, mean=0, sd=1)
		}
		
		RIResult$PointEst <- RI
	
		# reference intervals for Bootstrap samples
		if (BSPerformed) {
			
			for (i in 1:length(RIperc)) {
											
				# formula for truncated normal distribution
				RIBS <- pnorm(-1/x$LambdaBS, mean=x$MuBS, sd=x$SigmaBS) + RIperc[i]*(1 - pnorm(-1/x$LambdaBS, mean=x$MuBS, sd=x$SigmaBS))
				RIBS <- qnorm(RIBS, mean=x$MuBS, sd=x$SigmaBS)					
								
				for (l in 1:length(RIBS)) 
				{					
					RIBS[l]	<- max(0, invBoxCox(RIBS[l], x$LambdaBS[l]) + x$ShiftBS[l], na.rm = TRUE)					
				}				
				
				if (Scale == "transformed" | Scale == "zScore" ) {					
			
					RIBS <- suppressWarnings(BoxCox(RIBS-shift, lambda = lambda))				
					
					if(Scale == "zScore")
						RIBS <- (RIBS - mu) / sigma					
				} 
				
				RIResult$CILow[i]  <- as.numeric(quantile(x = RIBS, probs = (1-CIprop)/2, na.rm = TRUE))
				RIResult$CIHigh[i] <- as.numeric(quantile(x = RIBS, probs = 1-(1-CIprop)/2, na.rm = TRUE))				
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
#' 						(2) calculating the median model from all bootstrap samples ("medianBS"), (2) works only if NBootstrap > 0
#' @param ...			additional arguments passed forward to other functions.
#' 
#' @return				No return value. Instead, a summary is printed.
#' 
#' @author Christopher Rank \email{christopher.rank@@roche.com}
#' 
#' @method print RWDRI

print.RWDRI <- function(x, RIperc = c(0.025, 0.975), CIprop = 0.95, pointEst = c("fullDataEst", "medianBS"), ...) {
	
	stopifnot(class(x) == "RWDRI")
	stopifnot(is.numeric(RIperc) & min(RIperc)>=0 & max(RIperc)<=1)
	stopifnot(is.numeric(CIprop) & length(CIprop)==1 & CIprop>=0 & CIprop<=1)
	
	pointEst <- match.arg(pointEst[1], choices = c("fullDataEst", "medianBS"))
	
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
	cat(paste0("     N data: ", length(x$Data), "\n"))
	
	if(x$NBootstrap > 0)
		cat(paste0("N bootstrap: ", x$NBootstrap, "\n"))
	
	if(!is.na(x$roundingBase))
		cat(paste0("    rounded: yes (base: ", x$roundingBase, ")\n"))
	else
		cat(paste0("    rounded: no\n"))
	
	if (!is.null(x$AgeMin) & !is.null(x$AgeMax))
		cat(paste0("  Age range: ", x$AgeMin, " to ", x$AgeMax, " years\n"))
	
	if (!is.null(x$Group))
		cat(paste0("     Gender: ", paste(x$Group, collapse=", "), "\n"))		
	
	BSPerformed  <- (length(x$MuBS)>0 & length(x$SigmaBS)>0 & length(x$LambdaBS)>0 & length(x$ShiftBS)>0)
	
	cat(paste0("  point est: ", ifelse(pointEst=="medianBS" & BSPerformed, "medianBS", "fullDataEst"), "\n"))	
	cat(paste0("     lambda: ", ifelse(pointEst=="medianBS" & BSPerformed, signif(x$LambdaMed, 3), signif(x$Lambda, 3)), "\n"))
	cat(paste0("         mu: ", ifelse(pointEst=="medianBS" & BSPerformed, signif(x$MuMed, 3), 	   signif(x$Mu, 3)), 	 "\n"))
	cat(paste0("      sigma: ", ifelse(pointEst=="medianBS" & BSPerformed, signif(x$SigmaMed, 3),  signif(x$Sigma, 3)),  "\n"))
	cat(paste0("      shift: ", ifelse(pointEst=="medianBS" & BSPerformed, signif(x$ShiftMed, 3),  signif(x$Shift, 3)),  "\n"))
	cat(paste0("       cost: ", ifelse(pointEst=="medianBS" & BSPerformed, signif(x$CostMed, 3),   signif(x$Cost, 3)),   "\n"))
	cat(paste0("NP fraction: ", ifelse(pointEst=="medianBS" & BSPerformed, signif(x$PMed, 3), 	   signif(x$P, 3)), 	 "\n"))	
}
