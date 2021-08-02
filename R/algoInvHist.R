#' Generate list with histogram data.
#' 
#' @param x			  (numeric) vector of data points
#' @param ab		  (numeric) vector of lower and higher limit embedding appropriate region with the main peak
#' 
#' @return (list) with histogram data used in the calculation of cost. 
#'
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

generateHistData <- function(x, ab) {
	
	HistData <- list()	
	HistData$counts <- HistData$breakL <- HistData$breakR <- NULL
	
	# calculate number of bins (without overlapping) depending on size of data
	NBins <- min(40, max(9, 11*(sum(x >=ab[1] & x <= ab[2])/5000)^(1/5)))
	
	# calculate number of bins and number of total bins (with overlapping)
	overlapFactor <- round(256/NBins)
	
	NBinsTotal    <- round(NBins)*overlapFactor
	
	# generate table with difference between neighboring values 
	diffVal  <- table(round(diff(sort(unique(x[x>=ab[1] & x<=ab[2]]))), digits=10))				
	# select differences that occur in >= 10% of cases and that have a base 10
	logNames <- log10(as.numeric(names(diffVal)))	
	diffVal  <- diffVal[diffVal>=0.1*sum(diffVal) & logNames==round(logNames)]
	
	roundingBase <- NA
	
	# if dataset has only finite number of unique values (e.g. rounded data)
	if(length(diffVal) > 0)
	{		
		# determine rounding base (maximum step size of discrete values)
		roundingBase <- max(as.numeric(names(diffVal)))	
		
		# determine bin size, which is a multiple of rounding base
		abStepSize <- ceiling(diff(ab)/256/roundingBase)*roundingBase		
		
		# determine overlap factor
		overlapFactor <- max(1, round(diff(ab)/abStepSize/NBins))
		
		# calculate number of unique values to adjust NBins 
		xAB <- unique(x[x>=ab[1] & x<=ab[2]])
		
		if(length(xAB) < NBins){
			NBins <- ceiling(length(xAB)/overlapFactor)
		}
				
		# calculate number of bins and number of total bins (with overlapping)	
		NBinsTotal <- round(NBins)*overlapFactor
		
		# adapt ab to discrete grid
		ab[1] <- max(0.5*roundingBase, round(ab[1]/roundingBase)*roundingBase - 0.5*roundingBase)		
		ab[2] <- ab[1] + (NBinsTotal+overlapFactor-1)*abStepSize			
	}
	
	# define vector of breaks and perform histogramming
	histBreaks <- seq(ab[1], ab[2], length.out = NBinsTotal+overlapFactor)		
	fullHist   <- hist(x[x > histBreaks[1] & x <= histBreaks[NBinsTotal+overlapFactor]], breaks = histBreaks, plot = FALSE)
	
	# combine bins...
	for (i in 1:NBinsTotal) {
		HistData$counts <- c(HistData$counts, sum(fullHist$counts[i:(i+overlapFactor-1)]))			
		HistData$breakL <- c(HistData$breakL, fullHist$breaks[i])
		HistData$breakR <- c(HistData$breakR, fullHist$breaks[i+overlapFactor])	
	}	
	
	# add bin for lower border
	HistData$counts <- c(HistData$counts, sum(x > 1e-20 & x <= ab[1]))
	HistData$breakL <- c(HistData$breakL, max(min(ab[1], x), 1e-20))
	HistData$breakR <- c(HistData$breakR, ab[1])
	
	# add bin for upper border
	HistData$counts <- c(HistData$counts, sum(x > ab[2] & x <= 1e20))
	HistData$breakL <- c(HistData$breakL, ab[2])
	HistData$breakR <- c(HistData$breakR, min(max(ab[2], x), 1e20))
	
	fullHist$counts <- c(sum(x > 1e-20 & x <= ab[1]), fullHist$counts, sum(x > ab[2] & x <= 1e20))
	fullHist$breaks <- c(max(min(ab[1], x), 1e-20), fullHist$breaks, min(max(ab[2], x), 1e20))
	fullHist$mids   <- 0.5*(fullHist$breaks[2:length(fullHist$breaks)] + fullHist$breaks[1:(length(fullHist$breaks)-1)]) 
	
	# calculate mids of histogram bins
	HistData$mids   <- 0.5*(HistData$breakL + HistData$breakR)	
	
	# sort bins in order of midpoints
	sortIndex <- sort(HistData$mids, index.return = TRUE)$ix 
	HistData$counts <- HistData$counts[sortIndex]
	HistData$mids   <- HistData$mids[sortIndex]
	HistData$breakL <- HistData$breakL[sortIndex]
	HistData$breakR <- HistData$breakR[sortIndex]
	
	# add further parameters to list
	HistData$overlapFactor <- overlapFactor	
	HistData$NBins 	 	   <- round(NBins)
	HistData$abOr		   <- ab
	HistData$roundingBase  <- roundingBase
	HistData$NData 		   <- length(x)
	
	HistData$fullHist <- fullHist
	
	return(HistData)
}


#' Function to estimate reference intervals for a single population
#' 
#' The function estimates the optimal parameters lambda, mu and sigma for a raw data set contatining pathological 
#' and non-pathological values. The optimization is carried out via a multi-level grid search to 
#' minimize the cost function (negative log-likelihood with regularization) and to find a model that fits the 
#' distribution of the physiological values and thus separates pathological from non-pathological values.
#'
#' @param Data			(numeric) values specifying data points comprising pathological
#' 						and non-pathological values 
#' @param NBootstrap	(integer) specifying the number of bootstrap repetitions
#' @param seed			(integer) specifying the seed used for bootstrapping
#' 
#' @return (object) of class "RWDRI" with parameters optimized
#'
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}
#' 
#'  
#' @examples
#' 
#' # first example
#' 
#' \donttest{
#' data(testcase1)
#' resRI <- findRI(Data = testcase1)
#' print(resRI)
#' plot(resRI, showPathol = FALSE)
#' }
#'
#' # second example
#' data(testcase2)
#' resRI <- findRI(Data = testcase2)
#' print(resRI)
#' 
#' \donttest{
#' # third example, with bootstrapping 
#' data(testcase3)
#' resRI <- findRI(Data = testcase3, NBootstrap = 30, seed = 123)
#' print(resRI)
#' getRI(resRI, RIperc = c(0.025, 0.5, 0.975), CIprop = 0.95, pointEst ="fullDataEst")
#' getRI(resRI, RIperc = c(0.025, 0.5, 0.975), CIprop = 0.95, pointEst ="medianBS")
#' plot(resRI)
#' # plot without showing values and pathological distribution
#' plot(resRI, showValue = FALSE, showPathol = FALSE)
#' plot(resRI,  RIperc = c(0.025, 0.5, 0.975), showPathol = FALSE, showCI = TRUE) 
#' 
#' # forth example, with bootstrapping
#' data(testcase4)
#' resRI <- findRI(Data = testcase4, NBootstrap = 30)
#' plot(resRI,  RIperc = c(0.025, 0.5, 0.975), showPathol = FALSE, showCI = TRUE)
#' }
#
findRI <- function(Data = NULL,  NBootstrap = 0, seed = 123) {

	stopifnot(!is.null(Data))
	stopifnot(is.numeric(Data))
	stopifnot(is.numeric(NBootstrap) & NBootstrap %%1 == 0 & NBootstrap >= 0)
	stopifnot(is.numeric(seed) & seed%%1 ==0)
	
	if(length(Data) > 50000000)
		warning(paste0("\n Warning: Extremely large dataset (> 50'000'000), might lead to performance and memory complications!\n\n"))
	
	obj	  	   <- list()		
	obj$Data   <- Data
	# remove NAs
	if(anyNA(obj$Data))
		obj$Data <- obj$Data[!is.na(obj$Data)]
		
	# set all obj values for bootstrapping to NA
	obj$LambdaBS <- obj$MuBS <- obj$SigmaBS <- obj$CostBS <- rep(NA, times = NBootstrap) 	
	
	# define iteration indices
	iterations <- 1:(NBootstrap+1)
	
	# apply parallel computing only if > 10 iterations
	if(length(iterations) > 10)
		plan(multisession, workers = min(floor((length(iterations)+0.1)/5), detectCores(logical = TRUE)))
	
	#use a pre-generated sequence of seeds 
	res <- future_lapply(iterations, FUN=function(i) {	
				
				if(i==1)
					Data <- obj$Data
				else 
					Data <- sample(x = obj$Data, size = length(obj$Data), replace = TRUE)
	
				# initial estimation of lambda...
				lambdaVec <- seq(0.0, 1.3, by = 0.1)^1.54542
				
				#estimate start values for mu, sigma, and get estimation of ab for all lambdas
				startVals <- estimateStartValues(Data, lambdaVec)
						
				#generate list with histogram data
				Hist      <- generateHistData(x = Data, ab = startVals$abOr)
								
				bestParam <- c(NA, NA, NA, NA, 1e9, NA, NA )
				bestParam <- testParam(lambdaVec = lambdaVec, bestParam = bestParam, Data = Data, HistData = Hist, 
									   startValues = startVals, NIter = 8, alpha = 0.01, alphaMcb = 0.1)
				
				# estimation with updated lambdaVec..	
				# refine search region of lambda...		
				minIndex  <- min(length(lambdaVec)-2, max(3, which.min(abs(bestParam[1]-lambdaVec))))
			
				lambdaVec <- c(seq(from = lambdaVec[minIndex-2], to = lambdaVec[minIndex-1], length.out = 6)[2:5],
						seq(from = lambdaVec[minIndex-1], to = lambdaVec[minIndex+0], length.out = 6)[2:5],
						seq(from = lambdaVec[minIndex+0], to = lambdaVec[minIndex+1], length.out = 6)[2:5],
						seq(from = lambdaVec[minIndex+1], to = lambdaVec[minIndex+2], length.out = 6)[2:5])			
				
				startValsRefined <- estimateStartValues(Data, lambdaVec)
				bestParam <- testParam(lambdaVec = lambdaVec, bestParam = bestParam, Data = Data, HistData = Hist, startValues = startValsRefined,
									   NIter=8, alpha = 0.01, alphaMcb = 0.1)
					
				return(c(bestParam, Hist$abOr, Hist$roundingBase))		
				
			}, future.seed=future_lapply(1:length(iterations), FUN = function(x) .Random.seed, future.chunk.size = Inf, future.seed = seed), future.packages = c("ash", "refineR"))
			
	if(length(iterations) > 10)
		plan(sequential)
	
	for (i in iterations) {	
		if (i == 1) {
			obj$Lambda       <- res[[i]][1]
			obj$Mu 	         <- res[[i]][2]
			obj$Sigma        <- res[[i]][3]
			obj$P 	         <- res[[i]][4]
			obj$Cost	     <- res[[i]][5]
			obj$abOr         <- res[[i]][6:7]
			obj$roundingBase <- res[[i]][8]
			obj$Method	     <- "refineR"
		} else {
			obj$LambdaBS[i-1] <- res[[i]][1]
			obj$MuBS[i-1] 	  <- res[[i]][2]
			obj$SigmaBS[i-1]  <- res[[i]][3]	
			obj$PBS[i-1]	  <- res[[i]][4]
			obj$CostBS[i-1]   <- res[[i]][5]
		}	
	}	
	
	# remove bootstrap samples with result NA
	selNotNA     <- !is.na(obj$LambdaBS)
	obj$LambdaBS <- obj$LambdaBS[selNotNA]
	obj$MuBS     <- obj$MuBS[selNotNA]
	obj$SigmaBS  <- obj$SigmaBS[selNotNA]
	obj$PBS		 <- obj$PBS[selNotNA]
	obj$CostBS   <- obj$CostBS[selNotNA]
	
	class(obj) <- "RWDRI"
	
	return(obj)
}


#' Helper function to estimate search regions for mu and sigma and to get the region around main peak 'ab'
#' 
#' The function estimates start search regions for mu and sigma for each lambda. Further it determines an appropriate 
#' region around the main peak 'ab' that is used for all lambdas. 
#'
#' @param Data			(numeric) values specifying data points comprising pathological
#' 						and non-pathological values 
#' @param lambdaVec		(numeric) transformation parameter for inverse Box-Cox transformation
#' @param useQuantiles	(logical) indicating if quantiles or raw data should be used (for more than 100000 data points 
#' 							quantiles are always used)
#'  
#' @return (list) with (abOriginal, search region for mu and sigma)
#' 
#' 
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}
#
estimateStartValues <- function(Data, lambdaVec, useQuantiles = FALSE) {
		
	if (useQuantiles | length(Data) > 100000) {
		# calculate 10000 quantiles type=1 for data (to get better performance)
		qData <- as.numeric(quantile(x = Data,probs = seq(0, 1, length.out = 50000), type = 7))	
	} else {
		qData <- Data
	}
	
	muEst 		<- matrix(nrow = length(lambdaVec), ncol = 3)
	sigmaEst 	<- matrix(nrow = length(lambdaVec), ncol = 3)
	
	#determine m for density estimation
	densM <- min(48, max(24, floor(72 - 3/50 * length(unique(qData)))))
	
	c <- 1
	# loop over lambda
	for (lambda in lambdaVec) {
		# BoxCox transformation and density estimation...
		# apply BoxCox transformation to quantiles/data and remove NAs an infinity values
		qDataTr <- suppressWarnings(BoxCox(x = qData, lambda = lambda))			
		qDataTr <- qDataTr[!is.na(qDataTr) & is.finite(qDataTr)]		
		
		#get range of transformed data for density estimation
		abTr 	 <- range(qDataTr)	
		abTrDiff <- diff(abTr)
		abTr[1]  <- abTr[1]-0.1*abTrDiff
		abTr[2]  <- abTr[2]+0.1*abTrDiff
		
		qDensTr <- ashDensity(x = qDataTr, ab = abTr, nbin = 512, m = densM)
		
		# find peak with largest area
		listPV <- findPeaksAndValleys(qDensTr)
		peaks 	<- listPV$peaks
		valleys <- listPV$valleys
		area 	<- vector(mode = "list")
		
		#get area under peak (between two valleys) 
		for (counter in 1:length(peaks)) {
			if (counter ==1) {
				vLow <- 1
			} else {
				vLow <- valleys[counter-1] #-1 as the function returns the value after the peak/valley
			}
			p <- peaks[counter]
			
			if (counter > length(valleys)) {
				vHigh <- length(qDensTr$y)
			} else {
				vHigh <- valleys[counter]
			}
			area[[as.character(p)]] <- sum(qDensTr$y[vLow:vHigh])
		}		
		
		#get index of peak in area with biggest surface
		peakInd <- as.numeric(names(which.max(unlist(area))))
		modEst 	<- sum(qDensTr$x[peakInd+c(-1:1)]*qDensTr$y[peakInd+c(-1:1)])/sum(qDensTr$y[peakInd+c(-1:1)])	
		
		# select region around peak for stronger smoothing
		abPeak	  <- NULL
		abPeak[1] <- max(qDensTr$x[which(qDensTr$x < modEst & qDensTr$y <= 0.35*qDensTr$y[peakInd])])
		abPeak[2] <- min(qDensTr$x[which(qDensTr$x > modEst & qDensTr$y <= 0.35*qDensTr$y[peakInd])])
		
		# apply more significant/stronger smoothing in selected region
		densMM <- min(48, max(16, floor(48 - 1/25 * length(unique(qDataTr[qDataTr >= abPeak[1] & qDataTr <= abPeak[2]])))))		
		qDensTr <- ashDensity(x = qDataTr, ab = abPeak, nbin = 256, m = densMM)
				
		# find peak with largest area
		listPV <- findPeaksAndValleys(qDensTr)
		peaks 	<- listPV$peaks
		valleys <- listPV$valleys
		area 	<- vector(mode = "list")
	
		#get area under peak (between two valleys) 
		for (counter in 1:length(peaks)) {
			if (counter ==1) {
				vLow <- 1
			} else {
				vLow <- valleys[counter-1] #-1 as the function returns the value after the peak/valley
			}
			p <- peaks[counter]
			
			if (counter > length(valleys)) {
				vHigh <- length(qDensTr$y)
			} else {
				vHigh <- valleys[counter]
			}
			area[[as.character(p)]] <- sum(qDensTr$y[vLow:vHigh])
		}		
		
		peakInd <- as.numeric(names(which.max(unlist(area))))
		modEst 	<- sum(qDensTr$x[peakInd+c(-1:1)]*qDensTr$y[peakInd+c(-1:1)])/sum(qDensTr$y[peakInd+c(-1:1)])	

		# find width of peak for densities 0.5 to 0.9 of max(density) in 0.01 steps
		DensRange <- seq(0.5, 0.95, by = 0.05)*qDensTr$y[peakInd]
		
		widthsL <- NULL
		widthsR <- NULL
		
		for (d in DensRange) {
			# get width (distance) from peak to the left to the first points that is smaller or equal to d
			widthsL <- c(widthsL, modEst - max(qDensTr$x[which(qDensTr$x < modEst & qDensTr$y <= d)]))
			# get width (distance) from peak to the right to the first point that is smaller or equal to d
			widthsR <- c(widthsR, (min(qDensTr$x[which(qDensTr$x > modEst & qDensTr$y <= d)]) - modEst))
			
		}

		# get points left and right of the peak 
		left 	<- modEst - widthsL
		right 	<- modEst + widthsR
		
		# leftRight Matrix: cols: 1=density range, 2=left value, 3=right value, 4=widthsL, 5=widthsR, 6=mean of left and right --> mu
		leftRight 	<- matrix(c(seq(0.5, 0.95, by = 0.05), left, right, widthsL, widthsR), nrow = length(left), ncol = 5, byrow = FALSE)
		leftRight 	<- cbind(leftRight, apply(leftRight[, 2:3], 1, mean))
						
		muTmp 	<- c(leftRight[, 6], modEst)
			
		# estimate range of muEst as range of centers of c(widths, modEst)				
		muRange 	<- range(muTmp)
		#add margins to muRange	
		diffMuRange	<- diff(muRange)
		muEstRange	<- c(muRange[1] - 0.5*diffMuRange, muRange[2] + 0.5*diffMuRange)		
		muEst[c,] 	<- c(lambda, muEstRange)		
						
		# estimate range of sigmaEst as range of sds calculated from widths
		# minimal width (as width is highly influenced by pathological samples), denominator is derived from HWHM formula
		leftRight <- cbind(leftRight, apply(leftRight[, c(1, 4:5)], 1, function(x) {
							return(min(x[c(2, 3)])/(sqrt((-2)*log(x[1]))))
						}))
		sigmaTmp  <- leftRight[, 7]
						
		# add margins to sigmaRange on both sides
		sigmaRange 		<- range(sigmaTmp)
		diffSigmaRange 	<- diff(sigmaRange)
		sigmaEstRange	<- c(max(0, sigmaRange[1] - 2.5*diffSigmaRange), sigmaRange[2]+1*diffSigmaRange)	
		sigmaEst[c,] 	<- c(lambda, sigmaEstRange)
	
		c = c +1
	}
	
	# select lambda with lowest normalized ranges of muEst and sigmaEst => range(muEst)/mean(muEst) and range(sigmaEst)/mean(sigmaEst)
	normRangeMu 	<- abs((muEst[, 3] - muEst[, 2])/rowMeans(muEst[, 2:3]))
	normRangeSigma 	<- (sigmaEst[, 3] - sigmaEst[, 2])/rowMeans(sigmaEst[, 2:3])
	
	normDf <- matrix(c(normRangeMu, normRangeSigma), ncol = 2, nrow = length(lambdaVec), byrow = FALSE)
	normDf <- cbind(normDf, sqrt(normDf[, 1]^2 + normDf[, 2]^2))
	
	lambdaInd 	<- which.min(normDf[, 3])
	lambda 		<- lambdaVec[lambdaInd]
	
	# for that lambda, calculate mu +/- 3 sigma and transform to original space
	abTr 	<- NULL
	abTr[1] <- mean(muEst[lambdaInd, 2:3]) - 3*mean(sigmaEst[lambdaInd, 2:3])
	abTr[2] <- mean(muEst[lambdaInd, 2:3]) + 3*mean(sigmaEst[lambdaInd, 2:3])
	
	ab 	 <- invBoxCox(abTr, lambda)
		
	# add small margins to get ab, which is used for all lambda values in the optimization
	if (is.na(ab[1])) ab[1] <- 1e-20
	ab[1] <- ab[1] - 0.05*diff(ab)
	ab[2] <- ab[2] + 0.05*diff(ab)
		
	# improve determination of ab by considering the min/max of the data, the measuring range, and ensuring a minimum density 
	abDens 	 	<- range(qData)
	abDensDiff 	<- diff(abDens)
	abDens[1] 	<- abDens[1] - 0.1*abDensDiff
	abDens[2] 	<- abDens[2] + 0.1*abDensDiff
	qDensOr 	<- ashDensity(x = qData, ab = abDens, nbin = 512, m = 48)	
	muEstLambda <- invBoxCox(max(muEst[lambdaInd, 2:3]), lambda)
	
	# get points where density is below 5% or above 95%
	abDensL  <- qDensOr$x[which(qDensOr$x < muEstLambda & qDensOr$y < 0.05*max(qDensOr$y))]
	abDensR  <- qDensOr$x[which(qDensOr$x > muEstLambda & qDensOr$y < 0.05*max(qDensOr$y))]		

	# determine ab 
	perc85 <-  as.numeric(quantile(x = Data, probs = 0.85, na.rm = TRUE))
	abOr <- NULL
	abOr[1] <- max(min(Data), ab[1], abDensL, 1e-20, na.rm = TRUE)
	abOr[2] <- max(perc85, min(max(Data), ab[2], abDensR, na.rm = TRUE))
	
	return(list("abOr" = abOr, "muEst" = muEst, "sigmaEst" = sigmaEst))
	
}


#' Helper function to find optimal parameters lambda, mu and sigma.
#' 
#' @param lambdaVec		(numeric) transformation parameter for inverse Box-Cox transformation
#' @param bestParam		(numeric) vector containing best guess for lambda, mu, sigma, P, cost
#' @param Data			(numeric) values specifying percentiles or data points comprising pathological
#' 						and non-pathological values 
#' @param HistData		(list) with histogram data
#' @param startValues	(list) with start search regions for mu and sigma
#' @param NIter			(integer) specifying the number of iterations for optimized grid-search
#' @param alpha			(numeric) specifying the confidence region used for selection of histogram bins in cost calculation
#' @param alphaMcb		(numeric) specifying the confidence level defining the maximal allowed counts below the asymmetric confidence region
#' 
#' @return (numeric) vector with best parameters for lambda, mu, sigma, P, cost. 
#' 
#'
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

testParam <- function(lambdaVec, bestParam, Data, HistData, startValues, NIter, alpha = 0.01, alphaMcb = 0.1) {		
	
	h <- 0.002
	pNormLookup <- list(h = h, val = pnorm(q = seq(from = -5.2, to = 5.2+h, by = h)))
	
	c <- 1
	# loop over all lambda values
	for (lambda in lambdaVec) {
		# generate mu and sigma vectors 
		muVec <- startValues$muEst[c, 2:3]
		muVec <- seq(muVec[1], muVec[2], length.out = 9)
		
		sigmaVec <- startValues$sigmaEst[c, 2:3]
		sigmaVec <- seq(sigmaVec[1], sigmaVec[2], length.out = 13)
					
		muUnique 		<- sigmaUnique <- NULL			
		currentBestMu   <- currentBestSigma <- NA	
		currentBestCost <- 1e9
				
		# loop for optimization using adapted grid search
		iter = 1
		continueIteration = TRUE
		while (continueIteration & iter <= NIter) {				
			#testedParam <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 1e9)
			testedParam <- c(NA, NA, NA, NA, 1e9)
			muUnique    <- sort(unique(c(muUnique,    muVec)),    method="quick")
			sigmaUnique <- sort(unique(c(sigmaUnique, sigmaVec)), method="quick")			
			
			if (iter > 1) {	
				# if in second iteration mu and sigma are still NA, no value in the grid will fit -> stop iteration
				# else get optimized grid 
				if (is.na(currentBestMu) & is.na(currentBestSigma) & iter == 2) {
					continueIteration = FALSE
				} else {
					muVec 	 <- optimizeGrid(currentBestMu, muUnique, iter, sigmLimit = FALSE)
					sigmaVec <- optimizeGrid(currentBestSigma, sigmaUnique, iter, sigmLimit = TRUE)
				}
			}	
			
			if ((iter == 1) | (iter > 1 & !is.na(currentBestMu) & !is.na(currentBestSigma))) {
				testedParam <- calculateCostHist(lambda = lambda, muVec = muVec, sigmaVec = sigmaVec, HistData = HistData, alpha = alpha, alphaMcb = alphaMcb, pNormLookup = pNormLookup)				
			}			
			
			# update currentBest Parameters is better parameters (smaller costs) are found
			if (testedParam[5] < currentBestCost) {
				currentBestMu    <- testedParam[2]
				currentBestSigma <- testedParam[3]
				currentBestCost  <- testedParam[5]
			}	
			
			# update overall bestParam if tested parameters are superior
			if (testedParam[5] < bestParam[5]) {
				bestParam    <- testedParam
				
			}	
			iter = iter + 1
		}	
		c = c + 1			
	}	
	
	return(bestParam)
}


#' Calculate costs for a specific combinations of lambda, muVec and sigmaVec.
#' 
#' @param lambda		(numeric) transformation parameter for inverse Box-Cox transformation
#' @param muVec			(numeric) vector of mean values of non-pathological Gaussian distribution in transformed space
#' @param sigmaVec		(numeric) vector of sd values of non-pathological Gaussian distribution in transformed space
#' @param HistData		(list) with histogram data generated by function \code{\link{generateHistData}}
#' @param alpha			(numeric) specifying the confidence region used for selection of histgram bins in cost calculation
#' @param alphaMcb	 	(numeric) specifying the confidence level defining the maximal allowed counts below asymmetric confidence region
#' @param pNormLookup 	(list) with lookup table for pnormApprox function \code{\link{pnormApprox}}
#' 
#' @return (numeric) vector with (lambda, mu, sigma, P, cost). 
#'
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

calculateCostHist <- function(lambda, muVec, sigmaVec, HistData, alpha = 0.01, alphaMcb = 0.1, pNormLookup) {	

	N <- HistData$NData
	countsData 	   <- HistData$counts
	countsFullHist <- HistData$fullHist$counts 
	
	# Box Cox transformation of histogram breaks and histogram range
	breakL <- suppressWarnings(BoxCox(HistData$breakL, lambda))		
	breakR <- suppressWarnings(BoxCox(HistData$breakR, lambda))		
	
	# parameters for pnormApprox function
	oneOverH <- 1/pNormLookup$h
	pNormVal <- pNormLookup$val
	
	# limits for surrounding peak areas needed fo calculation of P 
	pLimitMin <- c(0.5, 0.55, 0.6, 0.65, 0.7)
	pLimitMax <- c(0.9, 0.94, 0.95, 0.96, 0.97) 

	# precalculated factor for confidence region
	qNormFactor <- -qnorm(p = alpha/2)
	
	# maximum number of bins expected to be below confidence region
	maxCountsBelow <- as.integer(ceiling(length(countsData)*alphaMcb/2))
	# initial bestParam
	bestParam <- c(NA, NA, NA, NA, 1e9)
	
	for (mu in muVec) {
		for (sigma in sigmaVec)	{	
			# theoretical prediction of bin counts
			countsPred <- HistData$NData*(pnormApprox(q = breakR, mean = mu, oneOverSd = 1/sigma, oneOverH = oneOverH, pNormVal = pNormVal) - pnormApprox(q = breakL, mean = mu, oneOverSd = 1/sigma, oneOverH = oneOverH, pNormVal = pNormVal))
			
			countsPred[countsPred < 0] <- 0
			
			# selection of bins that can be approximated by normal distribution
			selectionCounts <- (countsPred >= 40)
			# selection of peak bins
			selectionPeak   <- (countsPred >= 0.95*max(countsPred))				
			
			# lower 0.5% should not be negative, thus check if lowerBound is NA
			RRLowerBound <- invBoxCox(mu-qNormFactor*sigma, lambda = lambda)	
			
			if (!is.na(sum(selectionCounts)) & sum(selectionCounts) > (HistData$NBins*HistData$overlapFactor)/16 &  !is.na(RRLowerBound)) {					
					
				# calculation of slope serving as estimation of P
				cData <- countsData[selectionCounts]
				cPred <- countsPred[selectionCounts]		
											
				# get index of peak 
				indexPeak <- which.max(countsPred)
			
				# get sum of data points in specified region next to the peak for the actual data and the prediction
				peakArea  <- getSumForPArea(pLimitMin = pLimitMin, pLimitMax = pLimitMax, countsPred = countsPred, indexPeak = indexPeak, HistData = HistData, countsFullHist = countsFullHist, lambda = lambda, mu = mu, sigma = sigma)
				sumData11 <- peakArea$sumData11
				sumPred11 <- peakArea$sumPred11	

				ratio <- sumData11/sumPred11
				
				#increase ratio until all exceed data points 
				while (all(ratio <= 1.0) & any(ratio*sumPred11-qNormFactor*sqrt(ratio*sumPred11) <= sumData11)) {
					rInd <- ratio*sumPred11-qNormFactor*sqrt(ratio*sumPred11) <= sumData11
					ratio[rInd] <- ratio[rInd]+0.001
				}	
				
				PMax <- min(max(0.401, min(ratio)), 1.0)
				
				ratio <- sumData11/sumPred11
				#decrease ratio until all fall below the observed data point
				while(all(ratio >= 0.4) & any(ratio*sumPred11+qNormFactor*sqrt(ratio*sumPred11) >= sumData11)){ 
					rInd <- ratio*sumPred11+qNormFactor*sqrt(ratio*sumPred11) >= sumData11
					ratio[rInd] <- ratio[rInd] - 0.001
				}
				
				P <- max(0.4, min(ratio))
				
				# PStep is bounded between 0.0025 and 0.005 to ensure an appropriate increment step
				diffP <- PMax - P 
				PStep <- min(max(0.0025, diffP/8), 0.005)
				
				# adapt lower P such that PMax can be reached using calculated step size
				P <- min(PMax, 1.0-1e-20-ceiling((1.0-1e-20-P)/PStep)*PStep)

				continuePLoop <- TRUE
								
				# test different values for P and calculate costs...
				while (continuePLoop) {					
					# calculate actual predicted counts and corresponding confidence region...
					countsPredP      <- countsPred*P					
					confWidth 		 <- qNormFactor * sqrt(countsPredP)
					countsPredPLower <- countsPredP - confWidth
					countsPredPUpper <- countsPredP + P*confWidth
										
					# select bin within confidence region
					selectionEval    <- (selectionPeak | (selectionCounts & countsData <= countsPredPUpper))					
					sumSelectionEval <- sum(selectionEval)	
					
					# decide if loop shall be continued...
					continuePLoop <- (sum(countsData < countsPredPLower) <= maxCountsBelow & P <= PMax)
					# if expected number of bins are selected...
					if (sumSelectionEval > (HistData$NBins*HistData$overlapFactor)/16 & continuePLoop) {										
						# extract data of selection
						cData <- countsData[selectionEval]										
						cPred <- countsPredP[selectionEval]	
											
						# calculate negative log likelihood with regularization term
						pois 	  <- dpois(x = round(cData), lambda = cPred)*cPred
						nrRemove  <- length(pois[pois == 0])
						pois 	  <- pois[pois > 0]
						cost 	  <- -sum(log(pois))/sqrt(sumSelectionEval-nrRemove) 		#+1 as cData has an additional element	
								
						# update best cost
						if (!is.na(cost) & cost < bestParam[5]) {						
							bestParam <- c(lambda, mu, sigma, P, cost)													
						}						
					}	
					P = P + PStep
				}
			}
		}		
	}
	
	bestParam
}

#' Find the index of the peaks and valleys of the density estimation.
#' 
#' @param Dens		(list) with density estimation (x values, y values)
#' 
#' @return (list) specifying the index of the peaks and valleys of the density estimation. 
#'
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

findPeaksAndValleys <- function(Dens) {
	
	diffDens <- diff(Dens$y)
	peakCrit <- valleyCrit <- 0
	peaksV   <- valleysV  <- NULL 
	prev <- act <- diffDens[1]
	prevSgn <- ifelse(act < 0, "-", "+")
	
	lastP = FALSE
	lastV = FALSE

	for (i in 2:length(diffDens)) {
		# get current value and sign
		act <- diffDens[i]
		sgn <- ifelse(act < 0, "-", "+")
		
		# determine local minima and maxima depending on the change of signs (- to + is minima; + to - is maxima)
		if (prevSgn == '-' & sgn == '+') {
			valleyCrit <- 1
			prevSgn   <- sgn
			
		} else if(sgn == '+' & valleyCrit > 0) {
			valleyCrit <- valleyCrit + 1
			
		} else if (prevSgn == '+' & sgn == '-') {
			peakCrit <-  1
			prevSgn  <- sgn
			
		} else if (sgn == '-' & peakCrit > 0) {
			peakCrit <- peakCrit +1
			
		} else {
			prevSgn		<- sgn 
			prev		<- act
			peakCrit 	<- 0
			valleyCrit 	<- 0
		}
		#after a minima or maxima, the same sign has to appear at least 5 times, otherwise it is considered noise
		if (peakCrit == 5) {
			if (!lastP) {
				maxPeak <- which.max(Dens$y[(i-5):i])
				peaksV  <- c(peaksV, (i-5+maxPeak-1))
				lastP   <- TRUE
				lastV   <- FALSE
			}
			peakCrit 	<- 0
			valleyCrit 	<- 0 
			prevSgn 	<- sgn 
			prev 		<- act
			
		}
		if (valleyCrit == 5) {
			if (!lastV) {
				minValley <- which.min(Dens$y[(i-5):i])
				valleysV  <- c(valleysV, (i-5+minValley-1))
				lastP     <- FALSE
				lastV     <- TRUE
			}
			
			valleyCrit <- 0
			peakCrit   <- 0
			prevSgn    <- sgn 
			prev 	   <- act
		}
	}
		
	# if distribution starts with a valley, delete it
	if (length(valleysV)!= 0 && valleysV[1] < peaksV[1]) 
		valleysV <- valleysV[-1]
	
	# if there are no valleys found due to smoothed density function, add a valley at the end
	if(length(valleysV)==0)
		valleysV <- length(Dens$x)
	
	return (list("peaks" = peaksV, "valleys"= valleysV))
}


#' Helper function for grid search for mu and sigma.
#' 
#' @param currentBestParam		(numeric) value specifying the current best value for this parameter
#' @param paramUnique			(numeric) vector of possible values for this parameter
#' @param iter 					(integer) indicating the number of iteration, as in the first iteration the search region
#' 									is larger than in the following iterations 
#' @param sigmLimit 			(logical) specifiying if parameter is sigma and thus minimum is 0
#'  
#' @return (vector) specifying the new search region fo the parameter to be optimized
#'
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}
#' 
optimizeGrid <- function(currentBestParam, paramUnique, iter, sigmLimit = TRUE) {	

	optInd <- which.min(abs(paramUnique - currentBestParam))
	
	if (!is.na(currentBestParam) && optInd == 1) {
		# param is at left border, adapted grid exceeds original search space
		diff <- paramUnique[2]-paramUnique[1]
		# in first iteration, adapted grid is larger than in following iterations
		if (iter == 2) {
			# account for lower limit of sigma, cannot be < 0
			if (sigmLimit) 
				tmp <- c(max(paramUnique[optInd]-2*diff, 0), max(paramUnique[optInd]-diff, 0), paramUnique)
			else
				tmp <- c(paramUnique[optInd]-2*diff, paramUnique[optInd] - diff, paramUnique)
			optInd = optInd+2
		} else {
			# account for lower limit of sigma, cannot be < 0
			if (sigmLimit)
				tmp <- c(max(paramUnique[optInd]-diff, 0), paramUnique)
			else 
				tmp <- c(paramUnique[optInd]-diff, paramUnique)
			optInd = optInd + 1
		}
		
		OptL     <- max(optInd - ifelse(iter == 2, 2, 1), 1)
		OptR 	 <- min(optInd + ifelse(iter == 2, 2, 1), length(tmp))
		paramVec <- c(currentBestParam, 0.5*(tmp[OptL] + tmp[OptL+1]), 0.5*(tmp[OptR]+tmp[OptR-1]))
		
	} else if (!is.na(currentBestParam) && optInd == length(paramUnique)) {
		#param is at right border /upper margin of vector, adapted grid exceeds original search space
		diff <- paramUnique[optInd] - paramUnique[optInd-1]
		# in first iteration, adapted grid is larger than in following iterations
		if (iter == 2) {
			tmp <- c(paramUnique, paramUnique[optInd]+diff, paramUnique[optInd]+2*diff)
		} else {
			tmp <- c(paramUnique, paramUnique[optInd]+diff)
		}
		
		OptL 	 <- max(optInd - ifelse(iter == 2, 2, 1), 1)
		OptR 	 <- min(optInd + ifelse(iter == 2, 2, 1), length(tmp))
		paramVec <- c(currentBestParam, 0.5*(tmp[OptL]+ tmp[OptL+1]), 0.5*(tmp[OptR]+tmp[OptR-1]))
		
	} else if (!is.na(currentBestParam) && optInd == 2 & iter == 2) {
		#param is at position 2 but iter is 2, adapted grid exceeds original search space
		diff <- paramUnique[2]- paramUnique[1]
		# account for lower limit of sigma, cannot be < 0
		if (sigmLimit)
			tmp <- c(max(paramUnique[optInd]-2*diff, 0), paramUnique)
		else 
			tmp <- c(paramUnique[optInd]-2*diff, paramUnique)
		
		optInd = optInd+1
		
		OptL 	 <- max(optInd-2, 1)
		OptR 	 <- min(optInd+2, length(tmp))
		paramVec <- c(currentBestParam, 0.5*(tmp[OptL] + tmp[OptL+1]), 0.5*(tmp[OptR]+tmp[OptR-1]))
		
	} else if (!is.na(currentBestParam) && optInd == length(paramUnique)-1 & iter == 2) {
		#param is at position length -1 (penultimate position) and iter is 2, adapted grid exceeds original search space
		diff <- paramUnique[optInd]-paramUnique[optInd-1]
		tmp  <- c(paramUnique, paramUnique[optInd]+2*diff)
		OptL <- max(optInd-2, 1)
		OptR <- min(optInd+2, length(tmp))
		paramVec <- c(currentBestParam, 0.5*(tmp[OptL]+tmp[OptL+1]), 0.5*(tmp[OptR]+tmp[OptR-1]))
		
	} else {
		#param is in the middle of the vector/grid
		OptL <- max(optInd - ifelse(iter == 2, 2, 1), 1)
		OptR <- min(optInd + ifelse(iter == 2, 2, 1), length(paramUnique))					
		paramVec <- c(currentBestParam, 0.5*(paramUnique[OptL] + paramUnique[OptL+1]), 0.5*(paramUnique[OptR] + paramUnique[OptR-1]))
	}
	
	return(paramVec)
	
}


#' Helper function to calculate the amount of observed and estimated data points within specified regions around the peak.
#' 
#' The function helps to define the search region for P (fraction of non-pathological samples). 
#' 
#' @param pLimitMin			(numeric) vector specifying the lower limits for the regions next to the peak
#' @param pLimitMax			(numeric) vector specifying the upper limits for the regions next to the peak
#' @param countsPred		(numeric) vector with the predicted counts
#' @param indexPeak 		(integer) specifying the position of the peak
#' @param HistData			(list) with histogram data generated by function \code{\link{generateHistData}}
#' @param countsFullHist 	(integer) vector with the observed counts
#' @param lambda			(numeric) transformation parameter for inverse Box-Cox transformation
#' @param mu				(numeric) parameter of the mean of non-pathological distribution 
#' @param sigma				(numeric) parameter of the standard deviation of non-pathological distribution
#'  
#' @return (list) with two numeric vectors specifying the amount of observed and estimated data points surrounding the peak 
#'
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

getSumForPArea <- function(pLimitMin, pLimitMax, countsPred, indexPeak, HistData, countsFullHist, lambda, mu, sigma) {
	
	N <- length(countsPred)
	
	breaks    <- HistData$fullHist$breaks
	nr_limits <- length(pLimitMin)
	
	sumData11 <- rep(1, times = 2*nr_limits+1)
	sumPred11 <- rep(1, times = 2*nr_limits+1)
	
	borderL_vec <- rep(0, times = 2*nr_limits+1)
	borderR_vec <- rep(0, times = 2*nr_limits+1)
	
	peak <- countsPred[indexPeak]
	
	# selection of peak bins
	selectionPeak <- (countsPred >= 0.95*peak)				
	
	#get leftmost and rightmost value of peak bins
	borderL <- min(HistData$breakL[selectionPeak])
	borderR <- max(HistData$breakR[selectionPeak])
	borderL_vec[1] <- borderL
	borderR_vec[1] <- borderR
	
	#get observed counts at the peak
	sumData11[1] <- sum(countsFullHist[breaks >= borderL & breaks<borderR])
	
	if (length(indexPeak)>1)
		indexPeak <- indexPeak[1]
	
	for (i in 1:length(pLimitMin)) {
		
		#get index of area left from the peak
		selectionLeft <- (countsPred[1:indexPeak] >= pLimitMin[i]*peak & countsPred[1:indexPeak] < pLimitMax[i]*peak)
		
		if (sum(selectionLeft) > 0) {
			#get leftmost and rightmost value of selected region
			borderL <- min(HistData$breakL[1:indexPeak][selectionLeft])
			borderR <- max(HistData$breakR[1:indexPeak][selectionLeft])
			borderL_vec[2*i] <- borderL
			borderR_vec[2*i] <- borderR
			
			#get observed counts in the region next to the peak
			sumData11[2*i] <- sum(countsFullHist[breaks >= borderL & breaks<borderR])
		}	
		
		#get index of area right from the peak
		selectionRight <- (countsPred[indexPeak:N] >= pLimitMin[i]*peak & countsPred[indexPeak:N] < pLimitMax[i]*peak)
		
		if (sum(selectionRight) > 0) {
			#get leftmost and rightmost value of selected region
			borderL <- min(HistData$breakL[indexPeak:N][selectionRight])
			borderR <- max(HistData$breakR[indexPeak:N][selectionRight])
			borderL_vec[2*i+1] <- borderL
			borderR_vec[2*i+1] <- borderR
			
			#get observed counts in the region next to the peak
			sumData11[2*i+1] <- sum(countsFullHist[breaks >= borderL & breaks<borderR])
		}			
	}
	
	#get predicted counts at peak and in the regions next to the peak
	sumPred11 <- HistData$NData*(pnorm(q=BoxCox(borderR_vec, lambda), mean = mu, sd = sigma) - pnorm(q = BoxCox(borderL_vec, lambda), mean = mu, sd = sigma))
	sumPred11[sumPred11 == 0.0] <- 1
	
	return(list(sumData11 = sumData11, sumPred11 = sumPred11))	
}