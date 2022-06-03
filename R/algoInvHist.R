#' Generate list with histogram data.
#' 
#' @param x			   (numeric) vector of data points
#' @param ab		   (numeric) vector of lower and higher limit embedding appropriate region with the main peak
#' @param roundingBase (numeric) describing the rounding base of the dataset
#' 
#' @return (list) with histogram data used in the calculation of cost. 
#'
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

generateHistData <- function(x, ab, roundingBase) {
		
	HistData <- list()	
	HistData$counts <- HistData$breakL <- HistData$breakR <- NULL
	
	# calculate number of bins (without overlapping) depending on size of data
	NBins <- min(40, max(9, 11*(sum(x >=ab[1] & x <= ab[2])/5000)^(1/5)))
	
	# calculate number of bins and number of total bins (with overlapping)
	overlapFactor <- round(256/NBins)
	
	NBinsTotal    <- round(NBins)*overlapFactor
		
	# if dataset has only finite number of unique values (e.g. rounded data)
	if(!is.na(roundingBase))
	{	
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
	HistData$counts <- c(HistData$counts, sum(x > -1e-20 & x <= ab[1]))
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
	HistData$NData 		   <- length(x)
	
	HistData$fullHist <- fullHist
		
	return(HistData)
}


#' Function to estimate reference intervals for a single population
#' 
#' The function estimates the optimal parameters lambda, mu and sigma for a raw data set containing pathological 
#' and non-pathological values. The optimization is carried out via a multi-level grid search to 
#' minimize the cost function (negative log-likelihood with regularization) and to find a model that fits the 
#' distribution of the physiological values and thus separates pathological from non-pathological values.
#'
#' @param Data			(numeric) values specifying data points comprising pathological
#' 						and non-pathological values 
#' @param model			(character) specifying the applied model (can be either "BoxCox" (default), "modBoxCoxFast" or "modBoxCox"),
#' 							option "modBoxCoxFast" and "modBoxCox" first runs the original optimization using the Box-Cox transformation, 
#' 							afterwards the modified Box-Cox transformation is utilized and an optimal shift is identified 
#' 							('fast': only 1 iteration is carried out to find a shift)
#' @param NBootstrap	(integer) specifying the number of bootstrap repetitions
#' @param seed			(integer) specifying the seed used for bootstrapping
#' @param ...			additional arguments to be passed to the method
#' 
#' @return (object) of class "RWDRI" with parameters optimized
#'
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}
#' 
#'  
#' @examples
#' 
#' # first example
#' \donttest{
#' data(testcase1)
#' resRI <- findRI(Data = testcase1)
#' print(resRI)
#' plot(resRI, showPathol = FALSE)
#' 
#' # second example
#' data(testcase2)
#' resRI <- findRI(Data = testcase2)
#' print(resRI, RIperc = c(0.025, 0.5, 0.975))
#' plot(resRI, showPathol = FALSE)
#' 
#' # third example, with bootstrapping 
#' data(testcase3)
#' resRI <- findRI(Data = testcase3, NBootstrap = 30, seed = 123)
#' print(resRI)
#' getRI(resRI, RIperc = c(0.025, 0.5, 0.975), CIprop = 0.95, pointEst ="fullDataEst")
#' getRI(resRI, RIperc = c(0.025, 0.5, 0.975), CIprop = 0.95, pointEst ="medianBS")
#' plot(resRI)
#' 
#' # forth example, without values and pathological distribution in plot function 
#' data(testcase4)
#' resRI <- findRI(Data = testcase4)
#' print(resRI)
#' plot(resRI, showValue = FALSE, showPathol =FALSE) 
#' 
#' # fifth example, with bootstrapping
#' data(testcase5)
#' resRI <- findRI(Data = testcase5, NBootstrap = 30)
#' plot(resRI,  RIperc = c(0.025, 0.5, 0.975), showPathol = FALSE, showCI = TRUE)
#' }

findRI <- function(Data = NULL, model = c("BoxCox", "modBoxCoxFast", "modBoxCox"), NBootstrap = 0, seed = 123,...) {

	args = list(...)
	
	# check of input parameters 
	stopifnot(!is.null(Data))
	stopifnot(is.numeric(Data))
	stopifnot(length(Data) > 10) 
	stopifnot(is.numeric(NBootstrap) & NBootstrap %%1 == 0 & NBootstrap >= 0)
	stopifnot(is.numeric(seed) & seed%%1 ==0)
	stopifnot(is.character(model))
	stopifnot(names(args) %in% c("NCores"))	# stop if args contains something else than NCores 
	
	model <- match.arg(model[1], choices = c("BoxCox", "modBoxCoxFast", "modBoxCox"))
	
	NCores <- args$NCores
	stopifnot(is.null(NCores) | is.numeric(NCores))
	
	if(length(Data) > 50000000)
		cat(paste0("\n Warning: Extremely large dataset (> 50'000'000), might lead to performance and memory complications!\n\n"))
		
	obj	  	   <- list()		
	obj$Data   <- Data
	# remove NAs
	if(anyNA(obj$Data))
		obj$Data <- obj$Data[!is.na(obj$Data)]	
	
	# remove negative values 
	obj$Data <- obj$Data[obj$Data >= 0] 
		
	# set all obj values for bootstrapping to NA
	obj$LambdaBS <- obj$MuBS <- obj$SigmaBS <- obj$CostBS <- rep(NA, times = NBootstrap) 	
	
	# define iteration indices
	iterations <- 1:(NBootstrap+1)

	# apply parallel computing only if > 10 iterations
	if(length(iterations) > 10)
		plan(multisession, workers = max(1, min(floor((length(iterations)+0.1)/5), detectCores(logical = TRUE), NCores)))

	NIterMBC <- 1
	
	if(model=="modBoxCoxFast")
		NIterMBC <- 2
	else if(model=="modBoxCox")
		NIterMBC <- 6
	
	#use a pre-generated sequence of seeds 
	res <- future_lapply(iterations, FUN=function(i) {	
				
			if(i==1)
				Data <- obj$Data
			else 
				Data <- sample(x = obj$Data, size = length(obj$Data), replace = TRUE)			
					
			bestParamShifted <- c(NA, NA, NA, NA, 1e9, NA)
						
			# initialize parameters
			bestCost  <- 1e9
			bestShift <- 0
			shift     <- 0
			addShift  <- 0
			  
			# loop to determine optimal shift (first iteration has shift=0)
			for(iterMBC in 1:NIterMBC)
			{					
				if(iterMBC == 1 | addShift > 1e-20)
				{				
					# iteratively adapt shift and apply to data
					shift 		<- shift + max(0.6, 1.40-iterMBC*0.20)*addShift									
					DataShifted <- Data - shift
					   
					# estimate rounding base
					roundingBase <- findRoundingBase(x = DataShifted)
					   				   
					# initial estimation of lambda...
					lambdaVec <- (1/12*(0:15))^1.8170595
									
					#define search regions for mu, sigma, and get estimation of ab for all lambdas determining the area around the main peak
					startVals <- defineSearchRegions(x = DataShifted, lambdaVec = lambdaVec, roundingBase = roundingBase, abEst = NULL)
						
					#generate list with histogram data
					Hist      <- generateHistData(x = DataShifted, ab = startVals$abEst$abHist, roundingBase = roundingBase)
						
					if(iterMBC == 1)					
						abOr <- Hist$abOr				
																
					bestParam <- c(NA, NA, NA, NA, 1e9, NA)
								   
					bestParam <- testParam(lambdaVec = lambdaVec, bestParam = bestParam, Data = DataShifted, HistData = Hist, 
										   startValues = startVals, NIter = 8, alpha = 0.01, alphaMcb = 0.1)
										
					# estimation with updated lambdaVec..	
					# refine search region of lambda...		
					minIndex  <- min(length(lambdaVec)-2, max(3, which.min(abs(bestParam[1]-lambdaVec))))
						
					lambdaVec <- c(seq(from = lambdaVec[minIndex-2], to = lambdaVec[minIndex-1], length.out = 6)[2:5],
					   			   seq(from = lambdaVec[minIndex-1], to = lambdaVec[minIndex+0], length.out = 6)[2:5],
								   seq(from = lambdaVec[minIndex+0], to = lambdaVec[minIndex+1], length.out = 6)[2:5],
								   seq(from = lambdaVec[minIndex+1], to = lambdaVec[minIndex+2], length.out = 6)[2:5])			
						
					startValsRefined <- defineSearchRegions(x = DataShifted, lambdaVec = lambdaVec, roundingBase = roundingBase, abEst = startVals$abEst)
								   
					bestParam <- testParam(lambdaVec = lambdaVec, bestParam = bestParam, Data = DataShifted, HistData = Hist, 
										   startValues = startValsRefined, NIter = 8, alpha = 0.01, alphaMcb = 0.1)
					
					# check if costs have improved								
					if(bestParam[5] < bestParamShifted[5])
					{					
						bestParamShifted <- bestParam
					  	bestShift 		 <- shift					
					}	
					
					# update additional shift
					addShift <- ifelse(iterMBC>1, 0, max(0, 0.5*as.numeric(quantile(x = Data, probs = 0.005))))
										
					if(!is.na(bestParam[1]) & !is.na(bestParam[2]) & !is.na(bestParam[3]))
					{						
						addShift   <- max(0, na.omit(c(invBoxCox(qnorm(p=pnorm(-5), mean=bestParam[2], sd=bestParam[3]), lambda=bestParam[1]), addShift))[1])						
					}				
				}
			}										   
							   
			return(c(bestParamShifted, abOr, roundingBase, bestShift))		
				
		}, future.seed=future_lapply(1:length(iterations), FUN = function(x) .Random.seed, future.chunk.size = Inf, future.seed = seed), future.packages = c("ash", "refineR"))
			
	if(length(iterations) > 10)
		plan(sequential)
	
	for (i in iterations) {	
		if (i == 1) {
			obj$Lambda       <- res[[i]][1]
			obj$Mu 	         <- res[[i]][2]
			obj$Sigma        <- res[[i]][3]
			obj$P 	         <- res[[i]][4]
			obj$Shift		 <- res[[i]][10]
			obj$Cost	     <- res[[i]][5]
			obj$abOr         <- res[[i]][7:8]
			obj$roundingBase <- res[[i]][9]	
			obj$rf			 <- res[[i]][6]
			obj$Method	     <- "refineR"
			obj$PkgVersion 	 <- packageVersion("refineR")
			obj$Model		 <- model
		} else {
			obj$LambdaBS[i-1] <- res[[i]][1]
			obj$MuBS[i-1] 	  <- res[[i]][2]
			obj$SigmaBS[i-1]  <- res[[i]][3]	
			obj$PBS[i-1]	  <- res[[i]][4]
			obj$ShiftBS[i-1]  <- res[[i]][10]
			obj$CostBS[i-1]   <- res[[i]][5]			
		}	
	}	
	
	# remove bootstrap samples with result NA
	selNotNA     <- !is.na(obj$LambdaBS)
	obj$LambdaBS <- obj$LambdaBS[selNotNA]
	obj$MuBS     <- obj$MuBS[selNotNA]
	obj$SigmaBS  <- obj$SigmaBS[selNotNA]
	obj$PBS		 <- obj$PBS[selNotNA]
	obj$ShiftBS  <- obj$ShiftBS[selNotNA]
	obj$CostBS   <- obj$CostBS[selNotNA]	
	
	class(obj) <- "RWDRI"
	
	return(obj)
}


#' Helper function to define search regions for mu and sigma and to get the region around the main peak 'ab'
#' 
#' The function estimates the start search regions for mu and sigma for each lambda. Further it determines an appropriate 
#' region around the main peak 'ab' that is used for all lambdas. 
#'
#' @param x				(numeric) values specifying data points comprising pathological
#' 							and non-pathological values 
#' @param lambdaVec		(numeric) transformation parameter for inverse Box-Cox transformation
#' @param roundingBase 	(numeric) describing the rounding base of the dataset
#' @param abEst			(numeric) vector with already estimated abSearchReg and abHist for second definition of search regions 
#'  
#' @return (list) with (abEst, search region for mu and sigma)
#' 
#' 
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

defineSearchRegions <- function(x, lambdaVec, roundingBase, abEst = NULL) {
		
	if (length(x) > 100000) {
		# add jitter to the data in case of rounded data
		if (!is.na(roundingBase))
			x <- x + runif(n=length(x), min =-0.5*roundingBase, max=0.5*roundingBase-1e-20)
					
		# calculate 50000 quantiles type=7 for data (to get better performance)
		qData <- as.numeric(quantile(x = x, probs = seq(0, 1, length.out = 50000), type = 7))
		
	} else {				
		# duplicate dataset and add jitter to the data in case of rounded data
		if (!is.na(roundingBase)) {
			
			x <- rep(x, times=floor(100000/length(x)))			
			x <- x + runif(n=length(x), min =-0.5*roundingBase, max=0.5*roundingBase-1e-20)
			qData <- x
			
		}else {
			qData <- as.numeric(quantile(x = x, probs = seq(0, 1, length.out = 50000), type = 7))
		}			
	}	
	
	# only call estimate ab in first definition of search regions
	if(is.null(abEst))
		abEst <- estimateAB(x=qData) 	# restrict data for estimation of start values to region around main peak
	
	qData <- qData[qData >= abEst$abSearchReg[1] & qData <= abEst$abSearchReg[2]]
	
	muEst 	 <- matrix(nrow = length(lambdaVec), ncol = 3)
	sigmaEst <- matrix(nrow = length(lambdaVec), ncol = 3)
	
	#determine m for density estimation
	densM <- min(48, max(24, ceiling(50 - (1/10)*length(x)/1000)))
	
	c <- 1

	# loop over lambda
	for (lambda in lambdaVec) {
		
		# BoxCox transformation and density estimation...
		# apply BoxCox transformation to quantiles/data and remove NAs an infinity values
		qDataTr <- suppressWarnings(BoxCox(x = qData, lambda = lambda))			
		qDataTr <- qDataTr[!is.na(qDataTr) & is.finite(qDataTr)]		
		
		#get range of transformed data for density estimation
		abTr	 <- range(qDataTr)
		abTrDiff <- diff(abTr)
		abTr[1]  <- abTr[1]-0.1*abTrDiff
		abTr[2]  <- abTr[2]+0.1*abTrDiff
			
		# find main peak of distribution
		mainPeak <- findMainPeak(x=qDataTr, ab=abTr, mStart=densM)			
		peakInd  <- mainPeak$peakInd
		modEst   <- mainPeak$modEst
		qDensTr  <- mainPeak$Dens	
			
		# find width of peak for densities 0.5 to 0.9 of max(density) in 0.01 steps
		steps 	  <- seq(0.5, 0.95, by = 0.05)
		DensRange <- steps*qDensTr$y[peakInd]
		
		widthsL <- NULL
		widthsR <- NULL

		finalSteps <- NULL
		for (i in 1:length(DensRange)) {
			d <- DensRange[i]
			if(length(which(qDensTr$x < modEst & qDensTr$y <= d)) != 0 & (length(which(qDensTr$x > modEst & qDensTr$y <= d)) != 0)){
				# get width (distance) from peak to the left to the first points that is smaller or equal to d
				widthsL <- c(widthsL, modEst - max(qDensTr$x[which(qDensTr$x < modEst & qDensTr$y <= d)]))
				# get width (distance) from peak to the right to the first point that is smaller or equal to d
				widthsR <- c(widthsR, (min(qDensTr$x[which(qDensTr$x > modEst & qDensTr$y <= d)]) - modEst))
				finalSteps <- c(finalSteps, steps[i])
			}
		}
		# get points left and right of the peak 
		left 	<- modEst - widthsL
		right 	<- modEst + widthsR
		
		# leftRight Matrix: cols: 1=density range, 2=left value, 3=right value, 4=widthsL, 5=widthsR, 6=mean of left and right --> mu
		leftRight 	<- matrix(c(finalSteps, left, right, widthsL, widthsR), nrow = length(left), ncol = 5, byrow = FALSE)
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
	
		c = c + 1
	}
	
	return(list("abEst" = abEst, "muEst" = muEst, "sigmaEst" = sigmaEst))	
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
	
	countsData	<- HistData$counts
	
	# Box Cox transformation of histogram breaks and histogram range
	breakL <- suppressWarnings(BoxCox(HistData$breakL, lambda))		
	breakR <- suppressWarnings(BoxCox(HistData$breakR, lambda))		
	
	# parameters for pnormApprox function
	oneOverH <- 1/pNormLookup$h
	pNormVal <- pNormLookup$val
	
	# limits for surrounding peak areas needed fo calculation of P 
	pLimitMin <- c(0.5, 0.55, 0.60, 0.65, 0.70)
	pLimitMax <- c(0.9, 0.94, 0.95, 0.96, 0.97) 
	
	# precalculated factor for confidence region
	qNormFactor <- -qnorm(p = alpha/2)
	
	# maximum number of bins expected to be below confidence region
	maxCountsBelow <- as.integer(ceiling(length(countsData)*alphaMcb/2))
	# initial bestParam
	bestParam <- c(NA, NA, NA, NA, 1e9, NA)
	
	for (mu in muVec) {
		for (sigma in sigmaVec)	{	
			# theoretical prediction of bin counts
			countsPred <- HistData$NData*(pnormApprox(q = breakR, mean = mu, oneOverSd = 1/sigma, oneOverH = oneOverH, pNormVal = pNormVal) - pnormApprox(q = breakL, mean = mu, oneOverSd = 1/sigma, oneOverH = oneOverH, pNormVal = pNormVal))
			countsPred[countsPred < 0] <- 0
			
			# precalculate sqrt of predicted counts
			sqrtCountsPred <- sqrt(countsPred)
						
			# selection of bins that can be approximated by normal distribution
			selectionCounts <- (countsPred >= 40)
			# selection of peak bins (without first and last bin)
			maxCountsPred   <- max(countsPred[2:length(countsPred)-1])
			selectionPeak95 <- (countsPred >= 0.95*maxCountsPred)
			selectionPeak95[length(selectionPeak95)] <- FALSE
			selectionPeak95[1] <- FALSE
			
			selectionPeak80   <- (countsPred >= 0.80*max(countsPred[2:length(countsPred)-1])) 
			selectionPeak80[length(selectionPeak80)] <- FALSE
			selectionPeak80[1] <- FALSE
			
			# lower 0.5% should not be negative, thus check if lowerBound is NA
			RRLowerBound <- invBoxCox(mu-qNormFactor*sigma, lambda = lambda)		
			
			if (!is.na(sum(selectionCounts)) & sum(selectionCounts) > (HistData$NBins*HistData$overlapFactor)/16 &  !is.na(RRLowerBound)) {					
				
				# get sum of data points in specified region next to the peak for the actual data and the prediction
				peakArea  <- getSumForPArea(pLimitMin = pLimitMin, pLimitMax = pLimitMax, countsPred = countsPred, HistData = HistData, lambda = lambda, mu = mu, sigma = sigma)								
				sumDataPeak <- peakArea$sumDataPeak
				sumPredPeak <- peakArea$sumPredPeak	
				
				ratio <- sumDataPeak/sumPredPeak				
				#increase ratio until all exceed data points 
				while (all(ratio <= 1.0) & any(ratio*sumPredPeak-qNormFactor*sqrt(ratio*sumPredPeak) <= sumDataPeak)) {
					rInd <- ratio*sumPredPeak-qNormFactor*sqrt(ratio*sumPredPeak) <= sumDataPeak
					ratio[rInd] <- ratio[rInd] + 0.001
				}	
				PMax <- min(max(0.401, min(ratio)), 1.000)
				
				ratio <- sumDataPeak/sumPredPeak
				#decrease ratio until all fall below the observed data point
				while(all(ratio >= 0.4) & any(ratio*sumPredPeak+qNormFactor*sqrt(ratio*sumPredPeak) >= sumDataPeak)){ 
					rInd <- ratio*sumPredPeak+qNormFactor*sqrt(ratio*sumPredPeak) >= sumDataPeak
					ratio[rInd] <- ratio[rInd] - 0.001
				}		
				PMin <- min(max(0.400, min(ratio)), 0.999)
				
				PVec  <- c(seq(from=PMin, to=PMax, by=min(max(0.002, (PMax-PMin)/8), 0.005)), PMax)
				
				continuePLoop <- TRUE											
				PIndex 		  <- 1
				
				# test different values for P and calculate costs...
				while (continuePLoop) {					
					
					P <- PVec[PIndex]				
					
					# factors applied for relaxation of cost function
					rfVec <- c(5, 3, 1)/1000	
					relaxFactors <- P*maxCountsPred / (rfVec*P*maxCountsPred + sqrt(P*maxCountsPred))^2
					relaxFactors <- c(relaxFactors[relaxFactors>0.001 & relaxFactors<1], 1)								
					
					for(rf in relaxFactors)
					{			
						# calculate actual predicted counts and corresponding confidence region...
						countsPredP <- countsPred*P*rf					
						confWidth 	<- qNormFactor * sqrtCountsPred * sqrt(P*rf)
						countsPredPLower <- countsPredP - confWidth
						countsPredPUpper <- countsPredP + P*confWidth
						
						# select bins within confidence region
						selectionPredBand <- (rf*countsData <= countsPredPUpper & rf*countsData >= countsPredPLower)										
						selectionEval     <- (selectionCounts & selectionPredBand)			
						
						# count number of bins that are selected and not selected for evaluation
						sumSelectionEval  <- sum(selectionEval)
						
						# count number of bins that are at the peak and selected for evaluation
						sumSelectionPeakEval  <- sum(selectionEval & selectionPeak95)
						
						sumCountsEval    <- sum(countsData[ selectionEval & selectionPeak80])
						sumCountsNotEval <- sum(countsData[!selectionEval & selectionPeak80])
						
						ratioCounts80 <- max(0.01, min(1, sumCountsEval/sumCountsNotEval), na.rm=TRUE)^2
						
						# decide if loop shall be continued...
						continuePLoop <- (sum(rf*countsData < countsPredPLower) <= maxCountsBelow & PIndex < length(PVec))
						# if expected number of bins are selected...
						if (sumSelectionEval > (HistData$NBins*HistData$overlapFactor)/32 & (sumSelectionPeakEval>=(0.2*(sum(selectionPeak95))) | sum(selectionPeak95) ==1)) 
						{						
							# extract data of selection
							cData 	 <- rf*countsData[selectionEval]								
							cPred 	 <- countsPredP[selectionEval]						
							
							costPreSum <- sum(log(sqrt(P*rf/(2*pi))*sqrtCountsPred[selectionEval]))					
							
							# calculate costs (negative log likelihood with regularization term)											
							
							# identical to the approximation with dnorm(), but faster
							cost <- -(costPreSum + sumSelectionEval*log(ratioCounts80) + sum(-0.5*(cData-cPred)^2/cPred))/sqrt(sumSelectionEval) 
							
							### derivation of faster cost function calculation 
							
							# original costs with poisson distribution based on statistical assumption of histogram data, with additional regularization
							# pois 	    <- dpois(x = round(cData*rf), lambda = cPred*rf)*cPred*rf*ratioCounts80
							
							# very good approximation of calculation of likelihood (with regularization term)
							# pois <- dnorm(x=cData*rf, mean=cPred*rf, sd=sqrt(cPred*rf))*cPred*rf*ratioCounts80			
							
							# rearrange equation
							# cost <- -sum(log(dnorm(x=cData, mean=cPred, sd=sqrt(cPred))*cPred*ratioCounts80)/sqrt(sumSelectionEval))
							# cost <- -(  sumSelectionEval*log(ratioCounts80) + sum(log(sqrt((P*rf)/(2*pi))*sqrtCountsPred[selectionEval])) +sum((-0.5*(cData - cPred)^2)/cPred ))/sqrt(sumSelectionEval)
							# cost <- -(  sumSelectionEval*log(ratioCounts80) + costPreSum) +sum((-0.5*(cData - cPred)^2)/cPred ))/sqrt(sumSelectionEval)				
							
							# update best cost
							if (!is.na(cost) & cost < bestParam[5]) {						
								bestParam <- c(lambda, mu, sigma, P, cost, rf)							
							}						
						}
					}
					PIndex <- PIndex + 1
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
			
		}else if (peakCrit >= 1 & i == length(diffDens) & is.null(peaksV)){
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
	if (length(valleysV)!= 0 && !is.null(peaksV) && valleysV[1] < peaksV[1]) 
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
#' @param HistData			(list) with histogram data generated by function \code{\link{generateHistData}}
#' @param lambda			(numeric) transformation parameter for inverse Box-Cox transformation
#' @param mu				(numeric) parameter of the mean of non-pathological distribution 
#' @param sigma				(numeric) parameter of the standard deviation of non-pathological distribution
#'  
#' @return (list) with two numeric vectors specifying the amount of observed and estimated data points surrounding the peak 
#'
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}

getSumForPArea <- function(pLimitMin, pLimitMax, countsPred, HistData, lambda, mu, sigma) {
	
	N <- length(countsPred)
	
	countsFullHist <- HistData$fullHist$counts
	breaks    	   <- HistData$fullHist$breaks
	nLimits  	   <- length(pLimitMin)
	
	sumDataPeak <- rep(1, times = 2*nLimits+1)
	
	borderLvec <- rep(0, times = 2*nLimits+1)
	borderRvec <- rep(0, times = 2*nLimits+1)
	
	# determine peak
	indexPeak <- which.max(countsPred)	
	peak      <- countsPred[indexPeak]
	
	pLimitMin <- pLimitMin*peak
	pLimitMax <- pLimitMax*peak
	
	# selection of peak bins
	selectionPeak <- (countsPred >= 0.95*peak)				
	
	#get leftmost and rightmost value of peak bins
	borderL <- min(HistData$breakL[selectionPeak])
	borderR <- max(HistData$breakR[selectionPeak])	
	borderLvec[1] <- borderL
	borderRvec[1] <- borderR
	
	#get observed counts at the peak
	sumDataPeak[1] <- sum(countsFullHist[breaks >= borderL & breaks < borderR])
	
	# split data into left and right part
	countsPredLeft 	<- countsPred[1:indexPeak]	
	breakLLeft 		<- HistData$breakL[1:indexPeak]
	breakRLeft 		<- HistData$breakR[1:indexPeak]
	
	countsPredRight <- countsPred[indexPeak:N]
	breakLRight 	<- HistData$breakL[indexPeak:N]
	breakRRight 	<- HistData$breakR[indexPeak:N]
	
	for (i in 1:length(pLimitMin)) {
		
		#get index of area left from the peak
		selectionLeft <- (countsPredLeft >= pLimitMin[i] & countsPredLeft < pLimitMax[i])
		
		if (sum(selectionLeft) > 0) {
			#get leftmost and rightmost value of selected region
			borderL <- min(breakLLeft[selectionLeft])
			borderR <- max(breakRLeft[selectionLeft])
			borderLvec[2*i] <- borderL
			borderRvec[2*i] <- borderR
			
			#get observed counts in the region next to the peak
			sumDataPeak[2*i] <- sum(countsFullHist[breaks >= borderL & breaks < borderR])
		}	
		
		#get index of area right from the peak
		selectionRight <- (countsPredRight >= pLimitMin[i] & countsPredRight < pLimitMax[i])
		
		if (sum(selectionRight) > 0) {
			#get leftmost and rightmost value of selected region
			borderL <- min(breakLRight[selectionRight])
			borderR <- max(breakRRight[selectionRight])
			borderLvec[2*i+1] <- borderL
			borderRvec[2*i+1] <- borderR
			
			#get observed counts in the region next to the peak
			sumDataPeak[2*i+1] <- sum(countsFullHist[breaks >= borderL & breaks < borderR])
		}			
	}
	
	#get predicted counts at peak and in the regions next to the peak
	sumPredPeak <- HistData$NData*(pnorm(q=BoxCox(borderRvec, lambda), mean = mu, sd = sigma) - pnorm(q = BoxCox(borderLvec, lambda), mean = mu, sd = sigma))
	sumPredPeak[sumPredPeak == 0.0] <- 1
	
	return(list(sumDataPeak = sumDataPeak, sumPredPeak = sumPredPeak))	
}


#' Helper function to find the main peak of a distribution
#' 
#' The function uses a combination of the area under the curve between valleys and the peak height to detect the main peak. 
#' 
#' @param x				(numeric) vector of data points
#' @param ab			(numeric) vector specifying the lower and higher truncation limit of density estimation
#' @param mStart		(integer) specifying the width of the smoothing kernel(s) used for density estimation
#' @param withHeight	(logical) specifying if only the area under the curve (FALSE) or a combination of AUC and peak height (TRUE) should be used
#' 									to detect the main peak
#' @param prevPeak		(numeric) specifying the modEst of the previously estimated peak
#'  
#' @return (list) with the two numeric values peakInd, modEst, and a density list 
#'
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}
 
findMainPeak <- function(x, ab, mStart, withHeight = FALSE, prevPeak =NULL) {
	
	qDensFullTr <- ashDensity(x = x, ab = ab, nbin = 512, m = mStart)
	
	# find peak with largest area
	listPV   <- findPeaksAndValleys(qDensFullTr)
	peaks 	 <- listPV$peaks
	valleys  <- listPV$valleys
	area 	 <- rep(0, times = length(peaks))
	height95 <- rep(1, times = length(peaks))
	
	#get area under peak (between two valleys) 
	for (counter in 1:length(peaks)) {
		if (counter ==1) {
			vLow <- 1
		} else {
			vLow <- valleys[counter-1] #-1 as the function returns the value after the peak/valley
		}
		p <- peaks[counter]
		
		if (counter > length(valleys)) {
			vHigh <- length(qDensFullTr$y)
		} else {
			vHigh <- valleys[counter]
		}
		
		# normalize the two measures used for peak detection
		area[counter] 	  <- sum(qDensFullTr$y[vLow:vHigh])
			
		if(withHeight)
			height95[counter] <- quantile(qDensFullTr$y[vLow:vHigh], probs=0.95)
	}		
	
	area 	 <- area/max(area)	
	height95 <- height95/max(height95)
	locPrevP <- 0
	
	if(withHeight & length(height95)>1)
	{
		maxIndex <- which.max(height95)
		
		ratioHeight <- height95[maxIndex]/max(height95[-maxIndex])
		wHeight95	<- min(1, max(0, (ratioHeight-2)*5 )) 

		maxIndexA <- which.max(area)
		if(max(area[-maxIndexA]) >= 0.95 & wHeight95 == 0)
			wHeight95 <- min(1, max(0, (max(area[-maxIndexA])-0.95)*20)) #0.970 0.975 0.980 0.985 0.990 0.995 1.000 -->  0.4 0.5 0.6 0.7 0.8 0.9 1.0
		
		height95 <- height95*wHeight95
	}	
	
	if(!is.null(prevPeak) & length(peaks) > 1){
		locPrevP <- abs(qDensFullTr$x[peaks]-prevPeak)
		locPrevP <- locPrevP/max(locPrevP)
	}
	
	#get index of peak in area with biggest surface
	peakInd <- min(length(qDensFullTr$x)-1, max(2, peaks[which.max(area+height95-locPrevP)]))
	modEst 	<- sum(qDensFullTr$x[peakInd+c(-1:1)]*qDensFullTr$y[peakInd+c(-1:1)])/sum(qDensFullTr$y[peakInd+c(-1:1)])	
	
	# select region around peak for stronger smoothing
	abPeak	  <- NULL
	abPeak[1] <- max(ab[1], qDensFullTr$x[which(qDensFullTr$x < modEst & qDensFullTr$y <= 0.35*qDensFullTr$y[peakInd])])
	abPeak[2] <- min(ab[2], qDensFullTr$x[which(qDensFullTr$x > modEst & qDensFullTr$y <= 0.35*qDensFullTr$y[peakInd])])
	
	# apply more significant/stronger smoothing in selected region
	densMM <- ceiling(mStart/1.5)
	qDensTr <- ashDensity(x = x, ab = abPeak, nbin = 256, m = densMM)
	
	# find peak with largest area
	listPV   <- findPeaksAndValleys(qDensTr)
	peaks 	 <- listPV$peaks
	valleys  <- listPV$valleys
	area 	 <- rep(0, times=length(peaks))
	height95 <- rep(1, times=length(peaks))
	
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
		
		if(length(vLow) > 0 & length(vHigh) > 0 ){
				
			area[counter] <- sum(qDensTr$y[vLow:vHigh])
			
			if(withHeight)
				height95[counter] <- quantile(qDensTr$y[vLow:vHigh], probs=0.95)
		}
	}		
	
	# normalize the two measures used for peak detection
	area 	 <- area/max(area)
	height95 <- height95/max(height95)
	location <- 0

	if(withHeight & length(height95)>1)
	{
		maxIndex <- which.max(height95)
		
			
		ratioHeight <- height95[maxIndex]/max(height95[-maxIndex])
		
		wHeight95 <- min(1, max(0, (ratioHeight-2)*5 ))
		
		maxIndexA <- which.max(area)
		if(max(area[-maxIndexA]) >= 0.95 & wHeight95 == 0)
			wHeight95 <- min(1, max(0, (max(area[-maxIndexA])-0.95)*20)) #0.970 0.975 0.980 0.985 0.990 0.995 1.000 -->  0.4 0.5 0.6 0.7 0.8 0.9 1.0
		
		height95 <- height95*wHeight95
	}
	
	if(length(peaks) > 1){
		location <- abs(qDensTr$x[peaks]-modEst)
		location <- location/max(location)
	}
	
	#get index of peak in area with biggest surface
	peakInd <- min(length(qDensTr$x)-1, max(2, peaks[which.max(area+height95-location)]))	
	modEst 	<- sum(qDensTr$x[peakInd+c(-1:1)]*qDensTr$y[peakInd+c(-1:1)])/sum(qDensTr$y[peakInd+c(-1:1)])				
	
	return(list(peakInd=peakInd, modEst=modEst, Dens=qDensTr))
}


#' Helper function to find region around the main peak of a distribution
#' 
#' @param x				(numeric) vector of data points
#'  
#' @return (list) with two numeric vectors with lower and upper bound of region around the main peak used for 1) defining the search regions and
#' 			2) estimating the histogram with overlapping bins   
#'
#' @author Tatjana Ammer \email{tatjana.ammer@@roche.com}
 
estimateAB <- function(x) {

	QQ <- as.numeric(quantile(x, probs=c(0, 0.01, 0.95, 0.99, 0.999), na.rm=TRUE))
	
	# default ab used as a starting point for optimization
	abOpt    <- QQ[c(2, 4)]	
	rangeAB  <- diff(abOpt)	
	abOpt[1] <- max(abOpt[1]-0.1*rangeAB, QQ[1], 1e-20)
	
	abOpt[2] <- abHistMax <- min(abOpt[2]+0.1*rangeAB, QQ[5])
	
	# default ab for histogram
	abHist <- abSearchReg <- abOpt

	# ab for peak detection
	abPeakDet    <- QQ[c(2, 4)]		
	abPeakDet[1] <- max(abPeakDet[1], 1e-20)		
	abPeakDet[1] <- abPeakDet[1] - 0.01*diff(abPeakDet) 
	
	stepSize <- 0.5*diff(abOpt)
	peakM 	 <- 24
		
	iter <- 1
	continueABLoop <- TRUE
	
	while(continueABLoop & iter <= 30 & abOpt[2]>abOpt[1] & abOpt[2]>abPeakDet[1] ) {			

		abPeakTemp 		<- c(max(abPeakDet[1], 1e-20), min(abOpt[2], abPeakDet[2]))
		abPeakTemp[1] 	<- abPeakTemp[1] - 0.01*diff(abPeakTemp)

		# find main peak of distribution
		mainPeak  <- findMainPeak(x=x, ab=c(abPeakTemp[1], min(abOpt[2], abPeakDet[2])), mStart=peakM, withHeight=TRUE)
		Dens 	  <- ashDensity(x=x, ab=abOpt, nbin=512, m=24)			
		indexPeak <- which.min(abs(Dens$x-mainPeak$modEst))
		
		# left side of the peak
		yLeft <- rev(Dens$y[1:indexPeak])	
		
		# pre-process density that it becomes monotonic function
		if(length(yLeft) > 1)
		{			
			for(i in 2:length(yLeft))
			{					
				if(yLeft[i] > yLeft[i-1])
					yLeft[i] <- yLeft[i-1]						
			}	
		}

		# right side of the peak
		yRight <- Dens$y[indexPeak:length(Dens$x)]
		
		# pre-process density that it becomes monotonic function
		if(length(yRight) > 1)
		{			
			for(i in 2:length(yRight))
			{					
				if(yRight[i] > yRight[i-1])
					yRight[i] <- yRight[i-1]
			}								
		}
					
		# transform y into weights and calculate weighted proportion of low density bins
		y <- c(yLeft, yRight)
		y <- y/max(y)							
		y <- 2-y*5		
		y <- pmin(pmax(y, 0), 1)
		
		propLowDens <- sum(y)/length(y)
				
		if((iter==1 | propLowDens>0.85) & propLowDens<0.86)	{
			continueABLoop <- FALSE
			abHist 		   <- abOpt						
			abHist[2] 	   <- min(abHist[2] + 0.5*diff(abHist), abHistMax)				
			
		} else if(propLowDens>=0.86) {
			abOpt[2]    <- abOpt[2] - stepSize
			
		} else if(propLowDens<=0.85) {
			abOpt[2]    <- abOpt[2] + stepSize
			
		}				
		
		stepSize <- stepSize/2
		iter	 <- iter + 1	
	}
	
	abOpt <- abSearchReg <- abHist
	
	# find main peak of distribution
	mainPeak  <- findMainPeak(x=x, ab=abOpt, mStart=peakM, withHeight=TRUE, prevPeak = mainPeak$modEst)
	Dens 	  <- ashDensity(x=x, ab=abOpt, nbin=512, m=24)
	indexPeak <- which.min(abs(Dens$x-mainPeak$modEst))
		
	yMin <- yMax <- Dens$y[indexPeak]	
	
	indexL <- max(which((1:length(Dens$x))<indexPeak & Dens$y<0.5*yMax), 2)
	indexR <- min(which((1:length(Dens$x))>indexPeak & Dens$y<0.5*yMax), length(Dens$y)-1)
	
	area50 <- sum(Dens$y[indexL:indexR])
			
	# left side of the peak
	xLeft <- rev(Dens$x[1:(indexL-1)])
	yLeft <- rev(Dens$y[1:(indexL-1)])	
	
	if(length(yLeft) > 1) {			
		
		continueABLoop <- TRUE
		
		for(i in 1:length(yLeft)) {			
			
			if(yLeft[i]<0.5*yMax & yLeft[i]<yMin)
				yMin <- yLeft[i]
			
			if(yLeft[i] > yMin+0.5*yMax & continueABLoop) {
				
				abOpt[1] <- xLeft[i]	
				abSearchReg[1] <- xLeft[i]
				continueABLoop <- FALSE
			}
			
			if(sum(yLeft[1:i])<1.3*area50 & continueABLoop)
				abOpt[1] <- xLeft[i]
			
			if(sum(yLeft[1:i])<1.1*area50 & continueABLoop)
				abSearchReg[1] <- xLeft[i] 
		}	
	}
	
	yMin <- yMax
	
	# right side of the peak		
	xRight <- Dens$x[(indexR+1):length(Dens$x)]
	yRight <- Dens$y[(indexR+1):length(Dens$x)]
	
	if(length(yRight) > 1) {	
		
		continueABLoop <- TRUE
		
		for(i in 1:length(yRight)) {		
			
			if(yRight[i]<0.5*yMax & yRight[i]<yMin)
				yMin <- yRight[i]
			
			if(yRight[i] > yMin+0.5*yMax & continueABLoop) {
				abOpt[2] <- xRight[i]
				abSearchReg[2] <- xRight[i]
				continueABLoop <- FALSE
			}	
			
			if(sum(yRight[1:i])<1.3*area50 & continueABLoop)
				abOpt[2] <- xRight[i]	
			
			
			if(sum(yRight[1:i])<1.1*area50 & continueABLoop)
				abSearchReg[2] <- xRight[i] 
		}								
	}

	abHist <- abOpt
	
	return(list(abHist = abHist, abSearchReg = abSearchReg))
}
