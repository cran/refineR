#' Standard plot method for objects of class 'RWDRI'
#' 
#' @param x					(object) of class 'RWDRI'
#' @param Scale				(character) specifying if percentiles are shown on the original scale ("or") or the transformed scale ("tr") or the z-Score scale ("z")
#' @param RIperc			(numeric) value specifying the percentiles, which define the reference interval (default c(0.025, 0.975))
#' @param Nhist				(integer) number of bins in the histogram (derived automatically if not set)
#' @param showMargin			(logical) specifying if the specified margins, i.e. confidence intervals or uncertainty margins, shall be shown
#' @param showPathol		(logical) specifying if the estimated pathological distribution shall be shown
#' @param scalePathol		(logical) specifying if the estimated pathological distribution shall be weighted with the ration of pathol/non-pathol
#' @param showBSModels		(logical) specifying if the estimated bootstrapping models shall be shown
#' @param showValue			(logical) specifying if the exact value of the estimated reference intervals shall be shown above the plot
#' @param uncertaintyRegion	(character) specifying the type of the uncertainty region around point estimates 
#' @param CIprop			(numeric) value specifying the central region for estimation of confidence intervals
#' @param UMprop			(numeric) value specifying the central region for estimation of uncertainty margins
#' @param pointEst			(character) specifying the point estimate determination: (1) using the full dataset ("fullDataEst"),
#' 							(2) calculating the median model from the bootstrap samples ("medianBS"), (2) works only if NBootstrap > 0
#' @param colScheme			(character) specifying color scheme of the non-pathological distribution and reference interval; choices are "green" and "blue"
#' @param xlim				(numeric) vector specifying the limits in x-direction	
#' @param ylim				(numeric) vector specifying the limits in y-direction	
#' @param xlab				(character) specifying the x-axis label	
#' @param ylab				(character) specifying the y-axis label	
#' @param title				(character) specifying plot title
#' @param ...				additional arguments passed forward to other functions
#' 
#' @return					The applied plot limits in x-direction (xlim) are returned.
#' 
#' @author Christopher Rank \email{christopher.rank@@roche.com}, Tatjana Ammer \email{tatjana.ammer@@roche.com}
#' 
#' @method plot RWDRI

plot.RWDRI <- function(x, Scale = c("original", "transformed", "zScore"), RIperc = c(0.025, 0.975), Nhist = 60, showMargin = TRUE, showPathol = FALSE, scalePathol = TRUE, showBSModels = FALSE, showValue = TRUE,
					   uncertaintyRegion = c("bootstrapCI", "uncertaintyMargin"), CIprop = 0.95, UMprop = 0.9, pointEst = c("fullDataEst", "medianBS"), colScheme= c("green", "blue"),
					   xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, title = NULL, ...) {	
	
	stopifnot(class(x) == "RWDRI")
	stopifnot(!is.null(x$Data))		
	Scale    <- match.arg(Scale[1], choices = c("original", "transformed", "zScore"))
	uncertaintyRegion    <- match.arg(uncertaintyRegion[1], choices = c("bootstrapCI", "uncertaintyMargin"))
	stopifnot(is.numeric(RIperc) & min(RIperc)>=0 & max(RIperc)<=1)
	stopifnot(is.numeric(Nhist) & Nhist%%1==0 & Nhist>0)
	pointEst <- match.arg(pointEst[1], choices = c("fullDataEst", "medianBS"))
	colScheme <- match.arg(colScheme[1], choices = c("green", "blue"))
	stopifnot(is.null(xlim) | (is.numeric(xlim) & length(xlim)==2))
	stopifnot(is.null(ylim) | (is.numeric(ylim) & length(ylim)==2))
	stopifnot(is.logical(showMargin))
	stopifnot(is.logical(showPathol))
	stopifnot(is.logical(scalePathol))
	stopifnot(is.logical(showBSModels))
	stopifnot(is.logical(showValue))

	if(uncertaintyRegion == "bootstrapCI"){
		stopifnot(is.numeric(CIprop), length(CIprop)==1, CIprop>=0, CIprop<=1)
	}

	args <- list(...)
	asymmetryCorr <- TRUE
	if("asymmetryCorr" %in% names(args) && !(args$asymmetryCorr)){
		asymmetryCorr <- args$asymmetryCorr
	}

	n <- 120
	if("n" %in% names(args) && is.numeric(args$n) && length(args$n) == 1 && args$n > 0)
	{
		n <- args$n
	}

	modelFound   <- (!is.na(x$Mu) & !is.na(x$Sigma) & !is.na(x$Lambda) & !is.na(x$Shift))
	BSPerformed  <- (modelFound & length(x$MuBS) > 0 & length(x$SigmaBS) > 0 & length(x$LambdaBS) > 0 & length(x$ShiftBS) > 0)
	showBSModels <- ifelse(BSPerformed & Scale=="original", showBSModels, FALSE)

	# extract binned data
	Data <- x$Data

	# extract model parameters
	if(modelFound)
	{
		lambda <- ifelse(pointEst=="medianBS" & BSPerformed, x$LambdaMed, x$Lambda)
		mu 	   <- ifelse(pointEst=="medianBS" & BSPerformed, x$MuMed, 	  x$Mu)
		sigma  <- ifelse(pointEst=="medianBS" & BSPerformed, x$SigmaMed,  x$Sigma)
		shift  <- ifelse(pointEst=="medianBS" & BSPerformed, x$ShiftMed,  x$Shift)	
		P	   <- ifelse(pointEst=="medianBS" & BSPerformed, x$PMed,  	  x$P)
	}
	
	# calculate reference intervals
	if(uncertaintyRegion == "uncertaintyMargin"){
		stopifnot(is.numeric(UMprop), length(UMprop)==1, UMprop>=0, UMprop<=1)
		RI <- getRI(x = x, RIperc = RIperc, CIprop = CIprop, UMprop = UMprop, pointEst = pointEst, Scale = Scale, UCMargins=TRUE, asymmetryCorr = asymmetryCorr, n = n)
	} else if(uncertaintyRegion == "bootstrapCI"){
		RI <- getRI(x = x, RIperc = RIperc, CIprop = CIprop, pointEst = pointEst, Scale = Scale, UCMargins=FALSE)
	} 
	
	if (is.null(xlab) & Scale=="original")
		xlab <- "Concentration [Units]"
	else if (is.null(xlab) & Scale=="transformed")	
		xlab <- "Transformed Scale"
	else if (is.null(xlab) & Scale=="zScore")	
		xlab <- "z-Score"		
	
	if (is.null(ylab))	
		ylab <- "Frequency"
		
	# transform data to the correct scale
	if(modelFound & (Scale == "transformed" | Scale == "zScore"))
	{
		Data <- suppressWarnings(BoxCox(Data - shift, lambda))				
		Data <- Data[!is.na(Data) & is.finite(Data)]			
		
		if(Scale == "zScore")
			Data <- (Data - mu) / sigma		
	}	
	
	# determine reasonable xlim
	if (is.null(xlim))
	{
		if (!modelFound)
		{
			rangeData <- as.numeric(quantile(x = Data, probs = c(0.005, 0.995), na.rm = TRUE))
			rangeData <- rangeData + c(-0.02, 0.02)*diff(rangeData)
			
		} else if(modelFound & (Scale == "transformed" | Scale == "zScore"))
		{						
			rangeData <- as.numeric(quantile(x = Data, probs = c(0.001, 0.999), na.rm = TRUE))
			
		} else
		{				
			# calculate skewness of distribution
			skewnessRatio <- diff(getRI(x, RIperc=c(0.025, 0.5, 0.975), UCMargins=FALSE, pointEst = pointEst)$PointEst)
			skewnessRatio <- min(1, sqrt(skewnessRatio[1]/skewnessRatio[2]))
			
			# estimate appropriate concentration range for distribution
			perc595 <- getRI(x, RIperc=c(0.05, 0.95-0.04*(1-skewnessRatio)), UCMargins=FALSE, pointEst = pointEst)$PointEst
			rangeData <- perc595 + c(-1, 1)*1.05*diff(perc595)
			
			# determine appropriate min and max of dataset
			minData <- max(1e-20, quantile(x=Data, probs=0.005, na.rm=TRUE))
			maxData <- quantile(x=Data, probs=0.995, na.rm=TRUE)
			
			# shift range outside of the NP distribution to the right or left when covered by the dataset
			if(rangeData[1] < minData)
			{
				rangeData[2] <- min(rangeData[2]+minData-rangeData[1], maxData)
				rangeData[1] <- minData			
			}
			
			if(rangeData[2] > maxData)
			{
				rangeData[1] <- max(rangeData[1]+maxData-rangeData[2], minData)
				rangeData[2] <- maxData
			}			
		
			rangeData <- rangeData + c(-0.02, 0.02)*diff(rangeData)		
			rangeData[2] <- min(max(Data), rangeData[2])
		}
		
		rangePE <- range(RI$PointEst)
		rangePE <- rangePE + c(-0.06, 0.06)*diff(rangePE)
		
		if(uncertaintyRegion == "bootstrapCI")
			rangeCI	<- range(RI$CILow, RI$CIHigh)
		else
			rangeCI	<- range(RI$UMLow, RI$UMHigh)
				
		rangeCI <- rangeCI + c(-0.03, 0.03)*diff(rangeCI)
		
		xlim <- range(rangeData, rangePE, rangeCI, na.rm=TRUE)
		
		if(Scale == "original")
			xlim[1] <- max(xlim[1], 1e-20)		
	}		
	
	if (is.null(title))
		title <- paste0("Estimated Reference Interval")		
	
	if(is.na(x$roundingBase) | Scale == "transformed" | Scale == "zScore")
	{
		# generate histogram of data
		increment  <- diff(xlim)/Nhist	
		breaks1    <- seq(from = xlim[1] - Nhist*increment, to = xlim[2] + Nhist*increment, by = increment)		
		breaks2	   <- breaks1 + 0.5*increment
		
		if(Scale == "original")
		{
			breaks1    <- breaks1[breaks1 > 1e-20]
			breaks2    <- breaks2[breaks2 > 1e-20]
		}		
		
		hist1  	   <- hist(Data[Data >= min(breaks1) & Data <= max(breaks1)], breaks = breaks1, plot = FALSE)
		hist2  	   <- hist(Data[Data >= min(breaks2) & Data <= max(breaks2)], breaks = breaks2, plot = FALSE)	
		countsData <- c(hist1$counts, hist2$counts)
		mids 	   <- c(hist1$mids, hist2$mids)	
		
		# sort vectors in increasing order
		sortIndex  <- sort(mids, index.return = TRUE)$ix
		countsData <- countsData[sortIndex]
		mids	   <- mids[sortIndex]	
		
		# combine data from hist1 and hist2 that histograms overlap
		hist1$breaks  <- c(mids - 0.25*increment, mids[length(mids)] + 0.25*increment)
		hist1$counts  <- countsData
		hist1$density <- countsData/sum(countsData)
		hist1$mids    <- mids	
		
		breakL 	   <- breakLBS <- c(breaks1[1:(length(breaks1)-1)], breaks2[1:(length(breaks2)-1)])
		breakR 	   <- breakRBS <- c(breaks1[2:length(breaks1)], 	breaks2[2:length(breaks2)])
		
	} else
	{		
		xlimDiff <- diff(xlim)
		binSize <- x$roundingBase*max(1, round(xlimDiff/x$roundingBase/Nhist))
		
		# adapt xlim	
		xlim[1] <- max(0.5*x$roundingBase, round(xlim[1]/x$roundingBase)*x$roundingBase - 0.5*x$roundingBase)		
		xlim[2] <- xlim[1] + ceiling(xlimDiff/binSize)*binSize
		
		breaks1 <- seq(from=xlim[1], to=xlim[2], by=binSize)
		
		hist1  	   <- hist(Data[Data >= min(breaks1) & Data <= max(breaks1)], breaks = breaks1, plot = FALSE)
		
		sortIndex  <- 1:length(hist1$mids)
		mids	   <- hist1$mids 
		countsData <- hist1$counts
		
		breakL <- breakLBS <- breaks1[1:(length(breaks1)-1)]
		breakR <- breakRBS <- breaks1[2:length(breaks1)]			
	}	
	
	# Box Cox transformation of histogram breaks and histogram range
	if (Scale == "original" & modelFound) {		
		breakL 	  <- suppressWarnings(BoxCox(breakL-shift, lambda=lambda))		
		breakR 	  <- suppressWarnings(BoxCox(breakR-shift, lambda=lambda))	
	}		
	
	# calculate curves of BS models
	if (BSPerformed & showBSModels)	
	{
		countsPredBS <- matrix(NA, nrow=length(breakLBS), ncol=length(x$MuBS))
		
		for(i in 1:length(x$MuBS))
		{			
			breakLBS_i <- suppressWarnings(BoxCox(breakLBS-x$ShiftBS[i], lambda=x$LambdaBS[i]))		
			breakRBS_i <- suppressWarnings(BoxCox(breakRBS-x$ShiftBS[i], lambda=x$LambdaBS[i]))
										
			pCorrBS <- BoxCox(c(max(min(x$Data-x$ShiftBS[i]), 1e-20), min(max(x$Data-x$ShiftBS[i]), 1e20)), lambda=x$LambdaBS[i])				
			pCorrBS <- pnorm(q=pCorrBS, mean=x$MuBS[i], sd=x$SigmaBS[i])
			pCorrBS <- 1/(pCorrBS[2]-pCorrBS[1])
		
			tempCountsPredBS <- pCorrBS*length(Data)*x$PBS[i]*(pnorm(q = breakRBS_i, mean = x$MuBS[i], sd = x$SigmaBS[i]) - pnorm(q = breakLBS_i, mean = x$MuBS[i], sd = x$SigmaBS[i]))			
			tempCountsPredBS[tempCountsPredBS < 0] <- 0
			
			tempCountsPredBS  <- tempCountsPredBS[sortIndex]	
			
			countsPredBS[, i] <- tempCountsPredBS
		}		
	}
	
	maxPred <- NA 
	
	if (modelFound) {			
				
		pCorr <- BoxCox(c(max(min(x$Data-shift), 1e-20), min(max(x$Data-shift), 1e20)), lambda=lambda)				
		pCorr <- pnorm(q=pCorr, mean=mu, sd=sigma)
		pCorr <- 1/(pCorr[2]-pCorr[1])
			
		if(Scale == "zScore")
			countsPred <- pCorr*length(Data)*P*(pnorm(q = breakR, mean = 0, sd = 1) - pnorm(q = breakL, mean = 0, sd = 1))
		else
			countsPred <- pCorr*length(Data)*P*(pnorm(q = breakR, mean = mu, sd = sigma) - pnorm(q = breakL, mean = mu, sd = sigma))			
		
		countsPred[countsPred < 0] <- 0
		
		countsPred 	<- countsPred[sortIndex]		
		maxPred 	<- max(countsPred)
		
		# calculate difference of counts
		countsDiff <- countsData - countsPred
		countsDiff[countsDiff<0] <- 0
		
		if (scalePathol) {
			# calculate weighting of difference
			weightDiff <- countsDiff/(countsPred+1e-20)
			weightDiff <- pmin(weightDiff, 1.0)	
			
			countsDiff <- countsDiff*weightDiff				
		}			
	}
	
	if (is.null(ylim)) {
		ylim <- c(0, 1.03*max(countsData, maxPred, na.rm = TRUE))
		ylim[1] <- 0.03*ylim[2]
	}	
	
	if (colScheme == "green")
		nonPatholCol <- c(as.rgb("green3"), as.rgb("green2", 0.20))
	else if (colScheme == "blue")
		nonPatholCol <- c(as.rgb("blue2"), as.rgb("royalblue3", 0.25))
	
	plot(hist1, freq = TRUE, border = NA, col = "grey80", main = title, xaxt = 'n', yaxt = 'n', xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, cex.main = 1.3, cex.lab = 1.05)
	
	axis(1, at = pretty(xlim, n = 10), cex.axis = 0.85)
	axis(2, at = pretty(ylim), las = 1, cex.axis = ifelse(2/3*max(ylim) < 1e5, 0.75, 0.70), mgp = c(3, ifelse(2/3*max(ylim) < 1e4, 1.0, 0.8), 0))
	
	# add curves of BS models
	if (showBSModels) {
		for (i in 1:length(x$MuBS)) {
			
			lines(x = mids, y = countsPredBS[, i], lwd = 2, col = as.rgb(nonPatholCol[1], min(0.25, 4/length(x$MuBS))))			
		}		
	}
	
	# add curve of total distribution
	lines(x = mids, y = countsData, lty = 2, lwd = 1.5, col = "black")
	
	if (modelFound)	{				
		
		# add curve of non-pathological distribution
		if(!showBSModels)
			lines(x = mids, y = countsPred, lwd = 2, col = nonPatholCol[1])
		
		# add curve of pathological distribution
		if(showPathol & !showBSModels)
			lines(mids, countsDiff, col = "red3", lwd = 2)
		if(showPathol & showBSModels)
			warning("showPathol is set to TRUE, but showBSModels is also set to TRUE. The pathological distribution will not be shown in this case.")
	}
	
	addGrid(pretty(xlim, n = 10), pretty(ylim), col = "grey90")
	box()
	
	# add confidence intervals
	if (modelFound)	{
		for (i in 1:length(RIperc)) {
			if (showMargin & uncertaintyRegion=="bootstrapCI" & !is.na(RI$CILow[i]) & !is.na(RI$CIHigh[i]))		
				rect(RI$CILow[i],  -1e3, RI$CIHigh[i], 1e9, col = nonPatholCol[2], border = NA)	
			
			if (showMargin & uncertaintyRegion=="uncertaintyMargin" & !is.na(RI$UMLow[i]) & !is.na(RI$UMHigh[i]))		
				rect(RI$UMLow[i],  -1e3, RI$UMHigh[i], 1e9, col = nonPatholCol[2], border = NA)			
		}	
		
		abline(v = RI$PointEst, lwd = 2, lty = 2, col = nonPatholCol[1])	
			
		selection <- which(RI$PointEst>par("usr")[1] & RI$PointEst<par("usr")[2])
		
		if(showValue & length(selection)>0)
		{			
			adjust <- rep(0.5, times = length(RIperc))
			adjust[RI$Percentile < 0.5] <- 1
			adjust[RI$Percentile > 0.5] <- 0
			
			mtext(text = signif(RI$PointEst[selection], 3), at = RI$PointEst[selection], col = nonPatholCol[1], cex = 1.0, adj = adjust[selection])
		}			
	}
	
	return(invisible(xlim))
}


#' Add a grid to an existing plot.
#' 
#' It is possible to use automatically determined grid lines (\code{x=NULL, y=NULL}) or specifying the number 
#' of cells \code{x = 3, y = 4} as done by \code{grid}. Additionally, x- and y-locations of grid-lines can be specified,
#' e.g. \code{x = 1:10, y = seq(0,10,2)}.
#' 
#' @param x (integer, numeric) single integer specifies number of cells, numeric vector specifies vertical grid-lines
#' @param y (integer, numeric) single integer specifies number of cells, numeric vector specifies horizontal grid-lines
#' @param col (character) color of grid-lines
#' @param lwd (integer) line width of grid-lines
#' @param lty (integer) line type of grid-lines
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}

addGrid <- function(x = NULL, y = NULL, col = "lightgray", lwd = 1L, lty = 3L) {
	if (all(is.null(c(x,y))) || all(length(c(x,y))<2))               # call grid function
		grid(nx = x, ny = y, col = col, lwd = lwd, lty = lty)
	else {
		if (length(x) == 0)                                          # NULL
			xticks <- axTicks(side=1)
		else if (length(x) == 1) {
			U <- par("usr")
			xticks <- seq.int(U[1L], U[2L], length.out = x + 1)
		} else
			xticks <- x
		
		if (length(y) == 0)                                          # NULL
			yticks <- axTicks(side = 2)
		else if (length(y) == 1) {
			U <- par("usr")
			yticks <- seq.int(U[3L], U[4L], length.out = y + 1)
		}
		else
			yticks <- y
		
		abline(v = xticks, col = col, lwd = lwd, lty = lty)
		abline(h = yticks, col = col, lwd = lwd, lty = lty)
	}                                     
}


#' Convert color-names or RGB-code to possibly semi-transparent RGB-code.
#' 
#' Function takes the name of a color and converts it into the rgb space. Parameter "alpha" allows
#' to specify the transparency within (0,1), 0 meaning completey transparent and 1 meaning completey
#' opaque. If an RGB-code is provided and alpha != 1, the RGB-code of the transparency adapted color 
#' will be returned.
#' 
#' @param col (character) name of the color to be converted/transformed into RGB-space (code). Only
#'               those colors can be used which are part of the set returned by function colors(). Defaults
#'               to "black".
#' @param alpha (numeric) value specifying the transparency to be used, 0 = completely transparent, 
#'               1 = opaque.
#' 
#' @return RGB-code
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples 
#' \dontrun{
#' # convert character string representing a color to RGB-code using alpha-channel of .25 (75\% transparent)
#'      as.rgb("red", alpha = .25)
#' 
#' # same thing now using the RGB-code of red (alpha=1, i.e. as.rgb("red"))
#'      as.rgb("#FF0000FF", alpha = .25)
#' }

as.rgb <- function(col = "black", alpha = 1) {
	if (length(col) > 1 && (length(alpha) == 1 || length(alpha) < length(col))) {        # unclear which alpha to use or only one alpha specified
		
		if(length(alpha) < length(col) && length(alpha) > 1)
			warning("Multiple (but too few) 'alpha' specified! Only use 'alpha[1]' for each color!")
		return(sapply(col, as.rgb, alpha = alpha[1]))
	}
	
	if (length(col) > 1 && length(col) <= length(alpha)) {                                # process each color separately
		res <- character()
		for (i in 1:length(col))
			res <- c(res, as.rgb(col[i], alpha[i]))
		return(res)
	}
	
	if ( col %in% colors() )
		return( rgb(t(col2rgb(col))/255, alpha = alpha) )
	else {
		col <- sub("#", "", col)
		R <- as.numeric(paste("0x", substr(col, 1,2), sep = ""))
		G <- as.numeric(paste("0x", substr(col, 3,4), sep = ""))
		B <- as.numeric(paste("0x", substr(col, 5,6), sep = ""))
		return( rgb(R/255, G/255, B/255, alpha = alpha, maxColorValue = 1) )
	}        
}