#' Standard plot method for objects of class 'RWDRI'
#' 
#' @param x				(object) of class 'RWDRI'
#' @param Scale			(character) specifying if the plot is generated on the original scale ("Or") or the transformed scale ("Tr")
#' @param RIperc		(numeric) value specifying the percentiles, which define the reference interval
#' @param Nhist			(integer) number of bins in the histogram (derived automatically if not set)
#' @param showCI		(logical) specifying if the confidence intervals are shown
#' @param showPathol	(logical) specifying if the estimated pathological distribution shall be shown
#' @param showValue		(logical) specifying if the exact value of the estimated reference intervals shall be shown above the plot 
#' @param CIprop		(numeric) value specifying the central region for estimation of confidence intervals
#' @param pointEst		(character) specifying the point estimate determination: (1) using the full dataset ("fullDataEst"),
#' 						(2) calculating the median from the bootstrap samples ("medianBS"), (2) works only if NBootstrap > 0
#' @param xlim			(numeric) vector specifying the limits in x-direction	
#' @param ylim			(numeric) vector specifying the limits in y-direction	
#' @param xlab			(character) specifying the x-axis label	
#' @param ylab			(character) specifying the y-axis label	
#' @param title			(character) specifying plot title
#' @param ...			additional arguments passed forward to other functions
#' 
#' @return				No return value. Instead, a plot is generated.
#' 
#' @author Christopher Rank \email{christopher.rank@@roche.com}, Tatjana Ammer \email{tatjana.ammer@@roche.com}
#' 
#' @method plot RWDRI
#' 

plot.RWDRI <- function(x, Scale = c("original", "transformed"), RIperc = c(0.025, 0.975), Nhist = 60, showCI = TRUE, showPathol = TRUE, showValue = TRUE, 
					   CIprop = 0.95, pointEst = c("fullDataEst", "medianBS"), xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, title = NULL, ...) {	
				   
	stopifnot(class(x) == "RWDRI")
	stopifnot(!is.null(x$Data))	
	Scale    <- match.arg(Scale[1], choices = c("original", "transformed"))
	pointEst <- match.arg(pointEst[1], choices = c("fullDataEst", "medianBS"))
	
	modelFound <- (!is.na(x$Mu) & !is.na(x$Sigma) & !is.na(x$Lambda))
	
	# extract binned data
	Data <- x$Data		
	
	if (is.null(x$abOr) | !modelFound)
		ab <- as.numeric(quantile(x = Data, probs = c(0.005, 0.995)))
	else
		ab <- x$abOr
	
	if (Scale == "transformed" & modelFound) {
		Data <- suppressWarnings(BoxCox(Data, x$Lambda))		
		Data <- Data[!is.na(Data) & is.finite(Data)]			
		
		ab <- range(Data)
	}
	
	# calculate reference ranges
	RI <- getRI(x = x, RIperc = RIperc, CIprop = CIprop, pointEst = pointEst, Scale = Scale)	
		
	if (is.null(xlab))	
		xlab <- "Concentration [Units]"
	
	if (is.null(ylab))	
		ylab <- "Frequency"
	
	if (is.null(xlim))
		xlim <- range(c(ab, 0.98*min(RI$PointEst), 0.95*RI$CILow, 1.1*max(RI$PointEst), 1.02*RI$CIHigh), na.rm = TRUE)
	
	if (is.null(title))
		title <- paste0("Estimated Reference Interval (Costs: ", signif(x$Cost, 4), ")")	
	
	if(is.na(x$roundingBase) | Scale == "transformed")
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
			
		breakL 	   <- c(breaks1[1:(length(breaks1)-1)], breaks2[1:(length(breaks2)-1)])
		breakR 	   <- c(breaks1[2:length(breaks1)], 	breaks2[2:length(breaks2)])
	
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
		
		breakL <- breaks1[1:(length(breaks1)-1)]
		breakR <- breaks1[2:length(breaks1)]		
	}	
	
	# Box Cox transformation of histogram breaks and histogram range
	if (Scale == "original" & modelFound) {		
		breakL 	  <- suppressWarnings(BoxCox(breakL, x$Lambda))		
		breakR 	  <- suppressWarnings(BoxCox(breakR, x$Lambda))	
	}		
	
	maxPred <- NA 
	
	if (modelFound) {
		# theoretical prediction of bin counts
		countsPred <- length(Data)*x$P*(pnorm(q = breakR, mean = x$Mu, sd = x$Sigma) - pnorm(q = breakL, mean = x$Mu, sd = x$Sigma))			
		countsPred[countsPred < 0] <- 0
	
		countsPred 	<- countsPred[sortIndex]
		maxPred 	<- max(countsPred)
		
		# calculate difference of counts
		countsDiff <- countsData - countsPred
		countsDiff[countsDiff<0] <- 0
	}
	
	if (is.null(ylim)) {
		ylim <- c(0, 1.03*max(countsData, maxPred, na.rm = TRUE))
		ylim[1] <- 0.03*ylim[2]
	}	
	
	plot(hist1, freq = TRUE, border = NA, col = "grey65", main = title, xaxt = 'n', yaxt = 'n', xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, cex.main = 1.5, cex.lab = 1.25)

	axis(1, at = pretty(xlim, n = 10), cex.axis = 0.85)
	axis(2, at = pretty(ylim), las = 1, cex.axis = ifelse(2/3*max(ylim) < 1e5, 0.9, 0.9), mgp = c(3, ifelse(2/3*max(ylim) < 1e4, 1.0, 0.8), 0))
	
	#axis(2, at = pretty(ylim), las = 1, cex.axis = ifelse(2/3*max(ylim) < 1e5, 0.7, 0.65), mgp = c(3, ifelse(2/3*max(ylim) < 1e4, 1.0, 0.8), 0))
	
	lines(x = mids, y = countsData, lty = 2, lwd = 2, col = "dodgerblue4")
	
	if (modelFound)	{
		lines(x = mids, y = countsPred, lwd = 2, col = "green2")
		
		if(showPathol)
				lines(mids, countsDiff, col = "red3", lwd = 2)
	}
	
	addGrid(pretty(xlim, n = 10), pretty(ylim))
	box()
	
	if (modelFound)	{
		for (i in 1:length(RIperc)) {
			if (showCI & !is.na(RI$CILow[i]) & !is.na(RI$CIHigh[i]))		
				rect(RI$CILow[i],  -1e3, RI$CIHigh[i],  1e9, col = as.rgb("green2", 0.20), border = NA)		
		}	
		
		abline(v = RI$PointEst, lwd = 2, lty = 2, col = "green2")	
		
		adjust <- rep(0.5, times = length(RIperc))
		adjust[RI$Percentile < 0.5] <- 1
		adjust[RI$Percentile > 0.5] <- 0	
		if (showValue)
			mtext(text = signif(RI$PointEst, 3), at = RI$PointEst, col = "green3", cex = 1.3, adj = adjust)
	}
	
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
#' @return No return value, called for adding a grid to a plot
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
#' to specify the transparency within [0,1], 0 meaning completey transparent and 1 meaning completey
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