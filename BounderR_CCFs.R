# TIP FLUORESCENCE IN CCF - GENERAL SCRIPT

# rm(list = ls())

# This script is part of a suite of scripts for analysis of filopodia dynamics 
# using the Fiji plugin Bounder. The questions addressed here are whether the 
# accummulation of protein of interest in tips of filopodia correlates with their
# behaviour. This effect may occur either immediately (offset = 0) or with a delay 
# (offset > 0) if the protein requires time to activate other downstream effectors
# before exerting its effect on tip movement. For this reason we used a cross-correlation
# function to compute correlation (for each filopodium) at each value of the offset.
# We then looked at groups of filopodia that share a similar relationship between 
# fluorescence and movement (responding vs non-responding filopodia) using hierarchical 
# clustering, and compared the properties of those clusters.

# Data input: requires an .Rdata file from the upstream BounderR scripts.

# For more information contact the author (Vasja Urbancic) at vu203@cam.ac.uk.



# ---------------------------------------------------------------------------
# Required packages:

# install.packages("Hmisc", dependencies=TRUE, repos="http://cran.rstudio.com/")
# install.packages("RColorBrewer", dependencies=TRUE, repos="http://cran.rstudio.com/")
# install.packages("wavethresh", dependencies=TRUE, repos="http://cran.rstudio.com/")
library(Hmisc)
library(RColorBrewer)
library(wavethresh)

# Functions (general):

Count <- function(x) length(x[!is.na(x)])			 
SE <- function(x) sd(x, na.rm=TRUE)/sqrt(Count(x))	 							
CI <- function(x) 1.96*sd(x, na.rm=TRUE)/sqrt(Count(x))  

DrawErrorAsPolygon <- function(x, y1, y2, tt, col = 'grey') {
    polygon(c(x[tt], rev(x[tt])), c(y1[tt], rev(y2[tt])), 
    col = col,
    border = NA)			
    }

MovingAverage <- function(x, w = 5) {
		filter(x, rep(1/w, w), sides = 2)
}

# Functions (for block randomisation):

extractBlockIndex <- function(which.block, block.size, ...) {
	start <- ((which.block-1) * block.size) + 1
	end <- ((which.block) * block.size)
	c(start:end)
}


BlockReshuffle <- function(x, block.size = 12) {
	
	stopifnot(length(x) > block.size)
	
	n.blocks <- length(x) %/% block.size
	overhang <- length(x) %% block.size
			
	included <- 1:(block.size*n.blocks)
	excluded.overhang <- setdiff(seq_along(x), included) 
	
	x.in.blocks <- list()
	for(i in 1:n.blocks) {
		x.in.blocks[[i]] <- x[extractBlockIndex(i, 12)]
	}
	
	# which blocks to keep in place (full of NAs), which blocks to swap over?
	
	max.NA.per.block <- 0.25 * block.size 
	blocks.to.shuffle <- which(lapply(x.in.blocks, Count) > max.NA.per.block)
	blocks.to.keep <- which(lapply(x.in.blocks, Count) <= max.NA.per.block)	
	
	# generate permuted blocks, plus insert NA blocks into their respective positions

	set.seed(0.1)
	new.order <- c(sample(blocks.to.shuffle))
	for (j in blocks.to.keep) {
		new.order <- append(new.order, j, after = j-1)
	}
	
	# new vector
		
	for(k in new.order) {
		
		if(exists("z") == FALSE) {z <- c()}
		
		z <- append(z, x.in.blocks[[k]])
	}
	z <- append(z, x[excluded.overhang])
	z	
}




# ---------------------------------------------------------------------------
# 1. Load data from saved workspace

# Load data:

# TOCA (as metalist):
 load("/Users/Lab/Documents/Postdoc/ANALYSIS_local-files/ANALYSIS LOGS/s13_TOCA-in-tips_Preliminary/Analysis1_s13_Renyi4-1_TipFitting/LastWorkspace_TOCA_tracked.Rdata")


# ENA (as simple workspace):
# load("/Users/Lab/Documents/Postdoc/ANALYSIS_local-files/ANALYSIS LOGS/2016-08_s15_ENA-in-tips/LastWorkspace_ENA_tracked.Rdata")


# ---------------------------------------------------------------------------
# 2. Extract equivalent data equally either from simple workspace (ENA) or from within metalist:

if (exists("metalist")) {
	all.dS 	  <- metalist[[1]]$all.dS
	dS.vector <- metalist[[1]]$dS.vector
	bb        <- metalist[[1]]$bb
	max.t     <- metalist[[1]]$max.t
	spt       <- metalist[[1]]$spt
	
	if(exists("metalist[[1]]$all.th.tip.nor1")) {
		tip.f <- metalist[[1]]$all.th.tip.nor1
	} else {
		tip.f <- metalist[[1]]$all.tip.nor1
	}
			
	if(exists("metalist[[1]]$all.move")) {
		all.move <- metalist[[1]]$all.move	
	} else {
		all.move <- metalist[[1]]$all.dctm99  # Use raw measurements not 5-step averaged
	}

} else if (!exists("metalist")) {
	all.dS	  <- all.dS
	dS.vector <- dS.vector
	bb        <- bb
	max.t     <- max.t
	spt       <- spt

	if(exists("all.th.tip.nor1")) {
		tip.f <- all.th.tip.nor1
	} else {
		tip.f <- all.tip.nor1
	}
	
	if(exists("all.move")) {
		all.move <- all.move
	} else {
		all.move <- all.dctm99
	}
}


# Options for using FDCTM instead of raw DCTM, and smoothed tipF signal:

# If use.fdctm == TRUE?

use.fdctm = TRUE

if(use.fdctm == TRUE) {
	if(exists("metalist")) {
		all.move <- metalist[[1]]$all.fdctm
	} else if(exists("all.fdctm")) {
		all.move <- all.fdctm
	}	
}

use.ftip = FALSE

if(use.ftip == TRUE) {
	tip.f <- apply(tip.f, 2, MovingAverage)
}



# Use difference from last timepoint, instead of actual data? (Uncomment if yes.)

# all.move <- apply(all.move, 2, diff)
# tip.f <- apply(tip.f, 2, diff)

# Difference for tip F, raw for movement:
# all.move <- all.move[2:max.t, ]
# tip.f <- apply(tip.f, 2, diff)


# Difference for movement, raw for tip F:
# all.move <- apply(all.move, 2, diff)
# tip.f <- tip.f[2:max.t, ]

# ---------------------------------------------------------------------------
# 3. Necessary data restructuring:

#  3a) - shift up the all.move table by one timepoint: 

start.row <- bb+2
stop.row  <- max.t

if (bb > 0) {
	reshuffle.vec <- c(1:bb, start.row:stop.row, bb+1)
} else if (bb == 0) {
	reshuffle.vec <- c(start.row:stop.row, bb+1)	
}

all.move <- all.move[reshuffle.vec, ]; all.move[max.t, ] <- NA

#  3b) - check if any columns have zero DCTM measurements to remove from dataset 
#      (would trip CCF calculations and heatmaps):

n.timepoints <- colSums( !is.na(all.move)); n.timepoints  
zero.lengths <- which(n.timepoints == 0); zero.lengths

if (length(zero.lengths) > 0) {
	remove.cols <- zero.lengths
	all.move  <- all.move[, -zero.lengths]	
	tip.f         <- tip.f[, -zero.lengths]	
	all.dS	 <- all.dS[, -zero.lengths]
	n.timepoints <- n.timepoints[-zero.lengths]
	rm(remove.cols)
}

short.lengths <- which(n.timepoints < 17); short.lengths

if (length(short.lengths) > 0) {
	remove.cols <- short.lengths
	all.move  <- all.move[, -short.lengths]	
	tip.f         <- tip.f[, -short.lengths]
	all.dS	 <- all.dS[, -short.lengths]
	n.timepoints <- n.timepoints[-short.lengths]
	rm(remove.cols)
}


#  ---------------------------------------------------------------------------
# 4 Derived datasets:

# 4a) - create z scores


# z-scores are calculated with scale(x, scale = TRUE, center = TRUE)
# help(scale)

z.move <- scale(all.move, scale = TRUE, center = TRUE)
z.tip <- scale(tip.f, scale = TRUE, center = TRUE)


# 4b) - create randomised dataset as a control (randomised by shift): 

#all.move.rand <- data.frame(matrix(NA, nrow = nrow(all.move), ncol = ncol(all.move)))
#z.move.rand <- data.frame(matrix(NA, nrow = nrow(z.move), ncol = ncol(z.move)))

#set.seed(0.1)
#set.seed(0.3)
#offsets <- c()

# The way this randomised dataset works is such that auto-correlations in movement and tip f data are 
# preserved - each column is shifted by a random number of timepoints (within the scope of the length
# of the original dataset); shift each movement column by a random number of timepoints (using guyrot, rotational vector shift).
#scopes <- c()

#for (i in 1:n.filo) {
#	orig <- all.move[, i]
#	orig.z <- z.move[, i]
#	scope <- Count(orig)
#	rand.offset <- floor(runif(1, min = 8, max = scope-8))
#	randomised.short <- guyrot(orig[!is.na(orig)], n = rand.offset)
#	randomised.short.z <- guyrot(orig.z[!is.na(orig.z)], n = rand.offset) 
#	new <- orig; new[!is.na(new)] <- randomised.short
#	new.z <- orig.z; new.z[!is.na(new.z)] <- randomised.short.z
	
#	all.move.rand[, i] <- new
#	z.move.rand[, i] <- new.z	
#	offsets[i] <- rand.offset
#	scopes[i] <- scope
#}

# 4c) - create randomised dataset as a control (randomised by block permutation):


all.move.rand <- data.frame(matrix(NA, nrow = nrow(all.move), ncol = ncol(all.move)))
z.move.rand <- data.frame(matrix(NA, nrow = nrow(z.move), ncol = ncol(z.move)))

n.filo <- ncol(all.move)

for (i in 1:n.filo) {
	orig <- all.move[, i]
	orig.z <- z.move[, i]
	new <- BlockReshuffle(orig)
	new.z <- BlockReshuffle(orig.z)
	all.move.rand[, i] <- new
	z.move.rand[, i] <- new.z
	rm(orig, new, orig.z, new.z)
}
colnames(all.move.rand) <- colnames(all.move)
colnames(z.move.rand) <- colnames(z.move)


#  Split all.move into all.ext, all.retr, all.stall

# TipState.Abs <- function(x) {
#   tip.state <- cut(x,
#                    breaks = c(-Inf, threshold.retr.per.t, threshold.ext.per.t, Inf),
#                    labels = c("Retr", "Stall", "Ext"))
#   retrun(summary(tip.state))	
# }

# TipState.Rel <- function(x) {
#   tip.state <- cut(x,
#                    breaks = c(-Inf, threshold.retr.per.t, threshold.ext.per.t, Inf), 
#                    labels = c("Retr", "Stall", "Ext"))	
#   return(summary(tip.state)/Count(tip.state))
# }

all.states <- cut(all.move,
    breaks = c(-Inf, threshold.retr.per.t, threshold.ext.per.t, Inf), 
    labels = c("Retr", "Stall", "Ext"))	

all.ext <- all.move; all.ext[which(all.states != "Ext")] <- NA
all.retr <- all.move; all.retr[which(all.states != "Retr")] <- NA
all.stall <- all.move; all.stall[which(all.states != "Stall")] <- NA

# illustrate how this works:

data.frame("Movement" = all.move[, 2], 
           "Ext" = all.ext[, 2],
           "Stall" = all.stall[, 2],
           "Retr" = all.retr[, 2])[22:121, ]



head(all.states[20:25, 1:10])
  
  
# ---------------------------------------------------------------------------
# 5. Explore correlations (over whole population) with XY scatterplots

# 
dev.new()
par(mfrow = c(2,2))
par(mar = c(4,5,2,1)+0.1)

matplot(tip.f, all.move,
	pch = 16, cex = 0.8,
	col = "#41B6C460", 
	xlab = "Norm tip fluorescence [a.u.]",
#	xlab = expression(Delta * "Tip Fluorescence / Projection Fluorescence [a.u.]"),
	ylab = expression("Tip Movement [" * mu * "m]"),
#	ylab = expression(Delta * "Tip Movement [" * mu * "m]"),
	main = ""
)
	abline(h = 0, lty = 2, col = "grey")
#	abline(v = 1, lty = 2, col = "grey")
	abline(v = 0, lty = 2, col = "grey")

rho <- cor.test(unlist(as.data.frame(tip.f)), unlist(as.data.frame(all.move)), na.action = "na.exclude")$estimate

legend("bottomright", legend = paste("Pearson Rho =", signif(rho, 2)), cex= 0.8, bty = "n")


# As above, with z-scores

# dev.new()
matplot(z.tip, z.move,
	pch = 16, cex = 0.8,
	col = "#41B6C460", 
	xlab = "Norm tip fluorescence [z-score]",
#	xlab = expression(Delta * "Tip Fluorescence / Projection Fluorescence [a.u.]"),
	ylab = expression("Tip Movement [z-score]"),
#	ylab = expression(Delta * "Tip Movement [" * mu * "m]"),
	main = ""
)
	abline(h = 0, lty = 2, col = "grey")
#	abline(v = 1, lty = 2, col = "grey")
	abline(v = 0, lty = 2, col = "grey")

rho.z <- cor.test(unlist(as.data.frame(z.tip)), unlist(as.data.frame(z.move)), na.action = "na.exclude")$estimate

legend("bottomright", legend = paste("Pearson Rho =", signif(rho.z, 2)), cex= 0.8, bty = "n")


# Direct measurements with randomised dataset:

# dev.new()

matplot(tip.f, all.move.rand,
	pch = 16, cex = 0.8,
	col = "#41B6C460", 
	xlab = "Norm tip fluorescence [a.u.]",
#	xlab = expression(Delta * "Tip Fluorescence / Projection Fluorescence [a.u.]"),
	ylab = expression("Tip Movement RANDOMISED [" * mu * "m]"),
#	ylab = expression(Delta * "Tip Movement [" * mu * "m]"),
	main = ""
)
	abline(h = 0, lty = 2, col = "grey")
#	abline(v = 1, lty = 2, col = "grey")
	abline(v = 0, lty = 2, col = "grey")

rho.rand <- cor.test(unlist(as.data.frame(tip.f)), unlist(as.data.frame(all.move.rand)), na.action = "na.exclude")$estimate

legend("bottomright", legend = paste("Pearson Rho =", signif(rho.rand, 2)), cex= 0.8, bty = "n")



# z-scores with randomised dataset:

# dev.new()

matplot(z.tip, z.move.rand,
	pch = 16, cex = 0.8,
	col = "#41B6C460", 
	xlab = "Norm tip fluorescence [z-score]",
#	xlab = expression(Delta * "Tip Fluorescence / Projection Fluorescence [a.u.]"),
	ylab = expression("Tip Movement RANDOMISED [z-score]"),
#	ylab = expression(Delta * "Tip Movement [" * mu * "m]"),
	main = ""
)
	abline(h = 0, lty = 2, col = "grey")
#	abline(v = 1, lty = 2, col = "grey")
	abline(v = 0, lty = 2, col = "grey")

rho.z.rand <- cor.test(unlist(as.data.frame(z.tip)), unlist(as.data.frame(z.move.rand)), na.action = "na.exclude")$estimate

legend("bottomright", legend = paste("Pearson Rho =", signif(rho.z.rand, 2)), cex= 0.8, bty = "n")


# ---------------------------------------------------------------------------
# 6. Calculate CCFs from tip F and tip movement tables

maxlag = 20
lag.range <- -maxlag:maxlag
lag.in.s  <- lag.range * spt

ccf.tip.dctm <- data.frame(matrix(NA, ncol = ncol(all.move), nrow = 2*maxlag + 1))
ccf.tip.dctm.z <- data.frame(matrix(NA, ncol = ncol(all.move), nrow = 2*maxlag + 1))
# ccf.tip.dctm.zweight <- data.frame(matrix(NA, ncol = ncol(all.move), nrow = 2*maxlag + 1))
all.filo  <- seq_along(colnames(all.move))

ccf.tip.ext <- data.frame(matrix(NA, ncol = ncol(all.ext), nrow = 2*maxlag + 1))
ccf.tip.stall <- data.frame(matrix(NA, ncol = ncol(all.stall), nrow = 2*maxlag + 1))
ccf.tip.retr <- data.frame(matrix(NA, ncol = ncol(all.retr), nrow = 2*maxlag + 1))



for (i in all.filo) {
	ccf.i  <- ccf(tip.f[, i], all.move[, i], lag.max = 20, na.action = na.pass, plot = FALSE) 
	ccf.z.i <- ccf(z.tip[, i], z.move[, i], lag.max = 20, na.action = na.pass, plot = FALSE) 
##	ccf.zweight.i <- ccf(z.tip[, i], z.move[, i], lag.max = 20, na.action = na.pass, plot = FALSE) 
	ccf.tip.dctm[, i] <- ccf.i
	ccf.tip.dctm.z[, i] <- ccf.z.i
#	ccf.tip.dctm.zweight[, i] <- ccf.zweight.i
	rm(ccf.i, ccf.z.i)
}

# for (i in all.filo) {
#  ccf.ext.i  <- ccf(tip.f[, i], all.ext[, i], lag.max = 20, na.action = na.pass, plot = FALSE) 
#  ccf.stall.i <- ccf(tip.f[, i], all.stall[, i], lag.max = 20, na.action = na.pass, plot = FALSE) 
#  ccf.retr.i <- ccf(tip.f[, i], all.retr[, i], lag.max = 20, na.action = na.pass, plot = FALSE) 
  
#  ccf.tip.dctm[, i] <- ccf.i
  
  
#  rm(ccf.i, ccf.z.i)
#}



# Repeat for randomised dataset:

ccf.randomised <- data.frame(matrix(NA, ncol = ncol(all.move), nrow = 2*maxlag + 1))
ccf.randomised.z <- data.frame(matrix(NA, ncol = ncol(all.move), nrow = 2*maxlag + 1))

for (i in all.filo) {
	ccf.i  <- ccf(tip.f[, i], all.move.rand[, i], lag.max = 20, na.action = na.pass, plot = FALSE) 
	ccf.z.i <- ccf(z.tip[, i], z.move.rand[, i], lag.max = 20, na.action = na.pass, plot = FALSE) 
	ccf.randomised[, i] <- ccf.i
	ccf.randomised.z[, i] <- ccf.z.i
	rm(ccf.i, ccf.z.i)
}

# Repeat with weighting by z-score:





colnames(ccf.tip.dctm) <- colnames(all.move)
colnames(ccf.tip.dctm.z) <- colnames(all.move)
# colnames(ccf.tip.dctm.zweight) <- colnames(all.move)
row.names(ccf.tip.dctm)  <- lag.in.s
row.names(ccf.tip.dctm.z)  <- lag.in.s
# row.names(ccf.tip.dctm.zweight)  <- lag.in.s
colnames(ccf.randomised) <- colnames(all.move)
colnames(ccf.randomised.z) <- colnames(all.move)
row.names(ccf.randomised)  <- lag.in.s
row.names(ccf.randomised.z)  <- lag.in.s

# check all are there:
lapply(list(ccf.tip.dctm, ccf.tip.dctm.z, ccf.randomised, ccf.randomised.z), head)

# ---------------------------------------------------------------------------
# 7. Compute and plot weighted CCFs  (optional pre-clustering)

#  7a) - Compute weighted CCF metrics:

weights.vec      <- n.timepoints
mean.ccf    <- apply(ccf.tip.dctm, 1, mean, na.rm = TRUE)
w.mean.ccf  <- apply(ccf.tip.dctm, 1, weighted.mean, w = weights.vec, na.rm = TRUE)
w.var.ccf <- apply(ccf.tip.dctm, 1, wtd.var, weights = weights.vec); w.var.ccf
w.sd.ccf  <- sqrt(w.var.ccf); w.sd.ccf
counts.ccf  <- apply(ccf.tip.dctm, 1, Count); counts.ccf
w.ci.ccf  <- 1.96 * w.sd.ccf / sqrt(counts.ccf); w.ci.ccf
ci.ccf = apply(ccf.tip.dctm, 1, CI)


filo.ID.weights <- data.frame("Filo ID" = names(ccf.tip.dctm), "Timepoints" = weights.vec); filo.ID.weights


#  7b) - Plot weighted vs unweighted

dev.new()
matplot(lag.in.s, ccf.tip.dctm, type = "l",
	main = "Cross-correlation of tip fluorescence and movement",
	ylab = "CCF (Tip Fluorescence & DCTM (99%, smoothed))",
	xlab = "Lag [s]",
	col = rgb(0,0,0,0.12),
	lty = 1
	)
abline(v = 0, col = "black", lty = 3)
abline(h = 0, col = "black", lwd = 1)
lines (lag.in.s, w.mean.ccf, 								# RED: new mean (weighted)
	col = 'red',
	lwd = 4)
ci1 = w.mean.ccf + w.ci.ccf
ci2	= w.mean.ccf - w.ci.ccf
DrawErrorAsPolygon(lag.in.s, ci1, ci2, col = rgb(1,0,0,0.2))	
lines (lag.in.s, mean.ccf, 									# BLUE: old mean (unweighted)
	col = 'blue',
	lwd = 4)
ci1 = mean.ccf + ci.ccf
ci2	= mean.ccf - ci.ccf
DrawErrorAsPolygon(lag.in.s, ci1, ci2, col = rgb(0,0,1,0.2))	

text(-40, -0.5, "Mean and 95% CI", pos = 4, col = "blue")
text(-40, -0.6, "Weighted Mean and Weighted 95% CI", col = "red", pos = 4)


#  7c) -  Lines coloured according to weighting:

# ??colorRampPalette

weights.vec
weights.vec2 = weights.vec / max(weights.vec)
palette.Wh.Bu <- colorRampPalette(c("white", "midnightblue"))
palette.Wh.Cor <- colorRampPalette(c("white", "#F37370")) # coral colour palette for second dataset
palette.Wh.Bu(20)
palette.Wh.Cor(20)

# Vector according to which to assign colours:
weights.vec
weights.vec2
weight.interval <- as.numeric(cut(weights.vec, breaks = 10))
w.cols <- palette.Wh.Bu(60)[weight.interval] 
w.cols.Coral <- palette.Wh.Cor(60)[weight.interval] 

data.frame(weights.vec, weights.vec2, weight.interval, w.cols )

dev.new()
	matplot(lag.in.s, ccf.tip.dctm, type = "l",
		col = w.cols,
		lty = 1,
		main = "Cross-correlation of tip fluorescence and movement",
		ylab = "CCF (Tip Fluorescence & Movement)",
		xlab = "Lag [s]"
	)
	abline(v = 0, col = "black", lty = 3)
	abline(h = 0, col = "black", lwd = 1)
	lines(lag.in.s, w.mean.ccf, 								# MIDNIGHTBLUE: new mean (weighted)
		col = 'midnightblue',
		lwd = 4)
ci1 = w.mean.ccf + w.ci.ccf
ci2	= w.mean.ccf - w.ci.ccf
palette.Wh.Bu(20)[20]
palette.Wh.Bu(20)[20]

text(-40, -0.6, "Weighted Mean + 95% CI", col = 'midnightblue', pos = 4)

DrawErrorAsPolygon(lag.in.s, ci1, ci2, col = "#19197020")


# ---------------------------------------------------------------------------
# 8. Heatmaps and clustering

# display.brewer.all()
# ??heatmap

# Data handling: 

# all.ccf.tables <- list("ccf.tip.dctm", "ccf.randomised")
all.ccf.tables <- list(ccf.tip.dctm, ccf.randomised, ccf.tip.dctm.z, ccf.randomised.z)
 # Add derivatives and others
ccf.table.names <- list("tip.dctm", "randomised", "tip.dctm.z", "randomised.z")

names(all.ccf.tables) <- ccf.table.names


# This function creates n clusters from input table (based on euclid distance in rows 18:24) *** Generalise this!

GoCluster <- function(x, n.clusters) {
	map.input <- t(x)
	distance <- dist(map.input[, 18:24], method = "euclidean")
	cluster <- hclust(distance, method = "complete")
	cutree(cluster, k = n.clusters)
}

# e.g.:
# lapply(all.ccf.tables, function(x) GoCluster(x, n.clusters = 2))[[1]]
# lapply(all.ccf.tables, function(x) GoCluster(x, n.clusters = 3))

# This works! (At least some way:)
 sapply(all.ccf.tables, function(x) GoCluster(x, n.clusters = 2))

# This function extracts indices for filo of n-th subcluster within the cluster:

nthSubcluster <- function(x, n.clusters, nth) {	
	which(GoCluster(x, n.clusters = n.clusters) == nth)
}

nthSubclusterOthers <- function(x, n.clusters, nth) {
  which(GoCluster(x, n.clusters = n.clusters) != nth)
}

# e.g.:
# nthSubcluster(ccf.tip.dctm, n.clusters = 2, nth = 1)
# lapply(all.ccf.tables, function(x) nthSubcluster(x, 2, 1))

# ---------
# HEATMAPS:

myHeatmap <- function(x) {
	map.input = t(x)
	distance <- dist(map.input[, 18:24], method = "euclidean")
	cluster <- hclust(distance, method = "complete")
	heatmap(map.input, Rowv = as.dendrogram(cluster), Colv = NA, xlab = "Lag", col = brewer.pal(12, "YlGnBu"), 	scale = "none")	
}
dev.new()
myHeatmap(ccf.tip.dctm[, which(colSums(!is.na(ccf.tip.dctm)) != 0)])
dev.new()
myHeatmap(ccf.randomised[, which(colSums(!is.na(ccf.randomised)) != 0)])


table(GoCluster(ccf.tip.dctm, 5))
table(GoCluster(ccf.tip.dctm, 7))
table(GoCluster(ccf.tip.dctm, 8))
table(GoCluster(ccf.tip.dctm, 9))


# The z-score heatmaps below are identical to the original data:
# dev.new()
# myHeatmap(ccf.randomised.z[, which(colSums(!is.na(ccf.randomised.z)) != 0)])
# dev.new()
# myHeatmap(ccf.tip.dctm.z[, which(colSums(!is.na(ccf.tip.dctm.z)) != 0)])

# heatmaps look identical: check equality here
all.equal(ccf.tip.dctm, ccf.tip.dctm.z) # identical() trips
all(ccf.tip.dctm == ccf.tip.dctm.z) 
# matplot(unlist(ccf.tip.dctm), unlist(ccf.tip.dctm.z), pch = 4, type = "p") # removing any doubts...

# extract values for the heatmap scale min and max:
printEdges <- function(x) print(c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
edges <- lapply(all.ccf.tables, printEdges)
names(all.ccf.tables)

# ---------------------------------------------------------------------------
# 9. CCF LinePlots (per subcluster):

n.timepoints <- colSums( !is.na(all.move)); n.timepoints  
getWeights <- function(x) {colSums(!is.na(x))}
getWeights(all.move)


# *** carefully define weights and subclusters here ***
weights.vec <- getWeights(all.move)	
nthSubcluster(ccf.tip.dctm, 2, 2)

#apply(ccf.tip.dctm, 2, Count)
table(GoCluster(ccf.tip.dctm, 4))
#lapply(ccf.tip.dctm, function(x) nthSubcluster(x, 2, 2))

# e.g. s2.1 = indices of subcluster 1 if split into two clusters
s2.1 <- nthSubcluster(ccf.tip.dctm, 2, 1)
s2.2 <- nthSubcluster(ccf.tip.dctm, 2, 2)
s3.1 <- nthSubcluster(ccf.tip.dctm, 3, 1)
s3.2 <- nthSubcluster(ccf.tip.dctm, 3, 2)
s3.3 <- nthSubcluster(ccf.tip.dctm, 3, 3)
s4.1 <- nthSubcluster(ccf.tip.dctm, 4, 1)
s4.2 <- nthSubcluster(ccf.tip.dctm, 4, 2)
s4.3 <- nthSubcluster(ccf.tip.dctm, 4, 3)
s4.4 <- nthSubcluster(ccf.tip.dctm, 4, 4)

# Most correlated subcluster?

GoCluster(ccf.tip.dctm, 9)
nthSubcluster(ccf.tip.dctm, 9, 3) # Top most correlated! 
s8.3 <- nthSubcluster(ccf.tip.dctm, 8, 3) # Top most correlated! 
s8.others <- nthSubclusterOthers(ccf.tip.dctm, 8, 3)
s8.3; s8.others

# Most correlated subcluster in randomised dataset? (n = 5)

table(GoCluster(ccf.randomised, 7))

GoCluster(ccf.randomised, 7)
nthSubcluster(ccf.randomised, 7, 2)

s7.2 <- nthSubcluster(ccf.randomised, 7, 2)
s7.others <- nthSubclusterOthers(ccf.randomised, 7, 2)

# input: x is a ccf table

CcfLinePlot_2clust <- function(x, ...) {
	w.mean.1 <- apply(x[, s2.1], 1, weighted.mean, w = weights.vec[s2.1], na.rm = TRUE)
	w.mean.2 <- apply(x[, s2.2], 1, weighted.mean, w = weights.vec[s2.2], na.rm = TRUE)
	w.sd.1 <- sqrt( apply(x[, s2.1], 1, wtd.var, w = weights.vec[s2.1], na.rm = TRUE)) 
	w.sd.2 <- sqrt( apply(x[, s2.2], 1, wtd.var, w = weights.vec[s2.2], na.rm = TRUE)) 
	count.1 <- apply(x[, s2.1], 1, Count)
	count.2 <- apply(x[, s2.2], 1, Count)
	w.ci.1 <- 1.96 * w.sd.1 / sqrt(count.1)
	w.ci.2 <- 1.96 * w.sd.2 / sqrt(count.2)
	
	#yhi = max(x, na.rm = TRUE)
	#ylo = min(x, na.rm = TRUE)
	
	dev.new()		
		matplot(lag.in.s, w.mean.1,	
			type = "l",
			main = "",
			xlab = "Lag [s]",
			cex.lab = 1.2,
			cex.axis = 1.2,
			ylab = "CCF (Tip Fluorescence & Movement)",
			col = "#225EA8",
			lwd = 4,
			...)
		matplot(lag.in.s, w.mean.2,	
			add = TRUE,
			type = "l",
			col = "#7FCDBB",
			lwd = 4,
			...)
		abline(v = 0, col = "black", lty = 3)
		abline(h = 0, col = "black", lwd = 1)	
		ci1 = w.mean.1 + w.ci.1
		ci2 = w.mean.1 - w.ci.1
		DrawErrorAsPolygon(lag.in.s, ci1, ci2, col = "#225EA820")
		ci1 = w.mean.2 + w.ci.2
		ci2 = w.mean.2 - w.ci.2
		DrawErrorAsPolygon(lag.in.s, ci1, ci2, col = "#7FCDBB20")
}

CcfLinePlot_2clust.manual <- function(x, c1, c2, ...) {
  w.mean.1 <- apply(x[, c1], 1, weighted.mean, w = weights.vec[c1], na.rm = TRUE)
  w.mean.2 <- apply(x[, c2], 1, weighted.mean, w = weights.vec[c2], na.rm = TRUE)
  w.sd.1 <- sqrt( apply(x[, c1], 1, wtd.var, w = weights.vec[c1], na.rm = TRUE)) 
  w.sd.2 <- sqrt( apply(x[, c2], 1, wtd.var, w = weights.vec[c2], na.rm = TRUE)) 
  count.1 <- apply(x[, c1], 1, Count)
  count.2 <- apply(x[, c2], 1, Count)
  w.ci.1 <- 1.96 * w.sd.1 / sqrt(count.1)
  w.ci.2 <- 1.96 * w.sd.2 / sqrt(count.2)
  
  #yhi = max(x, na.rm = TRUE)
  #ylo = min(x, na.rm = TRUE)
  
 # dev.new()		
  matplot(lag.in.s, w.mean.1,	
          type = "l",
          main = "",
          xlab = "Lag [s]",
          cex.lab = 1.2,
          cex.axis = 1.2,
          ylab = "CCF (Tip Fluorescence & Movement)",
          col = "#225EA8",
          lwd = 4,
          ...)
  matplot(lag.in.s, w.mean.2,	
          add = TRUE,
          type = "l",
          col = "#7FCDBB",
          lwd = 4,
          ...)
  abline(v = 0, col = "black", lty = 3)
  abline(h = 0, col = "black", lwd = 1)	
  ci1 = w.mean.1 + w.ci.1
  ci2 = w.mean.1 - w.ci.1
  DrawErrorAsPolygon(lag.in.s, ci1, ci2, col = "#225EA820")
  ci1 = w.mean.2 + w.ci.2
  ci2 = w.mean.2 - w.ci.2
  DrawErrorAsPolygon(lag.in.s, ci1, ci2, col = "#7FCDBB20")
}

CcfLinePlot_3clust <- function(x, ...) {
	w.mean.1 <- apply(x[, s3.1], 1, weighted.mean, w = weights.vec[s3.1], na.rm = TRUE)
	w.mean.2 <- apply(x[, s3.2], 1, weighted.mean, w = weights.vec[s3.2], na.rm = TRUE)
	w.mean.3 <- apply(x[, s3.3], 1, weighted.mean, w = weights.vec[s3.3], na.rm = TRUE)
	w.sd.1 <- sqrt( apply(x[, s3.1], 1, wtd.var, w = weights.vec[s3.1], na.rm = TRUE)) 
	w.sd.2 <- sqrt( apply(x[, s3.2], 1, wtd.var, w = weights.vec[s3.2], na.rm = TRUE)) 
	w.sd.3 <- sqrt( apply(x[, s3.3], 1, wtd.var, w = weights.vec[s3.3], na.rm = TRUE)) 
	count.1 <- apply(x[, s3.1], 1, Count)
	count.2 <- apply(x[, s3.2], 1, Count)
	count.3 <- apply(x[, s3.3], 1, Count)
	w.ci.1 <- 1.96 * w.sd.1 / sqrt(count.1)
	w.ci.2 <- 1.96 * w.sd.2 / sqrt(count.2)
	w.ci.3 <- 1.96 * w.sd.3 / sqrt(count.3)	
	
	#yhi = max(x, na.rm = TRUE)
	#ylo = min(x, na.rm = TRUE)
	
	curr.cols = brewer.pal(12, "YlGnBu")[c(5, 7, 2)] # "#41B6C4" "#225EA8" "#EDF8B1"
	
	dev.new()		
		matplot(lag.in.s, w.mean.1,	
			type = "l",
			main = "",
			xlab = "Lag [s]",
			cex.lab = 1.2,
			cex.axis = 1.2,
			ylab = "CCF (Tip Fluorescence & Movement)",
			col = "#225EA8",
			lwd = 4,
			...)
		matplot(lag.in.s, w.mean.2,	
			add = TRUE,
			type = "l",
			col = "#7FCDBB",
			lwd = 4,
			...)
		matplot(lag.in.s, w.mean.3,	
			add = TRUE,
			type = "l",
			col = "#EDF8B1",
			lwd = 4,
			...)

		abline(v = 0, col = "black", lty = 3)
		abline(h = 0, col = "black", lwd = 1)	
		ci1 = w.mean.1 + w.ci.1
		ci2 = w.mean.1 - w.ci.1
		DrawErrorAsPolygon(lag.in.s, ci1, ci2, col = "#225EA820")
		ci1 = w.mean.2 + w.ci.2
		ci2 = w.mean.2 - w.ci.2
		DrawErrorAsPolygon(lag.in.s, ci1, ci2, col = "#7FCDBB20")
		ci1 = w.mean.3 + w.ci.3
		ci2 = w.mean.3 - w.ci.3
		DrawErrorAsPolygon(lag.in.s, ci1, ci2, col = "#EDF8B120")
}

CcfLinePlot_4clust <- function(x, ...) {
	w.mean.1 <- apply(x[, s4.1], 1, weighted.mean, w = weights.vec[s4.1], na.rm = TRUE)
	w.mean.2 <- apply(x[, s4.2], 1, weighted.mean, w = weights.vec[s4.2], na.rm = TRUE)
	w.mean.3 <- apply(x[, s4.3], 1, weighted.mean, w = weights.vec[s4.3], na.rm = TRUE)
	w.mean.4 <- apply(x[, s4.4], 1, weighted.mean, w = weights.vec[s4.4], na.rm = TRUE)
	w.sd.1 <- sqrt( apply(x[, s4.1], 1, wtd.var, w = weights.vec[s4.1], na.rm = TRUE)) 
	w.sd.2 <- sqrt( apply(x[, s4.2], 1, wtd.var, w = weights.vec[s4.2], na.rm = TRUE)) 
	w.sd.3 <- sqrt( apply(x[, s4.3], 1, wtd.var, w = weights.vec[s4.3], na.rm = TRUE)) 
	w.sd.4 <- sqrt( apply(x[, s4.4], 1, wtd.var, w = weights.vec[s4.4], na.rm = TRUE)) 
	count.1 <- apply(x[, s4.1], 1, Count)
	count.2 <- apply(x[, s4.2], 1, Count)
	count.3 <- apply(x[, s4.3], 1, Count)
	count.4 <- apply(x[, s4.3], 1, Count)
	w.ci.1 <- 1.96 * w.sd.1 / sqrt(count.1)
	w.ci.2 <- 1.96 * w.sd.2 / sqrt(count.2)
	w.ci.3 <- 1.96 * w.sd.3 / sqrt(count.3)
	w.ci.4 <- 1.96 * w.sd.4 / sqrt(count.4)	
	
	#yhi = max(x, na.rm = TRUE)
	#ylo = min(x, na.rm = TRUE)
	
	curr.cols = brewer.pal(12, "YlGnBu")[c(5, 7, 2)] # "#41B6C4" "#225EA8" "#EDF8B1"
	
	dev.new()		
		matplot(lag.in.s, w.mean.1,	
			type = "l",
			main = "",
			xlab = "Lag [s]",
			cex.lab = 1.2,
			cex.axis = 1.2,
			ylab = "CCF (Tip Fluorescence & Movement)",
			col = "#225EA8",
			lwd = 4,
			...)
		matplot(lag.in.s, w.mean.2,	
			add = TRUE,
			type = "l",
			col = "#7FCDBB",
			lwd = 4,
			...)
		matplot(lag.in.s, w.mean.3,	
			add = TRUE,
			type = "l",
			col = "#EDF8B1",
			lwd = 4,
			...)
		matplot(lag.in.s, w.mean.4,	
			add = TRUE,
			type = "l",
			col = "#EDF8B1",
			lwd = 4,
			...)
		abline(v = 0, col = "black", lty = 3)
		abline(h = 0, col = "black", lwd = 1)	
		ci1 = w.mean.1 + w.ci.1
		ci2 = w.mean.1 - w.ci.1
		DrawErrorAsPolygon(lag.in.s, ci1, ci2, col = "#225EA820")
		ci1 = w.mean.2 + w.ci.2
		ci2 = w.mean.2 - w.ci.2
		DrawErrorAsPolygon(lag.in.s, ci1, ci2, col = "#7FCDBB20")
		ci1 = w.mean.3 + w.ci.3
		ci2 = w.mean.3 - w.ci.3
		DrawErrorAsPolygon(lag.in.s, ci1, ci2, col = "#EDF8B120")
		ci1 = w.mean.4 + w.ci.4
		ci2 = w.mean.4 - w.ci.4
		DrawErrorAsPolygon(lag.in.s, ci1, ci2, col = "#EDF8B120")
}




s2.1 <- nthSubcluster(ccf.tip.dctm, 2, 1)
s2.2 <- nthSubcluster(ccf.tip.dctm, 2, 2)
CcfLinePlot_2clust(ccf.tip.dctm,  ylim = c(-0.2, 0.65))

CcfLinePlot_2clust.manual(ccf.tip.dctm, c1 = s8.3, c2 = s8.others, ylim = c(-0.2, 0.65))
CcfLinePlot_2clust.manual(ccf.randomised, c1 = s8.3, c2 = s8.others, ylim = c(-0.2, 0.65))
par(mfrow = c(1,1))

CcfLinePlot_2clust.manual(ccf.randomised, c1 = s7.2, c2 = s7.others, ylim = c(-0.2, 0.65))
CcfLinePlot_2clust.manual(ccf.tip.dctm, c1 = s7.2, c2 = s7.others, ylim = c(-0.2, 0.65))



s2.1 <- nthSubcluster(ccf.randomised, 2, 1)
s2.2 <- nthSubcluster(ccf.randomised, 2, 2)
CcfLinePlot_2clust(ccf.randomised, ylim = c(-0.4, 0.45))

s3.1 <- nthSubcluster(ccf.tip.dctm[, which(colSums(!is.na(ccf.tip.dctm)) != 0)], 3, 1)
s3.2 <- nthSubcluster(ccf.tip.dctm[, which(colSums(!is.na(ccf.tip.dctm)) != 0)], 3, 2)
s3.3 <- nthSubcluster(ccf.tip.dctm[, which(colSums(!is.na(ccf.tip.dctm)) != 0)], 3, 3)
CcfLinePlot_3clust(ccf.tip.dctm[, which(colSums(!is.na(ccf.tip.dctm)) != 0)],  ylim = c(-0.8, 0.8))

s3.1 <- nthSubcluster(ccf.randomised[, which(colSums(!is.na(ccf.randomised)) != 0)], 3, 1)
s3.2 <- nthSubcluster(ccf.randomised[, which(colSums(!is.na(ccf.randomised)) != 0)], 3, 2)
s3.3 <- nthSubcluster(ccf.randomised[, which(colSums(!is.na(ccf.randomised)) != 0)], 3, 3)
CcfLinePlot_3clust(ccf.randomised[, which(colSums(!is.na(ccf.randomised)) != 0)],  ylim = c(-0.8, 0.8))

# Repeat for z-scores - not much use, these are identical!
s3.1 <- nthSubcluster(ccf.tip.dctm.z[, which(colSums(!is.na(ccf.tip.dctm.z)) != 0)], 3, 1)
s3.2 <- nthSubcluster(ccf.tip.dctm.z[, which(colSums(!is.na(ccf.tip.dctm.z)) != 0)], 3, 2)
s3.3 <- nthSubcluster(ccf.tip.dctm.z[, which(colSums(!is.na(ccf.tip.dctm.z)) != 0)], 3, 3)
CcfLinePlot_3clust(ccf.tip.dctm.z[, which(colSums(!is.na(ccf.tip.dctm.z)) != 0)],  ylim = c(-0.8, 0.8))

s3.1 <- nthSubcluster(ccf.randomised.z[, which(colSums(!is.na(ccf.randomised.z)) != 0)], 3, 1)
s3.2 <- nthSubcluster(ccf.randomised.z[, which(colSums(!is.na(ccf.randomised.z)) != 0)], 3, 2)
s3.3 <- nthSubcluster(ccf.randomised.z[, which(colSums(!is.na(ccf.randomised.z)) != 0)], 3, 3)
CcfLinePlot_3clust(ccf.randomised.z[, which(colSums(!is.na(ccf.randomised.z)) != 0)],  ylim = c(-0.8, 0.8))


s4.1 <- nthSubcluster(ccf.tip.dctm[, which(colSums(!is.na(ccf.tip.dctm)) != 0)], 4, 1)
s4.2 <- nthSubcluster(ccf.tip.dctm[, which(colSums(!is.na(ccf.tip.dctm)) != 0)], 4, 2)
s4.3 <- nthSubcluster(ccf.tip.dctm[, which(colSums(!is.na(ccf.tip.dctm)) != 0)], 4, 3)
s4.4 <- nthSubcluster(ccf.tip.dctm[, which(colSums(!is.na(ccf.tip.dctm)) != 0)], 4, 4)

CcfLinePlot_4clust(ccf.tip.dctm[, which(colSums(!is.na(ccf.tip.dctm)) != 0)],  ylim = c(-0.6, 0.6))

s4.1 <- nthSubcluster(ccf.randomised[, which(colSums(!is.na(ccf.randomised)) != 0)], 4, 1)
s4.2 <- nthSubcluster(ccf.randomised[, which(colSums(!is.na(ccf.randomised)) != 0)], 4, 2)
s4.3 <- nthSubcluster(ccf.randomised[, which(colSums(!is.na(ccf.randomised)) != 0)], 4, 3)
s4.4 <- nthSubcluster(ccf.randomised[, which(colSums(!is.na(ccf.randomised)) != 0)], 4, 4)
CcfLinePlot_4clust(ccf.randomised[, which(colSums(!is.na(ccf.randomised)) != 0)],  ylim = c(-0.6, 0.6))


# ---------------------------------------------------------------------------
# 10. XY scatterplots (per subcluster): 


TripleGraphXY <- function(x, y, n.cluster = 3, ...) {
	
	# Clustering now is done using the original dataset! change this to clustering on the input dataset...
	
	ccf.table <- data.frame(matrix(NA, ncol = ncol(x), nrow = 2*maxlag + 1))
	all.filo  <- seq_along(colnames(x))
	for (i in all.filo) {
		ccf.i  <- ccf(x[, i], y[, i], lag.max = 20, na.action = na.pass, plot = 	FALSE) 
		ccf.table[, i] <- ccf.i
		rm(ccf.i)
	}
	colnames(ccf.table) <- colnames(all.move)
	row.names(ccf.table) <- lag.in.s

	curr.cluster.1 <- nthSubcluster(ccf.table, n.cluster, 1)
	curr.cluster.2 <- nthSubcluster(ccf.table, n.cluster, 2)
	curr.cluster.3 <- nthSubcluster(ccf.table, n.cluster, 3)

###	
# x = z.tip
# y = z.move

	rho.1 <- cor.test(unlist(as.data.frame(x[, curr.cluster.1])), unlist(as.data.frame(y[, curr.cluster.1])), na.action = "na.exclude")$estimate
	rho.2 <- cor.test(unlist(as.data.frame(x[, curr.cluster.2])), unlist(as.data.frame(y[, curr.cluster.2])), na.action = "na.exclude")$estimate
	rho.3 <- cor.test(unlist(as.data.frame(x[, curr.cluster.3])), unlist(as.data.frame(y[, curr.cluster.3])), na.action = "na.exclude")$estimate
		
	# dev.new()
	par(mfrow = c(1,3))
		matplot(x[, curr.cluster.1], y[, curr.cluster.1], ...); 
		legend("bottomright", legend = paste("Pearson Rho =", signif(rho.1, 2)), cex= 1, bty = "n")
		
		matplot(x[, curr.cluster.2], y[, curr.cluster.2], ...)	
		legend("bottomright", legend = paste("Pearson Rho =", signif(rho.2, 2)), cex= 1, bty = "n")

		matplot(x[, curr.cluster.3], y[, curr.cluster.3], ...)	
		legend("bottomright", legend = paste("Pearson Rho =", signif(rho.3, 2)), cex= 1, bty = "n")
}


# With paired datasets (as measured and as z-scores):

TripleGraphXY(tip.f, all.move, 3, 
	pch = 4, cex = 0.4,
	xlab = "Norm tip fluorescence [a.u.]", 
	ylab = expression("Tip movement [" * mu * "m]") )
	
TripleGraphXY(z.tip, z.move, 3, 
	pch = 4, cex = 0.4, 
	xlab = "Norm tip fluorescence [z-score]", 
	ylab = "Tip movement [z-score]")
	
	
# With randomised datasets (as measured and as z-scores):	
	
TripleGraphXY(tip.f, all.move.rand, 3, 
	pch = 4, cex = 0.4,
	xlab = "Norm tip fluorescence [a.u.]", 
	ylab = expression("Tip movement [" * mu * "m]") )

TripleGraphXY(z.tip, z.move.rand, 3, 
	pch = 4, cex = 0.4,
	xlab = "Norm tip fluorescence [z-score]", 
	ylab = "Tip movement [z-score]")



DoubleGraphXY <- function(x, y, n.cluster = 2, ...) {
	
	# Clustering now is done using the original dataset! change this to clustering on the input dataset...
	
	ccf.table <- data.frame(matrix(NA, ncol = ncol(x), nrow = 2*maxlag + 1))
	all.filo  <- seq_along(colnames(x))
	for (i in all.filo) {
		ccf.i  <- ccf(x[, i], y[, i], lag.max = 20, na.action = na.pass, plot = 	FALSE) 
		ccf.table[, i] <- ccf.i
		rm(ccf.i)
	}
	colnames(ccf.table) <- colnames(all.move)
	row.names(ccf.table) <- lag.in.s

	curr.cluster.1 <- nthSubcluster(ccf.table, n.cluster, 1)
	curr.cluster.2 <- nthSubcluster(ccf.table, n.cluster, 2)

	rho.1 <- cor.test(unlist(as.data.frame(x[, curr.cluster.1])), unlist(as.data.frame(y[, curr.cluster.1])), na.action = "na.exclude")$estimate
	rho.2 <- cor.test(unlist(as.data.frame(x[, curr.cluster.2])), unlist(as.data.frame(y[, curr.cluster.2])), na.action = "na.exclude")$estimate
			
	# dev.new()
	par(mfrow = c(1,2))
		matplot(x[, curr.cluster.1], y[, curr.cluster.1], ...); 
		legend("bottomright", legend = paste("Pearson Rho =", signif(rho.1, 2)), cex= 1, bty = "n")
		abline(h = 0, lty = 2, col = "grey"); abline(v = 0, lty = 2, col = "grey")
		
		matplot(x[, curr.cluster.2], y[, curr.cluster.2], ...)	
		legend("bottomright", legend = paste("Pearson Rho =", signif(rho.2, 2)), cex= 1, bty = "n")
		abline(h = 0, lty = 2, col = "grey"); abline(v = 0, lty = 2, col = "grey")
}

# DoubleGraph with MANUAL SUBCLUSTER selection

DoubleGraphXY.manual <- function(xtable, ytable, c1, c2, ...) {
  
  x1 <- xtable[, c1]
  y1 <- ytable[, c1]
  x2 <- xtable[, c2]
  y2 <- ytable[, c2]
  
  rho.1 <- cor.test(unlist(as.data.frame(x1)), unlist(as.data.frame(y1)), na.action = "na.exclude")$estimate
  rho.2<- cor.test(unlist(as.data.frame(x2)), unlist(as.data.frame(y2)), na.action = "na.exclude")$estimate
  
  par(mfrow = c(1,2))
  
  matplot(x1, y1, ...);
    legend("bottomright", legend = paste("Pearson Rho =", signif(rho.1, 2)), cex= 1, bty = "n")
    abline(h = 0, lty = 2, col = "grey"); abline(v = 0, lty = 2, col = "grey")
  matplot(x2, y2, ...);
    legend("bottomright", legend = paste("Pearson Rho =", signif(rho.2, 2)), cex= 1, bty = "n")
    abline(h = 0, lty = 2, col = "grey"); abline(v = 0, lty = 2, col = "grey")
}


DoubleGraphXY.manual(z.tip, z.move, c1 = s8.3, c2 = s8.others, 
                    type = "p", pch = 4, cex = 0.4,
                    xlab = "Norm tip fluorescence [z-score]",
                    ylab = expression("Tip movement [z-score]")
                    )

DoubleGraphXY.manual(z.tip, z.move.rand, c1 = s7.2, c2 = s7.others, 
                     type = "p", pch = 4, cex = 0.4,
                     xlab = "Norm tip fluorescence [z-score]",
                     ylab = expression("Tip movement [z-score]")
)

# With paired datasets (as measured and as z-scores):

DoubleGraphXY(tip.f, all.move, 2, 
	pch = 4, cex = 0.4,
	xlab = "Norm tip fluorescence [a.u.]", 
	ylab = expression("Tip movement [" * mu * "m]") )
	
DoubleGraphXY(z.tip, z.move, 2, 
	pch = 4, cex = 0.4, 
	xlab = "Norm tip fluorescence [z-score]", 
	ylab = "Tip movement [z-score]")
	
# With randomised datasets (as measured and as z-scores):	
	
DoubleGraphXY(tip.f, all.move.rand, 2, 
	pch = 4, cex = 0.4,
	xlab = "Norm tip fluorescence [a.u.]", 
	ylab = expression("Tip movement [" * mu * "m]") )

DoubleGraphXY(z.tip, z.move.rand, 2, 
	pch = 4, cex = 0.4,
	xlab = "Norm tip fluorescence [z-score]", 
	ylab = "Tip movement [z-score]")


# dev.new()
# par(mfrow = c(1,3))


# ---------------------------------------------------------------------------
# 11. Properties per cluster? 
#	- No. timepoints per track?
#	- "Eventfulness"?


 s2.1 <- nthSubcluster(ccf.tip.dctm, 2, 1)
 s2.2 <- nthSubcluster(ccf.tip.dctm, 2, 2)

# s3.1 <- nthSubcluster(ccf.tip.dctm, 3, 1)  # Uncomment for 3 clusters
# s3.2 <- nthSubcluster(ccf.tip.dctm, 3, 2)
# s3.3 <- nthSubcluster(ccf.tip.dctm, 3, 3)

# --- TRACK DURATION

# timepoints.1 <- apply(tip.f[, s2.1], 2, Count)
# timepoints.2 <- apply(tip.f[, s2.2], 2, Count)

 timepoints.1 <- apply(all.move[, s3.1], 2, Count)  # Uncomment for 3 clusters
 timepoints.2 <- apply(all.move[, s3.2], 2, Count)
 timepoints.3 <- apply(all.move[, s3.3], 2, Count)



timepoints.by.cluster <- data.frame(
	# Concatenate all timepoints measurements (counting movement measurements)
	"Timepoints" = c(timepoints.1, timepoints.2
	, timepoints.3    ### Uncomment for 3 clusters
	),
	# New column giving category
	"Subcluster" = c(
		rep("Subcluster1", length(timepoints.1)),
		rep("Subcluster2", length(timepoints.2))
		, rep("Subcluster3", length(timepoints.3))  ### Uncomment for 3 clusters
)
)

dev.new(width = 7, height = 3)
par(mfrow = c(1,3))

boxplot(Timepoints ~ Subcluster, timepoints.by.cluster,
	main = "Timepoints per subcluster",
	col = c("#225EA850","#7FCDBB50", "#EDF8B150"),
	outpch = NA)
stripchart(Timepoints ~ Subcluster, timepoints.by.cluster,
	add = TRUE,
	vertical = TRUE,
    method = "jitter",
	pch = 4,
	cex = 1)

 kruskal.test(Timepoints ~ Subcluster, data = timepoints.by.cluster) # p = 0.024 for ENA
 kruskal.p <- kruskal.test(Timepoints ~ Subcluster, data = timepoints.by.cluster)$p.value
 legend("bottomright", legend = paste("Kruskal P =", signif(kruskal.p, 2)), cex= 1, bty = "n")


# WITH MANUAL SUBCLUSTER (8.3)

timepoints.8.3 <- apply(all.move[, s8.3], 2, Count)
timepoints.8.others <- apply(all.move[, s8.others], 2, Count)

timepoints.by.cluster8.3 <- data.frame(
  "Timepoints" = c(timepoints.8.3, timepoints.8.others),
  "Subcluster" = c(
    rep("s8.3", length(timepoints.8.3)),
    rep("Other", length(timepoints.8.others))
  )
)
par(mfrow = c(1,3))
boxplot(Timepoints ~ Subcluster, timepoints.by.cluster8.3,
        main = "Timepoints",
        col = c("#7FCDBB50","#225EA850", "#EDF8B150"),
        outpch = NA)
stripchart(Timepoints ~ Subcluster, timepoints.by.cluster8.3,
           add = TRUE,
           vertical = TRUE,
           method = "jitter",
           pch = 4,
           cex = 1)

kruskal.test(Timepoints ~ Subcluster, data = timepoints.by.cluster8.3) # p = 0.024 for ENA
kruskal.p <- kruskal.test(Timepoints ~ Subcluster, data = timepoints.by.cluster8.3)$p.value
legend("bottomright", legend = paste("Kruskal P =", signif(kruskal.p, 2)), cex= 1, bty = "n")


# wilcox.test(Timepoints ~ Subcluster, data = timepoints.by.cluster) 
# mw.p <- wilcox.test(Timepoints ~ Subcluster, data = timepoints.by.cluster)$p.value
# legend("bottomright", legend = paste("Mann-Wh P =", signif(mw.p, 2)), cex= 1, bty = "n")


nrow(subset(timepoints.by.cluster, Timepoints > 20))  # 41 tracks have more than 20 timepoints
nrow(subset(timepoints.by.cluster, Timepoints < 20))  # 13 tracks have less than 20 timepoints


# --- TIP F VARIANCE


# tipvar.1 <- apply(tip.f[, s2.1], 2, var, na.rm = TRUE)
# tipvar.2 <- apply(tip.f[, s2.2], 2, var, na.rm = TRUE)
tipvar.1 <- apply(tip.f[, s3.1], 2, var, na.rm = TRUE)
tipvar.2 <- apply(tip.f[, s3.2], 2, var, na.rm = TRUE)
tipvar.3 <- apply(tip.f[, s3.3], 2, var, na.rm = TRUE)

# --- or TIP MOVEMENT VARIANCE:

# tipvar.1 <- apply(all.move[, s3.1], 2, var, na.rm = TRUE)
# tipvar.2 <- apply(all.move[, s3.2], 2, var, na.rm = TRUE)
# tipvar.3 <- apply(all.move[, s3.3], 2, var, na.rm = TRUE) 

 



tipvar.by.cluster <- data.frame(
	# Concatenate all timepoints measurements (counting movement measurements)
	"Tipvar" = c(tipvar.1, tipvar.2
	, tipvar.3    ### Uncomment for 3 clusters
	),
	# New column giving category
	"Subcluster" = c(
		rep("Subcluster1", length(tipvar.1)),
		rep("Subcluster2", length(tipvar.2))
	, rep("Subcluster3", length(timepoints.3))  ### Uncomment for 3 clusters
)
)

# dev.new()
boxplot(Tipvar ~ Subcluster, tipvar.by.cluster,
	main = "Tip F variance per subcluster",
	col = c("#225EA850","#7FCDBB50", "#EDF8B150"),
	outpch = NA)
stripchart(Tipvar ~ Subcluster, tipvar.by.cluster,
	add = TRUE,
	vertical = TRUE,
    method = "jitter",
	pch = 4,
	cex = 0.4)

tipvar.by.cluster

# kruskal.test(Timepoints ~ Subcluster, data = timepoints.by.cluster) 
# kruskal.p <- kruskal.test(Timepoints ~ Subcluster, data = timepoints.by.cluster)$p.value
# legend("bottomright", legend = paste("Kruskal P =", signif(kruskal.p, 2)), cex= 1, bty = "n")

kruskal.test(Tipvar ~ Subcluster, data = tipvar.by.cluster) 
kr.p <- kruskal.test(Tipvar ~ Subcluster, data = tipvar.by.cluster)$p.value
legend("topleft", legend = paste("Kruskal P =", signif(kr.p, 2)), cex= 1, bty = "n")


# WITH MANUAL SUBCLUSTER (8.3)

tipvar.1 <- apply(tip.f[, s8.3], 2, var, na.rm = TRUE)
tipvar.2 <- apply(tip.f[, s8.others], 2, var, na.rm = TRUE)

tipvar.by.cluster <- data.frame(
  # Concatenate all timepoints measurements (counting movement measurements)
  "Tipvar" = c(tipvar.1, tipvar.2),
  # New column giving category
  "Subcluster" = c(
    rep("s8.3", length(tipvar.1)),
    rep("Other", length(tipvar.2))
  )
)

boxplot(Tipvar ~ Subcluster, tipvar.by.cluster,
        main = "Tip F variance",
        col = c("#7FCDBB50","#225EA850", "#EDF8B150"),
        outpch = NA)
stripchart(Tipvar ~ Subcluster, tipvar.by.cluster,
           add = TRUE,
           vertical = TRUE,
           method = "jitter",
           pch = 4,
           cex = 0.4)

kruskal.test(Tipvar ~ Subcluster, data = tipvar.by.cluster) 
kr.p <- kruskal.test(Tipvar ~ Subcluster, data = tipvar.by.cluster)$p.value
legend("topleft", legend = paste("Kruskal P =", signif(kr.p, 2)), cex= 1, bty = "n")



# Tip F autocorrelation

MaxLag <- 120

AcfTable <- function(x, L) {
	y <- data.frame(matrix(NA, ncol = ncol(x), nrow = L + 1))
	colnames(y) <- colnames(x)
	for (i in 1:ncol(x)) {
		acf.i <- as.vector(acf(x[, i], na.action = na.pass, lag.max = MaxLag, plot = FALSE)[[1]])
		y[, i] <- acf.i
		rm(acf.i)
		}
	y	
}
acf.tip.f <- AcfTable(tip.f, 120) 

# Extract first negative in column (the more strongly correlated the timepoints, the longer it takes for ACF to decay to 0; the offset of decay to 0 is taken as a measure of consistency [" root of ACF (TipF)")]

FirstNegative <- function(x) {
	if(sum(!is.na(x)) > 0) {				
		min (which(x <= 0), na.rm = TRUE)	# requires at least one non-NA value in vector (if >3 elements in dctm column acf is not computed, so acf is NA throughout, which returns Inf) 
	} else {
		NA									# returns NA in place of Inf
	}
}

acf.tip.f.roots <- apply (acf.tip.f, 2, FirstNegative)

#	hist(acf.tip.f.roots, breaks = 10, col = 'salmon', border = 'white')

acf.tip.f.roots

tipacf.1 <- acf.tip.f.roots[s3.1]
tipacf.2 <- acf.tip.f.roots[s3.2]
tipacf.3 <- acf.tip.f.roots[s3.3]

tipacf.by.cluster <- data.frame(
	# Concatenate all timepoints measurements (counting movement measurements)
	"RootACF" = c(tipacf.1, tipacf.2
	, tipacf.3    ### Uncomment for 3 clusters
	),
	# New column giving category
	"Subcluster" = c(
		rep("Subcluster1", length(tipacf.1)),
		rep("Subcluster2", length(tipacf.2))
		, rep("Subcluster3", length(timepoints.3))  ### Uncomment for 3 clusters
)
)
# tipacf.by.cluster

# dev.new()
boxplot(RootACF ~ Subcluster, tipacf.by.cluster,
	main = "Persistence of TipF per subcluster",
	col = c("#225EA850","#7FCDBB50", "#EDF8B150"),
	outpch = NA)
stripchart(RootACF ~ Subcluster, tipacf.by.cluster,
	add = TRUE,
	vertical = TRUE,
    method = "jitter",
	pch = 4,
	cex = 0.4)

kruskal.test(RootACF ~ Subcluster, data = tipacf.by.cluster) 
kr.p <- kruskal.test(RootACF ~ Subcluster, data = tipacf.by.cluster)$p.value
legend("topleft", legend = paste("Kruskal P =", signif(kr.p, 2)), cex= 1, bty = "n")

# FOR MANUAL s8.3

tipacf.1 <- acf.tip.f.roots[s8.3]
tipacf.2 <- acf.tip.f.roots[s8.others]

tipacf.by.cluster <- data.frame(
  "RootACF" = c(tipacf.1, tipacf.2),
  "Subcluster" = c(
    rep("s8.3", length(tipacf.1)),
    rep("Other", length(tipacf.2))
    )
  )
# tipacf.by.cluster

# dev.new()

boxplot(RootACF ~ Subcluster, tipacf.by.cluster,
        main = "Persistence of TipF",
        col = c("#7FCDBB50","#225EA850", "#EDF8B150"),
        outpch = NA)
stripchart(RootACF ~ Subcluster, tipacf.by.cluster,
           add = TRUE,
           vertical = TRUE,
           method = "jitter",
           pch = 4,
           cex = 0.4)

kruskal.test(RootACF ~ Subcluster, data = tipacf.by.cluster) 
kr.p <- kruskal.test(RootACF ~ Subcluster, data = tipacf.by.cluster)$p.value
legend("topleft", legend = paste("Kruskal P =", signif(kr.p, 2)), cex= 1, bty = "n")





# --- EVENTFULNESS


# Generalise the above in a function (for other purposes):
# DataStack <- function(name.cat, x, y, z) {
#	coltitle = as.character(name.cat)
#	rep
#}


timepoint.counts <- apply(all.length, 2, Count, na.rm = TRUE)

	df1 <- data.frame(timepoint.counts[cluster1], rep("Cluster1", n1))
	df2 <- data.frame(timepoint.counts[cluster2], rep("Cluster2", n2))
	df3 <- data.frame(timepoint.counts[cluster3], rep("Cluster3", n3))
	col.names <- c("Value", "Source")
	colnames(df1) <- col.names
	colnames(df2) <- col.names
	colnames(df3) <- col.names
	rbind(df1, df2, df3)
	
	by.clust.timepoint.counts <- rbind(df1, df2, df3)
	


# setwd("/Users/Lab/Documents/Postdoc/ANALYSIS_local-files/ANALYSIS LOGS/2016-11_CCFs_Improvements")
# save.image("LastWorkspace_CCF_TOCA.Rdata")
# save.image("LastWorkspace_CCF_ENA.Rdata")
	
	
	
##### PLOTTING PROPERTIES PER (MANUAL) SUBCLUSTER:
	
library(RColorBrewer)
curr.cols <- brewer.pal(12, "YlGnBu")[c(7, 5, 1)]	

CurrData <- function(x, c1, c2, ...) {
  n1 <- Count(c1)
  n2 <- Count(c2)
  df1 <- data.frame(x[c1], rep("s8.3", n1))
  df2 <- data.frame(x[c2], rep("Other", n2))
  col.names <- c("Value", "Source")
  colnames(df1) <- col.names
  colnames(df2) <- col.names
  z <- rbind(df1, df2)
  # print(df1)
  # print(df2)
  # print(z)
  z
  # df1
}


CurrData <- function(x, c1, c2, ...) {
  n1 <- Count(c1)
  n2 <- Count(c2)
  df1 <- data.frame(x[c1], rep("Cluster 1", n1))
  df2 <- data.frame(x[c2], rep("Cluster 2", n2))
  col.names <- c("Value", "Source")
  colnames(df1) <- col.names
  colnames(df2) <- col.names
  z <- rbind(df1, df2)
  # print(df1)
  # print(df2)
  # print(z)
  z
  # df1
}


Boxplot <- function(col = c("red", "white", "grey"), ...) {
  # Boxplot <- function(col = c("grey80", "grey90", "grey70"), ...) {
  
  #dev.new()
  ylo = min(curr.data$Value, na.rm = TRUE)
  yhi = 1.1 * max(curr.data$Value, na.rm = TRUE)
  boxplot(Value ~ Source, data = curr.data,
          outpch = NA,
          col = col,
          #  main = curr.title,
          # ylab = curr.Ylab,
          #  ylim = c(ylo, yhi),
          ...
  )          
  stripchart(Value ~ Source, data = curr.data,
             add = TRUE,
             vertical = TRUE,
             method = "jitter",
             pch = 4, cex = 0.5,
             col = rgb(0,0,0,0.4)   
  )
  
  mw.p <- wilcox.test(Value ~ Source, curr.data)$p.value
  legend("topright", legend = paste("M-Wh P =", signif(mw.p, 2)), bty = "n")
}

dev.new(width = 8, height = 3)
par(mfrow = c(1,3))

# curr.data <- CurrData(max.lengths, s8.3, s8.others)
curr.data <- CurrData(metalist[[1]]$max.lengths, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = expression("Max length [" * mu * "m]"))
	
# curr.data <- CurrData(med.rate.extens/spt, s8.3, s8.others)
curr.data <- CurrData(metalist[[1]]$med.rate.extens/spt, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = expression("Med tip extension [" * mu * "m/s]"))

# curr.data <- CurrData(-med.rate.retract/spt, s8.3, s8.others)
curr.data <- CurrData(-metalist[[1]]$med.rate.retract/spt, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = expression("Med tip retraction [" * mu * "m/s]"))


# curr.data <- CurrData(straightness.mean, s8.3, s8.others)
curr.data <- CurrData(metalist[[1]]$straightness.mean, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = expression("Mean straightness [au]"))

# curr.data <- CurrData(med.fdcbm.invas/spt, s8.3, s8.others)
curr.data <- CurrData(metalist[[1]]$med.dB.invas/spt, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = expression("Med base invasion [" * mu * "m/s]"))

# curr.data <- CurrData(-med.fdcbm.retract/spt, s8.3, s8.others)
curr.data <- CurrData(-metalist[[1]]$med.dB.retract/spt, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = expression("Med base retraction [" * mu * "m/s]"))


# curr.data <- CurrData(all.time.ext, s8.3, s8.others)
curr.data <- CurrData(metalist[[1]]$all.time.ext, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = expression("Time extending"))
# curr.data <- CurrData(all.time.retr, s8.3, s8.others)
curr.data <- CurrData(metalist[[1]]$all.time.retr, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = expression("Time retracting"))
# curr.data <- CurrData(all.time.stall, s8.3, s8.others)
curr.data <- CurrData(metalist[[1]]$all.time.stall, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = expression("Time stalling"))

# curr.data <- CurrData(all.tip.nor1, s8.3, s8.others)
curr.data <- CurrData(metalist[[1]]$all.tip.nor1, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = expression("Norm tip fluorescence [au]"))


# curr.data <- CurrData(mean.tip.nor, s8.3, s8.others)
curr.data <- CurrData(metalist[[1]]$mean.tip.nor, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = expression("Mean norm tip fluorescence [au]"))

# mean.f.proj <- apply(all.proj, 2, mean, na.rm = TRUE)
# mean.f.body <- apply(all.body, 2, mean, na.rm = TRUE)
# mean.f.tip.raw <- apply(all.tip, 2, mean, na.rm = TRUE)

mean.f.proj <- apply(metalist[[1]]$all.proj, 2, mean, na.rm = TRUE)
mean.f.body <- apply(metalist[[1]]$all.body, 2, mean, na.rm = TRUE)
mean.f.tip.raw <- apply(metalist[[1]]$all.tip, 2, mean, na.rm = TRUE)



# curr.data <- CurrData(mean.f.proj, s8.3, s8.others)
curr.data <- CurrData(mean.f.proj, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = expression("Mean raw filo fluorescence [au]"))

# curr.data <- CurrData(mean.f.body, s8.3, s8.others)
curr.data <- CurrData(mean.f.body, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = expression("Mean raw body fluorescence [au]"))

# curr.data <- CurrData(mean.f.tip.raw, s8.3, s8.others)
curr.data <- CurrData(mean.f.tip.raw, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = expression("Mean raw tip fluorescence [au]"))
#

# PERSISTENCE OF TIP MOVEMENT:
#curr.data <- CurrData(acf.fdctm.roots*spt, s8.3, s8.others)
curr.data <- CurrData(metalist[[1]]$acf.dctm.crosspoints*spt, s2.1, s2.2) #!!!
Boxplot(curr.data, col = curr.cols, ylab = expression("Tip movement persistence [s]"))

objectnames

# Variance of tip movement
var.tipmove <- apply(all.move, 2, var, na.rm = TRUE)
var.tipmove <- apply(all.move, 2, var, na.rm = TRUE)
# curr.data <- CurrData(var.tipmove/spt, s8.3, s8.others)
curr.data <- CurrData(var.tipmove/spt, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = expression("Tip movement variance [" * mu * "m/s]"))

# Timepoints
timepoints <- apply(all.move, 2, Count)
# curr.data <- CurrData(timepoints, s8.3, s8.others)
curr.data <- CurrData(timepoints, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = expression("Timepoints"))

# PERSISTENCE OF TIP FLUORESCENCE:
# curr.data <- CurrData(acf.tip.f.roots*spt, s8.3, s8.others)
curr.data <- CurrData(acf.tip.f.roots*spt, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = expression("Tip F persistence [s]"))

# Variance of tip fluorescence
var.tipf <- apply(tip.f, 2, var, na.rm = TRUE)
curr.data <- CurrData(var.tipf, s8.3, s8.others)
curr.data <- CurrData(var.tipf, s2.1, s2.2)
Boxplot(curr.data, col = curr.cols, ylab = "Tip F variance")




