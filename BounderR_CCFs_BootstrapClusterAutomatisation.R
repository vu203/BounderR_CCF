"Script obsolete, replaced by BounderR_CCF_Randomisations.R"

## NEW R SCRIPT FOR DATASET RANDOMISATIONS AND AUTOMATED HIERARCHICAL CLUSTERING:


# QUESTION: 
# "A subset of filopodia in our dataset have a positive cross-correlation 
# between their tip fluorescence and tip movement. What is the likelihood that 
# this level of cross-correlation in a subcluster of the same size would be 
# observed by random if there was no relationship between fluorescence and 
# movement?"

# APPROACH:
# Take the original dataset and permute the order of observations for one of the 
# two parameters (tip movement) using a block reshuffling method (preserving the 
# distributions for the two parameters for each filopodium; blocks at least 
# roughly preserve auto-correlation within the data). Recluster the 
# randomised dataset, identify top correlating subcluster and record its cross-
# correlation properties. Repeat 500 times, then estimate the likelihood of getting 
# the observed cross-correlation between the two parameters by chance 
# from a randomly reshuffled dataset by comparing how many of the randomisations
# recapitulated the originally observed cross correlation.

# AIMS:

# 1: For a given clustered dataset (CCFs over lag for many filopodia), identify a subcluster 
# (monophyletic and within a specific range of element number) with the highest mean CCF value 
# at lag0, then record the mean CCF value at each lag of the "top subcluster" (if the
# top subcluster satisfies the condition of containing the correct number of elements,
# which must be similar as in the original dataset).

# 2: Loop the above operation across all randomisations

# 3: Analyse the CCFs of top correlating subclusters in all randomised datasets
# and compare to the top subcluster of the original dataset.


# ----------------------------------------------------------------------------
# Dependencies: 

# install.packages('purrr', dependencies= TRUE, repos='http://cran.rstudio.com/')
library(purrr)

# Requirements, from parent script ('BounderR_CCFs.R'):

# - Packages: 
library(RColorBrewer)

# - Data derived from the parent script ('BounderR_CCFs.R'):
#     dataframe 'all.move' (movement over time)
#     dataframe 'tip.f' (fluorescence over time)
#     dataframe 'ccf.tip.dctm' (cross-correlation table from the two dataframes above)

stopifnot(exists('all.move'))
stopifnot(exists('tip.f'))

# - Functions from parent script ('BounderR_CCFs.R'):
#     'Count', CI', 'ExtractBlockIndex', 'BlockReshuffle', 'DrawErrorAsPolygon'


# ----------------------------------------------------------------------------
# Original clustering functions (as in parent script):

GoCluster <- function(x, n.clusters) {
  map.input <- t(x)
  distance <- dist(map.input[, 18:24], method = "euclidean")
  cluster <- hclust(distance, method = "complete")
  cutree(cluster, k = n.clusters)
}
myHeatmap <- function(x) {
  map.input = t(x)
  distance <- dist(map.input[, 18:24], method = "euclidean")
  cluster <- hclust(distance, method = "complete")
  heatmap(map.input, Rowv = as.dendrogram(cluster), Colv = NA, xlab = "Lag", col = brewer.pal(12, "YlGnBu"), 	scale = "none")	
}
# example use: 
# myHeatmap(ccf.tip.dctm)

# nthSubcluster <- function(x, n.clusters, nth) {	
#   which(GoCluster(x, n.clusters = n.clusters) == nth)
# }
# nthSubclusterOthers <- function(x, n.clusters, nth) {
#   which(GoCluster(x, n.clusters = n.clusters) != nth)
# }

# ----------------------------------------------------------------------------
# FUNCTIONS FOR ITERATIVE SEARCH THROUGH NUMBERS OF CLUSTERS (until Top Cluster ~ Guide Size):

# Requirements: identify "top cluster" 
# loop k of cluster numbers until top cluster == guide size

ListSubclusters <- function(x, n.clusters = 4) {
  
  # Takes data frame x, outputs list of subtables as per specified number of clusters
  
  # Create a list of all subclusters:
    for (i in 1:n.clusters) {
      if (i == 1) {all.subcl = list()}
      cluster.i.index <- which(GoCluster(x, n.clusters) == i)
      all.subcl[[i]] <- x[, cluster.i.index]
    }
  return(all.subcl)
}

whichMax <-  function(x) which(x == max(x))

SizeMaxSubcluster <- function(x, n.clusters = 4) {
  
  # Takes data frame x, returns the size of the subcluster with highest mean at lag = 0 (at specified n clusters)
  #   - Calculate means at lag = 0
  #   - Find subcluster with highest mean at 0
  #   - Return the size of this top subcluster
  
  all.subcl <- ListSubclusters(x, n.clusters)
  clust.means <- lapply(all.subcl, function(x) apply(x["0", ], 1, mean, na.rm = TRUE)) # This trips single-element subclusters (number of dimensions)
  top.cluster <- whichMax(unlist(clust.means))
  top.size <- ncol(all.subcl[[top.cluster]])
  return(top.size)  
}

TopSubcluster <- function(x, n.clusters = 4) {
  
  # Takes data frame x, returns top most correlated cluster
  
  all.subcl <- ListSubclusters(x, n.clusters)
  clust.means <- lapply(all.subcl, function(x) apply(x["0", ], 1, mean, na.rm = TRUE)) # This trips single-element subclusters (number of dimensions)
  top.cluster <- whichMax(unlist(clust.means))
  output <- list(
    "top.cluster.ID" = top.cluster,
    "top.cluster"   = all.subcl[[top.cluster]]
  )
  return(output)
}

# e.g.
# ListSubclusters(ccf.tip.dctm, 7)[[7]]
# TopSubcluster(ccf.tip.dctm, 5)

SearchClusterSpace <- function(x, target.size = c(5,6,7)) {
  top.clust.size = NA
  initial.n.clusters = 1
  curr.n.clusters = initial.n.clusters
  while ((top.clust.size %in% target.size) == FALSE) {
    curr.n.clusters = curr.n.clusters + 1
    top.clust.size <- SizeMaxSubcluster(x, curr.n.clusters)
  }
  output <- list(
    "n.clusters" = curr.n.clusters,
    "top.cluster.ID" = TopSubcluster(x, curr.n.clusters)[["top.cluster.ID"]],
    "top.clust.size" = top.clust.size,
    "top.cluster" = TopSubcluster(x, curr.n.clusters)[["top.cluster"]]
  )
  return(output)
}


#------------------------------------------------------------------------------
#   BOOTSTRAPPING RANDOMISATIONS OF TIP MOVEMENT  
#------------------------------------------------------------------------------

### i.e. create 2500 randomised datasets not 1.
### (added to code on 8.3.2017) 


n.bootstrap = 2500  # <- Make this larger than the number of required clusters, not every 
                    #    dataset clusters well (allow 75% failure rate in the first instance)

#------------------------------------------------------------------------------
### 1. GENERATE RANDOMISATIONS:

# Create an empty randomisation table for each repeat (bootstrap) of randomisation,
# populate the table with results of BlockReshuffle on original data (filo by filo),
# and finally pass the resulting table for each boostrap to the list of randomisation tables.


CreateBootstrapRandomisations <- function(x, n.bootstrap) {
  
  for (boot in 1:n.bootstrap) {
    
    # Quality control (n.bootstrap 0 would cause failure):
    stopifnot(n.bootstrap >= 1)
    
    # In first run of the loop, create a list to store all randomisations 
    if (boot == 1) {list.all.rand <- list()}
    
    # Make new randomisation with each iteration of the loop:
    curr.randomisation <- matrix(NA, ncol = ncol(x), nrow = nrow(x))
    colnames(curr.randomisation) = colnames(x)
    for (i in seq_along(colnames(x))) {
      curr.randomisation[, i] <- BlockReshuffle(x[, i])
    }  
    
    # Add the new randomisation to the list of all: 
    list.all.rand[[boot]] <- curr.randomisation
    
    # Monitor loop progress:
    if (boot %% 25 == 0) {
      print(paste0(floor (boot/n.bootstrap * 100), "% completed"))
    }
   
  }
  stopifnot(length(list.all.rand) == n.bootstrap)
  return(list.all.rand)
}

list.all.move.rand <- CreateBootstrapRandomisations(all.move, n.bootstrap)

#------------------------------------------------------------------------------
### 2. CALCULATE CCF TABLES FOR ALL RANDOMISED DATASETS:

### Added to code ON 8.3.2017

list.ccf.randomised <- list()

for (boot in 1:n.bootstrap) {  
  curr.rand.move = list.all.move.rand[[boot]]
  
  # Create new table for CCFs of randomised dataset (using current randomisation only) 
  ccf.curr.randomisation = data.frame(matrix(NA, ncol = ncol(all.move), nrow = 2*maxlag + 1))
  colnames(ccf.curr.randomisation) = colnames(all.move)
  rownames(ccf.curr.randomisation) = lag.in.s
  
  stopifnot(n.filo == ncol(all.move))
  
  for (i in 1:n.filo) {
    
    # Create CCF for each filo, pass all into one table
    ccf.i <- ccf(tip.f[, i], curr.rand.move[, i], lag.max = 20, na.action = na.pass, plot = FALSE) 
    ccf.curr.randomisation[, i] <- ccf.i
    
    # Pass all CCFs for the current randomisation into a new list item
    list.ccf.randomised[[boot]] <- ccf.curr.randomisation
    rm(ccf.i)
  }
  rm(curr.randomisation)
 
  # Monitor loop progress:
  if (boot %% 25 == 0) {
    print(paste0(floor (boot/n.bootstrap * 100), "% completed"))
  } 
}

#------------------------------------------------------------------------------
### 3. CLUSTER ALL CCF TABLES  

# Run as many CCF tables through clustering as required to reach specified number of 
# successfully clustered randomised datasets

# 'new.list' stores result for all 2500 randomisations with NULL unsuccesful runs recorded; 
# 'boot.top.cluster.CCFs' stores results after cleaning up unsuccesful runs

n.bootstrap.clusters = 500  # <- HOW MANY CLUSTERS OF RANDOMISED DATASETS ARE REQUIRED?


#   This is the crux of the section! The loop underneath does the following: 
#   1. Creates an empty list to store new datasets on first iteration
#   2. Keeps count of 'successful' runs of the loop. A successful run means
#      that a randomised dataset could be clustered in such a way that the max
#      correlating cluster (top correlating subcluster, "TCS") had the required 
#      number of elements (as the original dataset, +- 1) (this number is specified 
#      in the SearchClusterSpace function as the 'target.size' argument).
#   3. While the count is less than specified (e.g. < 500 clusters), it keeps taking
#      new randomisations from the large pool of randomisations (list.ccf.randomised, 
#      e.g. 2500 randomisations), checks whether they can be 'successfully' clustered
#      - if yes, adds the result (of SearchClusterSpace) to the new list
#      - if no, adds a NULL list entry instead and moves to the next randomisation in
#      the pool.


for (i in 1:n.bootstrap) {
  if (i == 1) {new.list = list()}
  
  # Number of fruitful bootstrap clustering runs:
  n.count <- sum(unlist(lapply(new.list, function(x) !is.null(x))))
 
  # While count number is less than specified, keep running this loop:
  if (n.count < n.bootstrap.clusters) {
    try(SearchClusterSpace(list.ccf.randomised[[i]]), silent = TRUE)
    try(new.list[[i]] <- SearchClusterSpace(list.ccf.randomised[[i]]), silent = TRUE)
    n.count <- sum(unlist(lapply(new.list, function(x) !is.null(x))))
   # print(paste("n.count = ", n.count))
   # break
  } else {stop}
  # Monitor loop progress:
  if (n.count %% 10 == 0) {
    print(paste0(floor (n.count/n.bootstrap.clusters * 100), "% completed"))
  }
}

#------------------------ 
### A bit of cleaning up:

#new.list[[lapply(new.list, function (x) !is.null(x))]]

# non-null elements of list: 
notnull <- lapply(new.list, function (x) !is.null(x)) # (boolean vector)
which(notnull == TRUE) # (numeric vector)

# List of lists needs transposing
# Requirement: purrr package, function transpose()

new.list.non.null <- new.list[which(notnull == TRUE)]
new.list.non.null[1:10]
new.list.transposed <- transpose(new.list.non.null) # Transpose from PURRR package

boot.top.cluster.CCFs <- new.list.transposed[["top.cluster"]] # This works!

#------------------------- 
# Extract mean CCFs of the top correlating subcluster for each randomisation:

list.CCF.means <- lapply(boot.top.cluster.CCFs, function(x) apply(x, 1, mean, na.rm = TRUE))
list.CCF.ci <- lapply(boot.top.cluster.CCFs, function(x) apply(x, 1, CI))
df.CCF.means <- as.data.frame(list.CCF.means); colnames(df.CCF.means) <- as.character(1:500)
#df.CCF.CIs <- as.data.frame(list.CCF.ci); colnames(df.CCF.means) <- as.character(1:500)

# Mean and CI of CCFs of all 500 bootstrap repeats:
mean.CCF.boot <- apply(df.CCF.means, 1, mean, na.rm = TRUE) 
CI.CCF.boot <- apply(df.CCF.means, 1, CI) 

# Implement weighting according to length of time series?
# MeanByFiloWeight <- function(x) weighted.mean(x, w = weights.vec)
# list.CCF.w.means <- lapply(boot.top.cluster.CCFs, function(x) apply(x, 1, MeanByFiloWeight)) 
# Complicated at this point, weights for each subcluster will be different. Would need to extract 
# this upstream when creating the lists in the first place! 


#------------------------------------------------------------------------------
# 4. PLOT ALL RANDOMISATIONS:

# 4a) Line plot of top subcluster CCF (v offset) for all randomisations, 
#     compared with real dataset

range(list.CCF.means)

dev.new()
par(mfrow = c(1,1))
plot(NULL, 
     xlim = range(lag.in.s),
     ylim = 1.1*range(list.CCF.means),
     xlab = "Offset [s]",
     ylab = "CCF")

stopifnot(length(list.CCF.means) > 0)
for (i in 1:length(list.CCF.means)) {
  print(i)
  lines(lag.in.s, list.CCF.means[[i]],
        #col = "grey")
        col = "#AAAAAA30")
}

# Add the top correlating subcluster from actual data:

  real.top.cluster <- TopSubcluster(ccf.tip.dctm, 5)["top.cluster"]
  real.top.cluster <- as.data.frame(real.top.cluster)
  mean.ccf.real.top.cluster <- apply(real.top.cluster, 1, mean, na.rm = TRUE)
  
  lines(lag.in.s, mean.ccf.real.top.cluster,
      col = "#225EA8",
      lwd = 3)
  text(-40, 0.6, "Observed", col = "#225EA8", pos = 4)
  text(-40, 0.5, "Randomised", col = "#AAAAAA", pos = 4)
 
  abline(h = 0); abline(v = 0, lty = 3)

# After computing the representative randomisation (see 4c below), to plot its 
# line separately on the CCF line plot, use:
  
  # lines(lag.in.s, list.CCF.means[[representative]], col = "black", lty = 3)
  
# OPTION for alternative plot: Add CI (instead of many lines)
  
  ci1 = mean.CCF.boot + CI.CCF.boot
  ci2 = mean.CCF.boot - CI.CCF.boot
  
  DrawErrorAsPolygon(lag.in.s, ci1, ci2, tt = 0:41, col = "white")
  
  lines(lag.in.s, mean.CCF.boot,
        lwd = 1,
        col = "black")  
  

# 4b) Plot the distribution of simulations (as histogram) at offset = 0

means.at.lag.0 <- as.data.frame(list.CCF.means)[21, ]
colnames(means.at.lag.0) <- as.character(1:500)
table(means.at.lag.0)
  
hist(unlist(means.at.lag.0), breaks = 50, plot = TRUE,
     col = "grey", border = "white",
     xlim = c(0,0.7),
     main = "Mean CCF of top subcluster at offset = 0",
     xlab = "Mean CCF"
     )
abline(v = mean.ccf.real.top.cluster["0"],
       col = "red", lty = 2)
text(0, 29, "Randomised", pos = 4, col = "grey")
text(0, 26, "Observed", pos = 4, col = "red")

# How many randomisations have subclusters with a higher CCF at offset 0 than observed dataset?
Count(which(means.at.lag.0 >= mean.ccf.real.top.cluster["0"])) # Answer (ENA): 2 out of 500 (equivalent to p = 0.004)

# Just for the record: print top CCFs
means.at.lag.0[order(means.at.lag.0, decreasing = TRUE)][1:5]
mean.ccf.real.top.cluster["0"]


# 4c) Choose a representative example from 500 bootstraps & show its heatmap:

# The mean at offset 0 is:
mean.at.0 <- mean(unlist(means.at.lag.0), na.rm = TRUE)
distance.from.mean <- abs(means.at.lag.0 - mean.at.0)
representative <- which(distance.from.mean == min(distance.from.mean))

# This is not a trivial plotting task, see next section!

# -----------------------------------------------------------------------------
# 5. HEATMAP OF A REPESENTATIVE RANDOMISATION:

# First need to trace back where in the null-included list (new.list) the representative
# randomisation originally came from.

#-----------------
# 5a) DATA CLEANUP: Identify the full randomised dataset corresponding to the representative 
# top-correlating subcluster

# Question: if our hit is no.452 in new.list.non.null, which no. does it have in the bigger new.list?
# method: count which is the 452th occurence of non-NULL
# which(notnull == TRUE) is your friend

# inputs: 
# new.list.non.null # only non-NULL lists (only 500 randomisation that lead to fruitful clusters);
#                     this list contains top.cluster for each randomisation but not the full set of 45 filo!
# new.list          # Contains the full stuff (2500 randomisations), but not in clustered format

which(notnull == TRUE)
which(notnull == TRUE)[representative]
representative.match.in.big.list <- which(notnull == TRUE)[representative]
# Quality check:
stopifnot(identical(new.list[representative.match.in.big.list], new.list.non.null[representative]))

# From list.ccf.randomised pull the representative dataset
representative.randomisation <- as.data.frame(list.ccf.randomised[representative.match.in.big.list])

#-----------------
# 5a) Plot representative randomisation as heatmap:

map.input <- representative.randomisation
myHeatmap(map.input)

printEdges <- function(x) print(c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
printEdges(representative.randomisation)
leg = paste(c("min", "max"), signif(printEdges(map.input), 2))
mtext(leg[1], side=1, line=4, at=-40)
mtext(leg[2], side=1, line=3, at=-40)


# What proportion of randomisations result in succesful clustering?
# quick inaccurate estimate: 
representative/representative.match.in.big.list*100
# Proper way:
100 * (1 - n.bootstrap.clusters / length(new.list))


#------------------------------------------------------------------------------
# 6 Save workspace

save.image("CCF_500bootstraps.Rdata")

    
