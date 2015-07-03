# LIBRARIES ====
# http://cran.r-project.org/web/packages/biclust/biclust.pdf
library(biclust)  # Quite a recent package (May 12, 2015)

# LOAD DATA, REMOVE PDONSET ====
source('./preprocessing.R')
raw.nopdonset <- raw.filtered
raw.nopdonset.unscaled <- raw.filtered.unscaled
raw.nopdonset$pdonset <- NULL
raw.nopdonset.unscaled$pdonset <- NULL

# CONSTANTS ====
# Number of biclusters
# Note - N can't be more than 16 for par, will need to fix that later
# (multiple plots/pagination)
N <- 4
SAVE.PLOTS <- FALSE

# BICLUSTERING ====
# Bimax algorithm
bcl <- biclust(as.matrix(raw.nopdonset), method=BCPlaid(), iter.startup = 50, iter.layer = 50)
bcls <- bicluster(raw.nopdonset.unscaled, bcl)

# Transform pdonset to number of years w/ parkinsons?
# Then we want to look at some kind of comparison between NONMOTOR...
# MOTOR...
# and TIME to identify different TRAJECTORIES across time

# BUBBLEPLOT (?) ====
bubbleplot(raw.filtered, bcl)
if (SAVE.PLOTS) {
  dev.copy(pdf, paste("../figures/biclust-bubbleplot-", N, ".pdf", sep=""))
  dev.off()
}

# CLUSTER HEATMAPS ====

# Optimal grid organization
par(mfrow = n2mfrow(N))
for (i in 1:N) {
  drawHeatmap(raw.filtered, bicResult=bcl, number=i)
}
if (SAVE.PLOTS) {
  dev.copy(pdf, paste("../figures/biclust-heatmaps-", N, ".pdf", sep=""))
  dev.off()
}