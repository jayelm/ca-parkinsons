# LIBRARIES ====
# http://cran.r-project.org/web/packages/biclust/biclust.pdf
library(biclust)  # Quite a recent package (May 12, 2015)

# LOAD DATA ====
source('./preprocessing.R')

# CONSTANTS ====
# Number of biclusters
# Note - N can't be more than 16 for par, will need to fix that later
# (multiple plots/pagination)
N <- 16
SAVE.PLOTS <- FALSE

# BICLUSTERING ====
# Bimax algorithm
bcl <- biclust(as.matrix(raw.filtered), method=BCBimax(), minr=2, minc=2, number=N)
if (SAVE.PLOTS) {
  dev.copy(pdf, paste("../figures/biclust-", N, ".pdf", sep=""))
  dev.off()
}

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
# if (SAVE.PLOTS) {
#   dev.copy(pdf, paste("../figures/biclust-heatmap-", N, ".pdf", sep=""))
#   dev.off()
# }