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

# Trim this ====
nmses <- sapply(1:9, function(i) paste('nms_d', i, sep=''))
extras <- c('age', 'durat_pd', 'cisitot')
raw.nms <- raw.nopdonset[, c(nmses, extras)]
raw.nms.unscaled <- raw.nopdonset.unscaled[, c(nmses, extras)]
raw.motor <- raw.nopdonset[, c('rigidity', 'axial', 'bradykin', 'tremor', 'pigd')]
raw.motor.unscaled <- raw.nopdonset.unscaled[, c('rigidity', 'axial', 'bradykin', 'tremor', 'pigd')]

# BICLUSTERING ====
# Bimax algorithm
set.seed(0)
bcl <- biclust(as.matrix(raw.nms.unscaled), method=BCCC())
bcls <- bicluster(raw.nms.unscaled, bcl)
for (i in 1:100) {
  cat(i, '\n')
  varname <- paste('Bicluster', i, sep='')
  bc <- bcls[[varname]]
  plot(bc)
  x <- readline()
  corrplot(cor(bc), method = 'number')
  x <- readline()
  if (x == 'q') break
}
plot(bcls$Bicluster7)
drawHeatmap(raw.nms.unscaled, bicResult = bcl, number = 3)
for (bc in bcls) {
  cat('dim:', dim(bc), '\n')
}

# TODO: CORRPLOT

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
