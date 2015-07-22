# Libraries ====
library(bnlearn)
library(gRain) # Note - depends on biocLite graph and RBGL packages

# Constants ====
SAVE.BNETS <- F

# Preprocessing ====
source('./preprocessing.R')
# Need stuff from kmeans dtree - so run that

clus4.numclus <- clus4.wide
# Convert this to numeric
clus4.numclus$cluster <- as.numeric(clus4.wide$cluster)
# I don't know why this is required
clus4.df <- as.data.frame(as.matrix(clus4.numclus))
head(clus4.df)

# Learning a bayesian network ====
clus4.iamb <- iamb(as.data.frame(as.matrix(clus4.numclus)))
# Tons of warnings
plot(clus4.iamb)
mb(clus4.iamb, 'cluster')
plot(hc(clus4.df))
mb(hc(clus4.df), 'cluster')

# Learning only in c3 ====
# These vars are set in exploration.R
# Standardize to force gaussian vars
c1.woclus <- c1
c1.woclus$class <- NULL
c1.std <- as.data.frame(scale(c1.woclus))
# Hill climbing ====
c1.hc <- bnlearn::hc(c1.std)
plot(c1.hc)
title('Hill-Climbing Algorithm on Cluster 1')
if (SAVE.BNETS) {
  dev.copy(pdf, '../figures/clus1-bnet-hc.pdf')
  dev.off()
}
# Grow shrink ====
c1.gs <- bnlearn::gs(c1.std)
plot(c1.gs)
title('Grow-Shrink Algorithm on Cluster 1')
if (SAVE.BNETS) {
  dev.copy(pdf, '../figures/clus1-bnet-gs.pdf')
  dev.off()
}
# Max-min hill-climbing ====
c1.gs <- bnlearn::mmhc(c1.std)
plot(c1.mmhc)
title('Min-Max Hill Climbing Algorithm on Cluster 1')
if (SAVE.BNETS) {
  dev.copy(pdf, '../figures/clus1-bnet-mmhc.pdf')
  dev.off()
}

mb(c1.hc, "cisitot")
