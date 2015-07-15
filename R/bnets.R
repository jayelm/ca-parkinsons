# Libraries ====
library(bnlearn)
library(gRain) # Note - depends on biocLite graph and RBGL packages

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
