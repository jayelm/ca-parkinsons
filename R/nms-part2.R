# Preprocessing for nms-specific analysis.

# Libraries ====
library(plyr)
# For PCA
library(FactoMineR)
library(psych)
library(mclust)
library(tidyr)
library(ggplot2)

# GLOBAL CONSTANTS ====
SAVE.PREPROCESSING.PLOTS <- FALSE
DB_FILE <- "../data/DATABASE_NMS Burden levels_15-4-2012.csv"
ALL.SYMPTOMS <- c(
  "nms1", "nms2", "nms3", "nms4", "nms5", "nms6",
  "nms7", "nms8", "nms9", "nms10", "nms11", "nms12",
  "nms13", "nms14", "nms15", "nms16", "nms17", "nms18",
  "nms19", "nms20", "nms21", "nms22", "nms23", "nms24",
  "nms25", "nms26", "nms27", "nms28", "nms29", "nms30"
)

INTERPRETED <- c("age", "sex", "pdonset", "durat_pd", "cisitot", ALL.SYMPTOMS)
# Change this variable to change symptoms
SYMPTOMS.TO.USE <- ALL.SYMPTOMS

# UTILITY FUNCTIONS ====
splitdf <- function(dataframe, seed=NULL, trainfrac=0.7) {
  if (trainfrac<=0 | trainfrac>=1) stop("Training fraction must be between 0 and 1, not inclusive")
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, trunc(length(index)/(1/trainfrac)))
  trainset <- dataframe[trainindex, ]
  testset <- dataframe[-trainindex, ]
  list(trainset=trainset,testset=testset)
}

# IMPORT DATA, PREPROCESSING ====
raw <- read.csv(DB_FILE)
# Study is irrelevant
raw$study <- NULL
# Get rid of NAs in required vals (later may use some kind of missing value compensation method
raw.omitted <- raw[, INTERPRETED]
raw.omitted <- na.omit(raw.omitted)

# raw with only select nms and motor symptoms
raw.filtered <- raw.omitted[, SYMPTOMS.TO.USE]

# DESCRIPTIVE STATISTICS (BEFORE STANDARDIZATION) ====
# TODO - change raw.filtered to raw.omitted.filtered, or something like that
raw.omitted.stats <- describe(raw.omitted)
raw.filtered.stats <- describe(raw.filtered)

# STANDARDIZATION ====
raw.filtered.unscaled <- raw.filtered
raw.filtered <- as.data.frame(scale(raw.filtered))

# PCA (Note: not helpful for identifying specific factors) ====
par(mfrow = c(1, 2))
pca <- PCA(raw.filtered)
pca$eig
summary(pca)
if (SAVE.PREPROCESSING.PLOTS) {
  dev.copy(pdf, "../figures/pca.pdf")
  dev.off()
}
# pca$eig
# pca$var$coord
# head(pca$ind$coord)
par(mfrow = c(1, 1))
plot(x = 1:length(rownames(pca$eig)), y = (pca$eig$eigenvalue), type="b",
     xlab="Factor", ylab="Eigenvalue", xaxp=c(0, 19, 19))
if (SAVE.PREPROCESSING.PLOTS) {
  dev.copy(pdf, "../figures/pca-eigenvalues.pdf")
  dev.off()
}

# Elbow appears around 3 factors
# Reset par for the rest of the script

# Optimal cluster evaluation via wss and pamk.best and gap statistic gives two clusters

set.seed(0)
var <- NULL
wss <- (nrow(raw.filtered)-1)*sum(apply(raw.filtered, 2, var))
for (i in 2:15) {
  wss[i] <- sum(kmeans(raw.filtered, i, nstart = 25)$withinss)
}

plot(1:15, wss, type="b",
     xlab="Number of Clusters", ylab="Within groups sum of squares")
if (SAVE.EXPLORE.PLOTS) {
  dev.copy(pdf, "../figures/kmeans-wss-error.pdf")
  dev.off()
}

# TODO: BIC ====
# See http://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
# NOW - does this have to do with EM/hierarchical models. not cluster-based?
# Run the function to see how many clusters
# it finds to be optimal, set it to search for
# at least 1 model and up 20.
# Traditionally this takes a while, so it's uncommented
# Number - 14!!! (?)
# d_clust <- Mclust(as.matrix(raw.filtered), G=1:20)
# m.best <- dim(d_clust$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")
# 3 clusters
plot(d_clust)

# NBCLUST ESTIMATION FOR OPTIMAL K (30 metrics) ====
# This is taking a long time
# nb <- NbClust(raw.filtered, distance = "euclidean",
#               min.nc=2, max.nc=15, method = "kmeans",
#               index = "alllong", alphaBeale = 0.1)
# hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,])))

#NOTE - changed 20151020 - commented to speed up calc
# PAM ESTIMATION FOR OPTIMAL K ====
# Estimates 2
pam_sils <- c()
for (i in 1:15) {
  pam_sils <- c(pam_sils, pam(raw.filtered, k=i)$silinfo$avg.width)
}
plot(x=1:14, y=pam_sils, xlab="Clusters", ylab="Average Silhouette Width", type="b")
if (SAVE.EXPLORE.PLOTS) {
  dev.copy(pdf, "../figures/asw.pdf")
  dev.off()
}

library(fpc)
pamk.best <- pamk(raw.filtered)
cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
plot(pam(raw.filtered, pamk.best$nc))

# GAP STATISTIC ESTIMATION ====
gaps <- clusGap(raw.filtered, kmeans, 14, B = 100)
plot(x=1:14, y=gaps$Tab[, "gap"], xlab="Clusters", ylab="Gap Statitsic", type="b")
gaps
if (SAVE.EXPLORE.PLOTS) {
  dev.copy(pdf, "../figures/gap-statistic.pdf")
  dev.off()
}

# Optimum ASW and WSS give 2 clusters, gap statistic gives 4 - what about MClust?

# TODO - model-based number of clusters
# d_clust <- Mclust(as.matrix(raw.filtered), G=1:20)
# m.best <- dim(d_clust$z)[2]
# cat("model-based optimal number of clusters:", m.best, "\n")
# 3 clusters
# plot(d_clust)

# model-based optimal number of clusters - 14 ====
raw.em <- cbind(raw.filtered, "cluster"=d_clust$classification)
raw.em$cluster <- as.factor(raw.em$cluster)

# kmeans optimal number of clusters - 2? 4? ====
set.seed(0)
cl.kmeans <- kmeans(raw.filtered, 4, nstart = 25)
raw.kmeans <- cbind(raw.filtered, "cluster"=cl.kmeans$cluster)
raw.kmeans$cluster <- as.factor(raw.kmeans$cluster)
raw.kmeans.long <- gather(raw.kmeans, variable, measurement, nms1:nms30)

p <- ggplot(raw.kmeans.long, aes(x = factor(cluster), y = measurement, fill = factor(cluster))) +
  geom_boxplot() +
  guides(fill = FALSE) +
  facet_wrap( ~ variable, scales = "free")
print(p)
