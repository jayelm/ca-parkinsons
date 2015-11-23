# This contains further exploration of the nms-dominated subgroup
# added 20151022: it's c2 now

# Libraries ====
library(plyr)
library(ggplot2)
library(mclust)
library(cluster)
library(tidyr)

# Preprocessing ====
# At this point not necessary since I'm doing ad hoc exploration
# source('./preprocessing.R')
green.palette = c("#E6FAE6", "#80E680", "#00CD00", "#006600")

# Constants ====
SAVE.C2PLOTS <- TRUE

# VISUALIZE WSS ERROR TO FIND OPTIMAL K ====
# NOTE: Doesn't work well, there isn't any elbow

# Just in case it's still defined from below
var <- NULL
wss <- (nrow(c2.woclass)-1)*sum(apply(c2.woclass, 2, var))
for (i in 2:15) {
  wss[i] <- sum(kmeans(c2.woclass, i, nstart = 25)$withinss)
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
d_clust <- Mclust(as.matrix(c2.woclass), G=1:20)
m.best <- dim(d_clust$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")
# 3 clusters
# plot(d_clust) #, this is interactive use 0 to quit

# NBCLUST ESTIMATION FOR OPTIMAL K (30 metrics) ====
# This is taking a long time
# nb <- NbClust(c2.woclass, distance = "euclidean",
#               min.nc=2, max.nc=15, method = "kmeans",
#               index = "alllong", alphaBeale = 0.1)
# hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,])))

# PAM ESTIMATION FOR OPTIMAL K ====
# Estimates 2
pam_sils <- c()
for (i in 1:15) {
  pam_sils <- c(pam_sils, pam(c2.woclass, k=i)$silinfo$avg.width)
}
plot(x=1:14, y=pam_sils, xlab="Clusters", ylab="Average Silhouette Width", type="b")
if (SAVE.EXPLORE.PLOTS) {
  dev.copy(pdf, "../figures/asw.pdf")
  dev.off()
}

library(fpc)
pamk.best <- pamk(c2.woclass)
cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
plot(pam(c2.woclass, pamk.best$nc))

# GAP STATISTIC ESTIMATION ====
gaps <- clusGap(c2.woclass, kmeans, 14, B = 100)
plot(x=1:14, y=gaps$Tab[, "gap"], xlab="Clusters", ylab="Gap Statitsic", type="b")
if (SAVE.EXPLORE.PLOTS) {
  dev.copy(pdf, "../figures/gap-statistic.pdf")
  dev.off()
}

# Affinity propagation ====
# Auto commented out because this takes a while and throws error
# library(apcluster)
# c2.woclass.apclus <- apcluster(negDistMat(r=2), c2.woclass,
#           details = TRUE, q = 0, lam = 0.98,
#           maxits=10000, convits=1000)
# cat("affinity propogation optimal number of clusters:",
#     length(c2.woclass.apclus@clusters), "\n")
# # 4
# heatmap(c2.woclass.apclus)
# plot(c2.woclass.apclus, c2.woclass)

plot(1:15, wss, type="b",
     xlab="Number of Clusters", ylab="Within groups sum of squares")
if (SAVE.EXPLORE.PLOTS) {
  dev.copy(pdf, "../figures/kmeans-wss-error.pdf")
  dev.off()
}


# Do kmeans clustering on subgroup ====
# c2 is what we're interested in, that's set in exploration.R
set.seed(0)
c2.woclass <- c2
c2.woclass$class <- NULL
c2.woclass.scaled <- as.data.frame(scale(c2.woclass))
c2.clus <- lapply(2:4, function(i) kmeans(c2.woclass.scaled, i, nstart = 25))
names(c2.clus) <- 2:4

# Bind cluster to orig data
c2.labeled <- lapply(c("2", "3", "4"), function(i) {
  cbind(c2.woclass, cluster = c2.clus[[i]]$cluster)
})
names(c2.labeled) <- 2:4

# Convert wide -> long ====
c2.labeled.long <- lapply(c("2", "3", "4"), function(i) {
  gather(c2.labeled[[i]], variable, measurement, age:axial)
})
names(c2.labeled.long) <- c("2", "3", "4")

# REORDER FACTORS BY INCREASING CISITOT ====
for (i in c("2", "3", "4")) {
  c2.cisitot <- c2.labeled.long[[i]][c2.labeled.long[[i]]$variable == "cisitot", ]
  c2.cisitot.means <- sapply(factor(1:as.integer(i)), function(i) {
    mean(c2.cisitot[c2.cisitot$cluster == i, ]$measurement)
  })
  c2.labeled.long[[i]]$cluster <- factor(c2.labeled.long[[i]]$cluster,
                                         levels = order(c2.cisitot.means))
}

# Plot boxplots of clustering ====
for (i in c("2", "3", "4")) {
  clus <- c2.labeled.long[[i]]
  p <- ggplot(clus, aes(x = factor(cluster), y = measurement, fill = factor(cluster))) +
    geom_boxplot() +
    # geom_jitter() + # For now.
    guides(fill = FALSE) +
    facet_wrap( ~ variable, scales = "free") +
    # Added 20151022: colors!
    scale_fill_brewer(palette=green.palette)
  print(p)
  if (SAVE.C2PLOTS) {
    ggsave(paste("../figures/c2-summaries-", i, ".pdf", sep=""))
  }
}

# PAM optimum silhouette width ====
library(fpc)
pamk.best <- pamk(c2.woclass)
cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
plot(pam(c2.woclass, pamk.best$nc))

# WSS scree plot ====
mydata <- c2.woclass
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# Try to figure out what kind of motor domination takes place ====

# Apclus ====
# Why did this happen twice???
# library(apcluster)
# d.apclus <- apcluster(negDistMat(r=2), c2.woclass)
# cat("affinity propogation optimal number of clusters:", length(d.apclus@clusters), "\n")
# # 4
# heatmap(d.apclus)
# plot(d.apclus, c2.woclass)


# write.csv(c2.woclass, '~/Downloads/nms-unscaled.csv')
