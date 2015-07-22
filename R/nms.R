# This contains further exploration of the nms-dominated subgroup

# Libraries ====
library(plyr)
library(ggplot2)

# Preprocessing ====
# At this point not necessary since I'm doing ad hoc exploration
# source('./preprocessing.R')

# Constants ====
SAVE.C1PLOTS <- FALSE


# Do kmeans clustering on subgroup ====
# c1 is what we're interested in, that's set in kmeans dtree
set.seed(0)
c1.woclass <- c1
c1.woclass$class <- NULL
c1.woclass.scaled <- as.data.frame(scale(c1.woclass))
c1.clus <- lapply(2:4, function(i) kmeans(c1.woclass.scaled, i, nstart = 25))
names(c1.clus) <- 2:4

# Bind cluster to orig data
c1.labeled <- lapply(c("2", "3", "4"), function(i) {
  cbind(c1.woclass, cluster = c1.clus[[i]]$cluster)
})
names(c1.labeled) <- 2:4

# Convert wide -> long ====
c1.labeled.long <- lapply(c("2", "3", "4"), function(i) {
  gather(c1.labeled[[i]], variable, measurement, age:axial)
})
names(c1.labeled.long) <- c("2", "3", "4")

# REORDER FACTORS BY INCREASING CISITOT ====
for (i in c("2", "3", "4")) {
  c1.cisitot <- c1.labeled.long[[i]][c1.labeled.long[[i]]$variable == "cisitot", ]
  c1.cisitot.means <- sapply(factor(1:as.integer(i)), function(i) {
    mean(c1.cisitot[c1.cisitot$cluster == i, ]$measurement)
  })
  c1.labeled.long[[i]]$cluster <- factor(c1.labeled.long[[i]]$cluster,
                                         levels = order(c1.cisitot.means))
}

# Plot boxplots of clustering ====
for (i in c("2", "3", "4")) {
  clus <- c1.labeled.long[[i]]
  p <- ggplot(clus, aes(x = factor(cluster), y = measurement, fill = factor(cluster))) +
    geom_boxplot() +
    geom_jitter() + # For now.
    guides(fill = FALSE) +
    facet_wrap( ~ variable, scales = "free")
  print(p)
  if (SAVE.C1PLOTS) {
    ggsave(paste("../figures/c1-summaries-", i, ".pdf", sep=""))
  }
}


# Try to figure out what kind of motor domination takes place ====
