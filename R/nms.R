# This contains further exploration of the nms-dominated subgroup

# Libraries ====
library(plyr)
library(ggplot2)

# Preprocessing ====
# At this point not necessary since I'm doing ad hoc exploration
# source('./preprocessing.R')

# Constants ====
SAVE.C4PLOTS <- FALSE


# Do kmeans clustering on subgroup ====
# c4 is what we're interested in, that's set in kmeans dtree
set.seed(0)
c4.woclass <- c4
c4.woclass$class <- NULL
c4.woclass.scaled <- as.data.frame(scale(c4.woclass))
c4.clus <- lapply(2:4, function(i) kmeans(c4.woclass.scaled, i, nstart = 25))
names(c4.clus) <- 2:4

# Bind cluster to orig data
c4.labeled <- lapply(c("2", "3", "4"), function(i) {
  cbind(c4.woclass, cluster = c4.clus[[i]]$cluster)
})
names(c4.labeled) <- 2:4

# Convert wide -> long ====
c4.labeled.long <- lapply(c("2", "3", "4"), function(i) {
  gather(c4.labeled[[i]], variable, measurement, age:pigd)
})
names(c4.labeled.long) <- c("2", "3", "4")

# REORDER FACTORS BY INCREASING CISITOT ====
for (i in c("2", "3", "4")) {
  c4.cisitot <- c4.labeled.long[[i]][c4.labeled.long[[i]]$variable == "cisitot", ]
  c4.cisitot.means <- sapply(factor(1:as.integer(i)), function(i) {
    mean(c4.cisitot[c4.cisitot$cluster == i, ]$measurement)
  })
  c4.labeled.long[[i]]$cluster <- factor(c4.labeled.long[[i]]$cluster,
                                         levels = order(c4.cisitot.means))
}

# Plot boxplots of clustering ====
for (i in c("2", "3", "4")) {
  clus <- c4.labeled.long[[i]]
  p <- ggplot(clus, aes(x = factor(cluster), y = measurement, fill = factor(cluster))) +
    geom_boxplot() +
    geom_jitter() + # For now.
    guides(fill = FALSE) +
    facet_wrap( ~ variable, scales = "free")
  print(p)
  if (SAVE.C4PLOTS) {
    ggsave(paste("../figures/c4-summaries-", i, ".pdf", sep=""))
  }
}


# Try to figure out what kind of motor domination takes place ====
