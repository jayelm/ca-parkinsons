# LIBRARIES ====
library(apcluster)
library(tidyr)

# PREPROCESSING ====
source('./preprocessing.R')

# CONSTANTS ====
# TODO: Because this is used by this and kmeans-dtree
# it should probably go in as a constant in preprocessing.R
INTERPRETED <- c("age", "sex", "pdonset", "durat_pd", "cisitot",
                 "nms_d1", "nms_d2", "nms_d3", "nms_d4", "nms_d5",
                 "nms_d6", "nms_d7", "nms_d8", "nms_d9",
                 "tremor", "bradykin", "rigidity", "axial", "pigd")
# save.plots for now since I don't know what else to plot!
SAVE.PLOTS <- FALSE

# CLUSTERING ====
# Lambda - dampening factor - usually 0.9, higher takes longer to converge
# but lets it converge in cases where it doesn't
# q = 0 - smallest input preference, lowest number of clusters
raw.filtered.apclus <- apcluster(negDistMat(r=2), raw.filtered,
                                 details = TRUE, q = 0, lam = 0.98,
                                 maxits=10000, convits=1000)
nclust <- length(raw.filtered.apclus@clusters)
apclusters <- raw.filtered.apclus@clusters
cat("Found ", nclust, " clusters\n", sep="")

# BIND TO ORIGINAL DATA ====

unscaled.apclus <- raw.omitted[, INTERPRETED]
unscaled.apclus$cluster <- 0
for (i in 1:nclust) {
  unscaled.apclus[apclusters[[i]], ]$cluster <- i
}

# WIDE -> LONG ====
unscaled.apclus.long <- gather(unscaled.apclus, variable, measurement, age:pigd)

# FACTOR CLUSTERS, ORDER BY BY CISITOT ====
# FIXME: Is it bad that the clusters are only factored here rather than
# in the above code block?
cisitot <- unscaled.apclus.long[unscaled.apclus.long$variable == "cisitot", ]
cisitot.means <- sapply(factor(1:nclust), function(i) {
  mean(cisitot[cisitot$cluster == i, ]$measurement)
})
unscaled.apclus.long$cluster <- factor(unscaled.apclus.long$cluster,
                                       levels = order(cisitot.means))

# BOXPLOT SUMMARIES (BORROWED FROM KMEANS-DTREE) ====
# TODO: Might want to find some reuse between this and kmeans-dtree
# Reorder factors by increasing cisitot (that's the most indicative scale we have of severity)
p <- ggplot(unscaled.apclus.long,
            aes(x = factor(cluster), y = measurement, fill = factor(cluster))) +
  geom_boxplot() +
  guides(fill = FALSE) +
  facet_wrap( ~ variable, scales = "free")
print(p)
if (SAVE.PLOTS) {
  ggsave("../figures/ap-summaries.pdf")
}

# SILHOUETTE PLOT ====
# Get a vector as long as nrow(raw.filtered) with clusters like kmeans
apclust.vector <- integer(nrow(raw.filtered))
for (i in 1:nclust) {
  apclust.vector[apclusters[[i]]] <- i
}
clusplot(raw.filtered, apclust.vector, main = paste("AP Silhouette Plot k =", nclust))
# This is hideous.
if (SAVE.PLOTS) {
  dev.copy(pdf, "../figures/ap-silhouette.pdf")
  dev.off()
}

# GET SIZES OF CLUSTERS ====
cat("cluster,size\n")
for (i in 1:nclust) {
  cat(i, ",", length(apclusters[[i]]), "\n", sep="")
}

# Visualization
# heatmap(raw.filtered.apclus)
# plot(raw.filtered.apclus, raw.filtered)

# NOTES ====
# Try modifying similarity matrix to not consider differences of nonmotor between motor?
# 7 is interesting
