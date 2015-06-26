# LIBRARIES ====
library(gclus)
library(ggplot2)
library(tidyr)
library(dynamicTreeCut)

# PREPROCESSING ====
source('./preprocessing.R')

# CONSTANTS ====
# TODO: Because this is used by this and kmeans-dtree
# it should probably go in as a constant in preprocessing.R
INTERPRETED <- c("age", "sex", "pdonset", "durat_pd", "cisitot",
                 "nms_d1", "nms_d2", "nms_d3", "nms_d4", "nms_d5",
                 "nms_d6", "nms_d7", "nms_d8", "nms_d9",
                 "tremor", "bradykin", "rigidity", "axial", "pigd")
SAVE.PLOTS = TRUE

# COMPUTE DISTANCE MATRIX (EUCLIDEAN) ====
hclust.dist <- dist(raw.filtered)

# HCLUST - MAX, AVERAGE, WARD, WARD.D2 ====
# Look at all 4 clustering methods side by side
par(mfrow=n2mfrow(4))

# Maximum (Complete Linkage)
hcl.complete <- hclust(hclust.dist)
plot(hcl.complete)

# Average
hcl.average <- hclust(hclust.dist, method = "average")
plot(hcl.average)

# Ward
hcl.ward.D <- hclust(hclust.dist, method = "ward.D")
plot(hcl.ward.D)

# Ward D2 (squared distance)
hcl.ward.D2 <- hclust(hclust.dist, method = "ward.D2")
plot(hcl.ward.D2)

if (SAVE.PLOTS) {
  dev.copy(pdf, "../figures/hc-dendrograms.pdf")
  dev.off()
}

# Reset
par(mfrow=c(1, 1))

# CUT TREES ====
# Let's pick Max for now

# 11 results in 4 clusters
HEIGHT = 11
K = 4

# Complete linkage cuts
# K = 4 cut
hcl.complete.mem <- cutree(hcl.complete, k = K)

# fixed height cut
# hcl.complete.mem <- cutree(hcl.complete, h = HEIGHT)
# dynamicTreeCut
# Be advised the default minimum cluster size is 20
# Hybrid (bottom up)
# K = 11
# hcl.complete.mem <- cutreeDynamic(hcl.complete, method = "hybrid",
#                                   distM = as.matrix(hclust.dist),
#                                   deepSplit = FALSE)
# Tree (top down)
# This one seems to create a horrible box plot
# hcl.complete.mem <- cutreeDynamic(hcl.complete, method = "tree",
#                                   deepSplit = FALSE)

# Ward (1963) linkage cuts
# hcl.complete.mem <- cutree(hcl.ward.D, h = 60)
# hcl.complete.mem <- cutree(hcl.ward.D, k = K)

# Then print out the number of clusters and frequencies
nclust <- length(unique(hcl.complete.mem))
cat("Number of clusters:", nclust, "\n")
print(count(hcl.complete.mem))

# Other optional methods
# hcl.average.mem <- cutree(hcl.average, k = 4)
# hcl.ward.D.mem <- cutree(hcl.ward.D, k = 4)
# hcl.ward.D2.mem <- cutree(hcl.ward.D2, k = 4)

# BIND CLUSTERS TO ORIGINAL DATA ====
labeled.complete <- cbind(cluster=hcl.complete.mem, raw.omitted[, INTERPRETED])

# WIDE -> LONG ====
labeled.complete.long <- gather(labeled.complete, variable, measurement, age:pigd)

# ORDER BY INCREASING CISITOT ====
# We need to shift by 1 if using dynamic, since dynamicTreeCut starts at 0.
DYNAMIC = FALSE
# DYNAMIC = TRUE
cisitot <- labeled.complete.long[labeled.complete.long$variable == "cisitot", ] 
anyNA(cisitot)
cisitot.means <- sapply(1:nclust, function(i) {
  mean(cisitot[cisitot$cluster == i - DYNAMIC, ]$measurement)
})
anyNA(cisitot.means)
labeled.complete.long$cluster <- factor(labeled.complete.long$cluster,
                                         levels = order(cisitot.means) - DYNAMIC)
anyNA(labeled.complete.long)

# VISUALIZE ====
p <- ggplot(labeled.complete.long,
            aes(x = factor(cluster), y = measurement, fill = factor(cluster))) +
  geom_boxplot() +
  guides(fill = FALSE) +
  facet_wrap( ~ variable, scales = "free")
print(p)

if (SAVE.PLOTS) {
  ggsave("../figures/hc-summaries-complete-k4.pdf")
  # ggsave("../figures/hc-summaries-ward-D-k4.pdf")
  # ggsave("../figures/hc-summaries-ward-D-h60.pdf")
  # ggsave("../figures/hc-summaries-complete-dynamic.pdf")
}