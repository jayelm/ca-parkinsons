# ==== LOAD LIBRARIES ====
library(cluster)
library(rpart)
library(rpart.plot)
library(plyr)
library(fpc)
library(NbClust)

# ==== LOAD DATA ====
source('./preprocessing.R')

# ==== VISUALIZE WSS ERROR TO FIND OPTIMAL K ====
# NOTE: Doesn't work well, there isn't any elbow

wss <- (nrow(raw.filtered)-1)*sum(apply(raw.filtered, 2, var))
for (i in 2:15) {
  wss[i] <- sum(kmeans(raw.filtered, i)$withinss)
}

plot(1:15, wss, type="b",
     xlab="Number of Clusters", ylab="Within groups sum of squares")

# ==== TODO: BIC ====
# See http://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
library(mclust)
# NOW - does this have to do with EM/hierarchical models. not cluster-based?
# Run the function to see how many clusters
# it finds to be optimal, set it to search for
# at least 1 model and up 20.
# Traditionally this takes a while, so it's uncommented
# d_clust <- Mclust(as.matrix(raw.filtered), G=1:20)
# m.best <- dim(d_clust$z)[2]
# cat("model-based optimal number of clusters:", m.best, "\n")
# # 3 clusters
# plot(d_clust)

# ==== NBCLUST ESTIMATION FOR OPTIMAL K (30 metrics) ====
# This is taking a long time
# nb <- NbClust(raw.filtered, distance = "euclidean",
#               min.nc=2, max.nc=15, method = "kmeans",
#               index = "alllong", alphaBeale = 0.1)
# hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,])))

# ==== PAM ESTIMATION FOR OPTIMAL K ====
# Estimates 2
pam_sils <- c()
for (i in 1:15) {
  pam_sils <- c(pam_sils, pam(raw.filtered, k=i)$silinfo$avg.width)
}
plot(x=1:14, y=pam_sils, xlab="Clusters", ylab="Average Silhouette Width", type="b")

# ==== GAP STATISTIC ESTIMATION ====
gaps <- clusGap(raw.filtered, kmeans, 14, B = 100)
plot(x=1:14, y=gaps$Tab[, "gap"], xlab="Clusters", ylab="Gap Statitsic", type="b")

# ==== INITIAL KMEANS CLUSTERING ====
splits = splitdf(raw.filtered)
str(splits)
# NOTE: We're not splitting yet! (Don't know what the learning task is!)
# so override this and use everything!
splits$trainset = raw.filtered

# FIXME: How are we identifying k = 4? (Arbitrary right now)
# Arbitrary 4 choice - variable is later
cl <- kmeans(splits$trainset, 4, nstart = 25)

# Add cluster label to original data, rename
trainset.labeled <- cbind(splits$trainset, cl$cluster)
trainset.labeled <- rename(trainset.labeled, c("cl$cluster"="cluster"))
# Convert to factor
trainset.labeled$cluster <- as.factor(trainset.labeled$cluster)

# ==== BUILD DECISION TREE (UNPRUNED) ====
t <- rpart(cluster ~ ., trainset.labeled)

# Print error statistics
printcp(t)
# Find resubstitution rate
pred.t <- table(predict(t, type="class"), trainset.labeled$cluster)
resub <- 1-sum(diag(pred.t))/sum(pred.t)
print(resub)

# ==== PLOT RESULTS (UNPRUNED TREE) ====
plot(t, uniform=TRUE, main="Decision Tree")
text(t, use.n=TRUE, all=TRUE, cex=.8)

# ==== PRUNE TREE ====
# first find minimum cp
cp <- t$cptable[which.min(t$cptable[,"xerror"]),"CP"]
print(cp)
pruned.t <- prune(t, cp = cp)
# printcp(pruned.t)

# ==== PLOT RESULTS (PRUNED TREE) ====
plot(pruned.t, uniform=TRUE, main="Decision Tree (Pruned)")
text(pruned.t, use.n=TRUE, all=TRUE, cex=.8, pos=1)
# stopping point - pruned tree

# ==== KMEANS CLUSTERING COMPRESSED FUNCTION ====
# Use for iteration!
kmeans_dtree <- function(data, k, save = FALSE, seed = 911) {
  # Reproducibility!
  set.seed(seed)
  cl <- kmeans(splits$trainset, k, nstart = 25)
  labeled_data <- cbind(data, cl$cluster)
  # Add cluster label to original data, rename
  labeled_data <- rename(labeled_data, c("cl$cluster"="cluster"))
  # Convert to factor
  labeled_data$cluster <- as.factor(labeled_data$cluster)
  t <- rpart(cluster ~ ., labeled_data)

  # Find 10-fold CV error rate for UNPRUNED tree
  # Old resub rate method
  # pred.t <- table(predict(t, type="class"), trainset.labeled$cluster)
  # resub <- 1-sum(diag(pred.t))/sum(pred.t)

  # From printcp source, this is how root node error is calculated
  # cat("Root node error: ", format(t$frame$dev[1L], digits = digits),
  #     "/", t$frame$n[1L], " = ", format(t$frame$dev[1L]/t$frame$n[1L],
  #                                     digits = digits), "\n\n", sep = "")
  # Specifically the root node error is
  root.node.error.unpruned <- t$frame$dev[1L] / t$frame$n[1L]
  # Then, multiply by minimum xerror. We know this is the last element of the cp table,
  # Since that's how we pruned it
  xerror.min.unpruned <- tail(t$cptable[, "xerror"], n=1)
  xv.error.unpruned <- xerror.min.unpruned * root.node.error.unpruned


  # PRUNE, find complexity parameter associated with minimum cross-validated error
  # TODO: Maybe select within 1 stddev?
  cp <- t$cptable[which.min(t$cptable[,"xerror"]),"CP"]
  pruned.t <- prune(t, cp = cp)

  # Calculate 10-fold CV error for PRUNED tree
  root.node.error.pruned <- pruned.t$frame$dev[1L] / pruned.t$frame$n[1L]
  xerror.min.pruned <- tail(pruned.t$cptable[, "xerror"], n=1)
  xv.error.pruned <- xerror.min.pruned * root.node.error.pruned

  # Plot and save trees
  # Also look at fancyRpartPlot()
  # fancyRpartPlot(tree.2)
  # plot(t, uniform=TRUE, main="Decision Tree")
  # text(t, use.n=TRUE, all=TRUE, cex=.8)
  # ?prp for extra keyword args - this one displays class proportions
  prp(t, extra = 1, main = paste("Unpruned Tree, ", i, " clusters", sep=""))
  # Get usr coordinates to add text to bottom right
  # NOTE: doesn't work, passing for now.
  # usr <- par("usr")
  # text(usr[2], usr[3], paste("10-fold CV error: ", xv.error.unpruned, sep=""), adj=c(1, 0))
  if (save) {
    dev.copy(pdf, paste('../figures/dtree-kmeans-unpruned-', k, '.pdf', sep=''))
    dev.off()
  }

  prp(pruned.t, extra = 1, main = paste("Pruned Tree, ", i, " clusters", sep=""))
  # usr <- par("usr")
  # text(usr[2], usr[3], paste("10-fold CV error: ", xv.error.pruned, sep=""), adj=c(1, 0))
  # plot(pruned.t, uniform=TRUE, main="Decision Tree (Pruned)")
  # text(pruned.t, use.n=TRUE, all=TRUE, cex=.8)
  if (save) {
    dev.copy(pdf, paste('../figures/dtree-kmeans-pruned-', k, '.pdf', sep=''))
    dev.off()
  }

  list(
    "data" = labeled_data,
    "clustering" = cl,
    "unpruned.tree" = t,
    # These are the same, you idiot.
    "xv.error.unpruned" = xv.error.unpruned,
    "pruned.tree" = pruned.t,
    "xv.error.pruned" = xv.error.pruned
  )
}

# ==== KMEANS DTREE OBTAINING PLOTS ====
trees <- vector(mode = "list", length = 3)
names(trees) <- c("clusters2", "clusters3", "clusters4")
for (i in 2:4) {
  istr <- paste("clusters", i, sep="")
  trees[[istr]] <-  kmeans_dtree(raw.filtered, i, save = TRUE)
}

# ==== PRINT GLOBAL TREE STATS ====
for (i in 2:4) {
  istr <- paste("clusters", i, sep="")
  t <- trees[[istr]]$pruned.tree
  cat("CLUSTERS: ", i, "\n", sep="")
  cat("================================\n")
  cat("Complexity Parameter: ", tail(t$cptable[, "CP"], n=1), "\n", sep="")
  cat("10-fold CV error: ", trees[[istr]]$xv.error.pruned, "\n", sep="")
  cat("Root node error: ", t$frame$dev[1L] / t$frame$n[1L], "\n", sep="")
}

# ==== PRINT CLUSTER STATS ====
for (i in 2:4) {
  istr <- paste("clusters", i, sep="")
  cl <- trees[[istr]]$clustering
  cat("CLUSTERS: ", i, "\n", sep="")
  cat("================================\n")
  cat("WithinSS:", cl$withinss, "\n", sep=" ")
  cat("Sizes:", cl$size, "\n", sep=" ")
}
# ==== NOTES ====
# Note - with clustering, indeed the main feature (root of the tree) doesn't carry much
# information
# e.g. splits like
# 90/214/56/421
# 89/421/57/214
# 214/421/57/89
# but nms_d3 < 17.5 seems to be a main classifier here, eeven though it's not very
# practically informative
# also nms_d7
