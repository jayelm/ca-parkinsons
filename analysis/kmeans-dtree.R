# LOAD LIBRARIES ====
library(cluster)
library(reshape2)
library(clValid)
library(rpart)
library(rpart.plot)
library(corrplot)
library(plyr)
library(fpc)
library(NbClust)
library(MASS)
library(scatterplot3d)
library(rgl)
library(car)
# library for TODO: BIC
library(mclust)
library(ggplot2)
library(tidyr)
library(FSelector)
# library(apclust)
library(RColorBrewer) # For colors for decision tree plots
my.palette = brewer.pal(8, "Pastel1")


# CONSTANTS ====
# Remember - PDFs can vary even if the (seeded) clusters don't
# So this should be false unless I've changed something about the
# kmeans analysis
# explore plots is the determining clusters plots
SAVE.EXPLORE.PLOTS <- FALSE
SAVE.DTREES <- FALSE
SAVE.OVA.DTREES <- FALSE
SAVE.BOXPLOTS <- FALSE
SAVE.ARFF <- FALSE
SAVE.2VA.DTREES <- FALSE
# TODO: Make k = 2, 3, 4 modifiable via constant

# LOAD DATA ====
source('./preprocessing.R')

rename.clusters <- function(cl) {
  # Does passing in the argument not pass by reference?
  # I'm just putting this (ugly) code everywhere else
  for (i in 1:length(cl$cluster)) {
    curr.cluster <- cl$cluster[[i]]

    if (curr.cluster == 1) { # nms
      cl$cluster[[i]] <- 2
    } else if (curr.cluster == 2) { # all severe
      cl$cluster[[i]] <- 4
    } else if (curr.cluster == 3) { # motor
      cl$cluster[[i]] <- 3
    } else if (curr.cluster == 4) { # mild
      cl$cluster[[i]] <- 1
    }

  }
}

# VISUALIZE WSS ERROR TO FIND OPTIMAL K ====
# NOTE: Doesn't work well, there isn't any elbow

# Just in case it's still defined from below
library(TeachingDemos)
seed <- char2seed("Mu")
set.seed(seed)
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
# d_clust <- Mclust(as.matrix(raw.filtered), G=1:20)
# m.best <- dim(d_clust$z)[2]
# cat("model-based optimal number of clusters:", m.best, "\n")
# 3 clusters
# plot(d_clust)

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
set.seed(21) # Randomly chosen seed, for reproducibility
gaps <- clusGap(raw.filtered, kmeans, 14, B = 500)
par(mar = c(4.3, 4.7, 0.5, 0.5), ps = 18)
plot(gaps, xlab=expression(k), ylab=expression(Gap(k)))
# Embellishments
bestK <- 4 # Figure this out with method Tibs2001SEmax
gapk <- as.numeric(gaps$Tab[bestK, "gap"])
gapk1.se <- as.numeric(gaps$Tab[bestK + 1, "gap"] - gaps$Tab[bestK + 1, "SE.sim"])
segments(x0 = 4, y0 = 0, x1 = 4, y1 = gapk, lty = 2)
segments(x0 = 0, y0 = gapk, x1 = 4, y1 = gapk, lty = 2)
segments(x0 = 5, y0 = 0, x1 = 5, y1 = gapk1.se, lty = 2)
segments(x0 = 0, y0 = gapk1.se, x1 = 5, y1 = gapk1.se, lty = 2)
axis(1, at = c(1:14))
# mtext(expression(paste(Gap(4), phantom("|"))), side = 2, at = c(0, gapk), las = 1)
# text(x = 5, y = gapk1.se, labels = expression(paste(Gap(5) - se[5], phantom("|"))), pos = 4)

if (SAVE.EXPLORE.PLOTS) {
  dev.copy(pdf, "../figures/gap-statistic.pdf", width = 8, height = 5)
  dev.off()
}
print(gaps, method = "Tibs2001SEmax", SE.factor = 1)


# Affinity propagation ====
# Auto commented out because this takes a while
# library(apcluster)
# raw.filtered.apclus <- apcluster(negDistMat(r=2), raw.filtered,
#           details = TRUE, q = 0, lam = 0.98,
#           maxits=10000, convits=1000)
# cat("affinity propogation optimal number of clusters:",
#     length(raw.filtered.apclus@clusters), "\n")
# # 4
# heatmap(raw.filtered.apclus)
# plot(raw.filtered.apclus, raw.filtered)

# INITIAL KMEANS CLUSTERING ====
splits = splitdf(raw.filtered)
# str(splits)
# NOTE: We're not splitting yet! (Don't know what the learning task is!)
# so override this and use everything!
splits$trainset = raw.filtered

# FIXME: How are we identifying k = 4? (Arbitrary right now)
# Arbitrary 4 choice - variable is later
set.seed(0)  # Important to let clusters be the same
cl <- kmeans(splits$trainset, 4, nstart = 25)
# Added 20151020 - and rename here
# Changed 20160613 - added scmmotcp, need to rename
# 3 (sever)
print (cl$cluster)
  for (i in 1:length(cl$cluster)) {
    curr.cluster <- cl$cluster[[i]]

    if (curr.cluster == 1) { # nonmotor
      cl$cluster[[i]] <- 2
    } else if (curr.cluster == 2) { # motor
      cl$cluster[[i]] <- 3
    } else if (curr.cluster == 3) { # severe
      cl$cluster[[i]] <- 4
    } else if (curr.cluster == 4) { # mild
      cl$cluster[[i]] <- 1
    } else {
      print("wat")
    }
  }
print(cl$cluster)

# Add cluster label to original data, rename
trainset.labeled <- cbind(splits$trainset, cl$cluster)
trainset.labeled <- rename(trainset.labeled, c("cl$cluster"="cluster"))
# Convert to factor
trainset.labeled$cluster <- as.factor(trainset.labeled$cluster)

# BUILD DECISION TREE (UNPRUNED) ====
t <- rpart(cluster ~ ., trainset.labeled)

# Print error statistics
printcp(t)
# Find resubstitution rate
pred.t <- table(predict(t, type="class"), trainset.labeled$cluster)
resub <- 1-sum(diag(pred.t))/sum(pred.t)
print(resub)

# PLOT RESULTS (UNPRUNED TREE) ====
plot(t, uniform=TRUE, main="Decision Tree")
text(t, use.n=TRUE, all=TRUE, cex=.8)

# PRUNE TREE ====
# first find minimum cp
cp <- t$cptable[which.min(t$cptable[,"xerror"]),"CP"]
print(cp)
pruned.t <- prune(t, cp = cp)
# printcp(pruned.t)

# PLOT RESULTS (PRUNED TREE) ====
plot(pruned.t, uniform=TRUE, main="Decision Tree (Pruned)")
text(pruned.t, use.n=TRUE, all=TRUE, cex=.8, pos=1)
# stopping point - pruned tree

# KMEANS CLUSTERING COMPRESSED FUNCTION ====
# Use for iteration!
kmeans.dtree <- function(data, data.unscaled, k, save = FALSE, seed = 0) {
  # Reproducibility!
  set.seed(seed)
  cl <- kmeans(splits$trainset, k, nstart = 25)
  # Visualize clustering with clusplot
  clusplot(raw.filtered, cl$cluster, main = paste("Silhouette plot k =", k))
  if (save) {
    # FIXME: This save boolean flag is SAVE.DTREES outside of this function
    # which is slightly ambiguous (these are silhouette plots!)
    dev.copy(pdf, paste('../figures/kmeans-silhouette-', k, '.pdf', sep=''))
    dev.off()
  }
  labeled.data <- cbind(data, cluster=cl$cluster)
  # Rename clusters (added 20151020, changed again 2016)
  if (k == 4) {
    for (i in 1:length(cl$cluster)) {
      curr.cluster <- cl$cluster[[i]]
  
      if (curr.cluster == 1) { # nonmotor
        cl$cluster[[i]] <- 2
      } else if (curr.cluster == 2) { # motor
        cl$cluster[[i]] <- 3
      } else if (curr.cluster == 3) { # severe
        cl$cluster[[i]] <- 4
      } else if (curr.cluster == 4) { # mild
        cl$cluster[[i]] <- 1
      } else {
        print("wat")
      }
    }
  }
  # Convert to factor
  labeled.data$cluster <- as.factor(labeled.data$cluster)

  # Do the same thing with unscaled data to compare
  labeled.data.unscaled <- cbind(data.unscaled, cluster=cl$cluster)
  labeled.data.unscaled$cluster <- as.factor(labeled.data.unscaled$cluster)

  t.unscaled <- rpart(cluster ~ ., labeled.data.unscaled)
  t <- rpart(cluster ~ ., labeled.data)

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

  # Do the same with the unscaled tree
  cp.unscaled <- t.unscaled$cptable[which.min(t.unscaled$cptable[,"xerror"]),"CP"]
  pruned.t.unscaled <- prune(t.unscaled, cp = cp)

  # Plot and save trees
  # Also look at fancyRpartPlot()
  # fancyRpartPlot(tree.2)
  # plot(t, uniform=TRUE, main="Decision Tree")
  # text(t, use.n=TRUE, all=TRUE, cex=.8)
  # ?prp for extra keyword args - this one displays class proportions
  prp(t, extra = 1, main = paste("Unpruned Tree, ", i, " clusters", sep=""),
      box.col = c('red', 'green', 'blue', 'pink')[t$frame$yval])
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

  # Plot the unscaled tree
  prp(pruned.t.unscaled, extra = 1,
      main = paste("UNSCALED Pruned Tree, ", i, " clusters", sep=""))
  if (save) {
    dev.copy(pdf, paste('../figures/dtree-kmeans-pruned-unscaled-', k, '.pdf', sep=''))
    dev.off()
  }

  list(
    "data" = labeled.data,
    "clustering" = cl,
    "unpruned.tree" = t,
    "xv.error.unpruned" = xv.error.unpruned,
    "pruned.tree" = pruned.t,
    "xv.error.pruned" = xv.error.pruned
  )
}

# KMEANS DTREE OBTAINING PLOTS ====
# FIXME: This can be vectorized with lapply
trees <- vector(mode = "list", length = 3)
names(trees) <- c("clusters2", "clusters3", "clusters4")
for (i in 2:4) {
  istr <- paste("clusters", i, sep="")
  # Remember - PDFs can vary even if the (seeded) clusters don't
  trees[[istr]] <-  kmeans.dtree(raw.filtered, raw.filtered.unscaled, i,
                                 save = SAVE.DTREES, seed = 0)
}

# PRINT GLOBAL TREE STATS ====
for (i in 2:4) {
  istr <- paste("clusters", i, sep="")
  t <- trees[[istr]]$pruned.tree
  cat("CLUSTERS: ", i, "\n", sep="")
  cat("================================\n")
  cat("Complexity Parameter: ", tail(t$cptable[, "CP"], n=1), "\n", sep="")
  cat("10-fold CV error: ", trees[[istr]]$xv.error.pruned, "\n", sep="")
  cat("Root node error: ", t$frame$dev[1L] / t$frame$n[1L], "\n", sep="")
}

# PRINT CLUSTER STATS ====
for (i in 2:4) {
  istr <- paste("clusters", i, sep="")
  cl <- trees[[istr]]$clustering
  cat("CLUSTERS: ", i, "\n", sep="")
  cat("================================\n")
  cat("Sizes:", cl$size, "\n", sep=" ")
  cat("WithinSS:", cl$withinss, "\n", sep=" ")
  cat("Sum WithinSS:", sum(cl$withinss), "\n", sep=" ")
}

# PRINT CENTERS ====
for (i in 2:4) {
  istr <- paste("clusters", i, sep="")
  cl <- trees[[istr]]$clustering
  cat("CLUSTERS: ", i, "\n", sep="")
  cat("================================\n")
  for (j in 1:i) {
    center <- cl$centers[j, ]
    cat("> ", j, sep="")
    print(center)
  }
}

# BIND CLUSTER TO ORIGINAL DATA ====
labeled <- lapply(1:3, function(i) cbind(raw.omitted, trees[[paste("clusters", i+1, sep="")]]$clustering$cluster))
lnames <- c("clusters2", "clusters3", "clusters4")
names(labeled) <- lnames
# I use 145 indexing because I'm not sure if tail() copies
for (name in lnames) {
  # A convoluted way to get the last element - tail was throwing an error
  names(labeled[[name]])[length(names(labeled[[name]]))] <- "cluster"
}

# COMBINE LABELS WITH ORIGINAL DATA ====
clusters <- lapply(2:4, function(i) {
  name <- paste("clusters", i, sep="")
  current <- labeled[[name]]
  separated <- lapply(1:i, function(c) {
    current[current$cluster == c, ]
  })
})
names(clusters) <- lnames

for (i in 2:4) {
  names(clusters[[paste("clusters", i, sep="")]]) <- 1:i
}

# COLLECTED CLUSTERS AND SUMMARY OF CLUSTER STATISTICS ====
# Get both the summary of the statistics and collect the raw clusters
# in the form clusters$clustersN$[[N]] for box plotting later on
# TODO: Find similarity among the variables "For the interpreter"
# Rounding constant
PRECISION = 3

# FIXME: The results from shortened list of interpreted variables below, which I've saved in
# the data folder (for now), are not automatically saved to the data or figures folders. If
# they're needed, we need to save them again.
# interpreted <- c("age", "sex", "pdonset", "durat_pd", "cisitot")
# NOTE: This is now set as a constant up above

# This is a global collection of clusters to be used
clusters.raw <- vector(mode = "list", length = 3)
names(clusters.raw) <- c("2", "3", "4")

for (i in 2:4) {
  name <- paste("clusters", i, sep="")
  current <- clusters[[name]]
  # Begin collecting clusters of csv
  summaries.csv <- c()
  clusters.csv <- c()
  for (clus in 1:i) {
    cname <- paste(clus)
    # WE ONLY CARE ABOUT THE INTERPRETED
    # FALSE - now we're going to just include everything
    current.filtered <- current[[cname]][, INTERPRETED]
    # Print out stuff first
    # NOTE: this format differs from what is written to csv (csv has cluster column)
    cat("k = ", i, ", ", "CLUSTER ", cname, '\n', sep="")
    filtered.description <- round(describe(current.filtered)[, c('mean', 'sd', 'min', 'max', 'se')],
                                  PRECISION)
    # Select a couple of statistics
    write.csv(filtered.description)
    # Attach variable rownames, number of cluster (c)
    filtered.description <- cbind(variable=rownames(filtered.description), filtered.description)
    filtered.description <- cbind(cluster=rep(clus, length(INTERPRETED)), filtered.description)
    # Then add to csv list
    summaries.csv <- rbind(summaries.csv, filtered.description)

    # Attach cluster id to all elements in the cluster
    filtered.with_clusters <- cbind(cluster=rep(clus, nrow(current.filtered)), current.filtered)
    clusters.csv <- rbind(clusters.csv, filtered.with_clusters)
  }
  # Write cluster summaries as csv to data folder
  # NOTE: This (should) stays the same because of the seeding, unlike the PDFs,
  # so we don't need a boolean flag here
  write.csv(summaries.csv, file = paste("../data/kmeans-summaries-", i, ".csv", sep=""),
            row.names = FALSE)
  write.csv(clusters.csv, file = paste("../data/kmeans-raw-", i, ".csv", sep=""),
            row.names = FALSE)
  clusters.raw[[paste(i)]] <- clusters.csv
}

clusters.raw.nmsonly <- lapply(c("2", "3", "4"), function(i) {
  clusters.raw[[i]][, c(NMS.30, 'cisitot', 'cluster')]
})
names(clusters.raw.nmsonly) <- c("2", "3", "4")

# CLUSTERS RAW WIDE -> LONG ====
# NOTE: This lapply conversion hasn't been tested fully yet!
clusters.raw.long <- lapply(c("2", "3", "4"), function(i) {
  gather(clusters.raw[[i]], variable, measurement, age:scmmotcp)
})
names(clusters.raw.long) <- c("2", "3", "4")
# ADDED 01/20/16: shouldn't be through axial now, should be through nms30
# But we also add cisitot so we can do some factor reordering
clusters.raw.long.nmsonly <- lapply(c("2", "3", "4"), function(i) {
  gather(clusters.raw.nmsonly[[i]], variable, measurement, nms1:cisitot)
})
names(clusters.raw.long.nmsonly) <- c("2", "3", "4")

# REORDER FACTORS BY INCREASING CISITOT ====
for (i in c("2", "3", "4")) {
  cisitot <- clusters.raw.long[[i]][clusters.raw.long[[i]]$variable == "cisitot", ]
  cisitot.means <- sapply(factor(1:as.integer(i)), function(i) {
    mean(cisitot[cisitot$cluster == i, ]$measurement)
  })
  clusters.raw.long[[i]]$cluster <- factor(clusters.raw.long[[i]]$cluster,
                                           levels = order(cisitot.means))
}

for (i in c("2", "3", "4")) {
  cisitot <- clusters.raw.long.nmsonly[[i]][clusters.raw.long.nmsonly[[i]]$variable == "cisitot", ]
  cisitot.means <- sapply(factor(1:as.integer(i)), function(i) {
    mean(cisitot[cisitot$cluster == i, ]$measurement)
  })
  clusters.raw.long.nmsonly[[i]]$cluster <- factor(clusters.raw.long.nmsonly[[i]]$cluster,
                                           levels = order(cisitot.means))
}

# PLOT CLUSTER RESULTS ====
# FIXME: By using standard numeric list names, I can't use summaries$2 or summaries[[2]],
# only summaries$"2" or summaries[["2"]]
# seems confusing and not standardized with the rest of the script

for (i in c("2", "3", "4")) {
  clus <- clusters.raw.long[[i]]
  p <- ggplot(clus, aes(x = factor(cluster), y = measurement, fill = factor(cluster))) +
    geom_boxplot() +
    guides(fill = FALSE) +
    facet_wrap( ~ variable, scales = "free")
  print(p)
  if (SAVE.BOXPLOTS) {
    ggsave(paste("../figuresekmeans-summaries-", i, ".pdf", sep=""))
  }
}

# PLOT NMSONLY PLOTS ====
# Try just k = 4
p <- ggplot(clusters.raw.long.nmsonly[["4"]],
            aes(x = factor(cluster), y = measurement, fill = factor(cluster))) +
  geom_boxplot() +
  guides(fill = FALSE) +
  facet_wrap( ~ variable, scales = "free")
print(p)


# Conclusions:
# Most prevalent differences seen along sleep/fatigue and mood/cognit
# Also interesting is nms23, nms24, *nms28??*
# nms24 especially is revelant - nocturia
# nms23 - frequent urination
# 28 - change in perceiving flavors?

# NMSD2 / NMSD3 plot ====
# Plot nmsd2 and nmsd3 against each other, first for all clusters,
# then cluster 2 only
# All clusters
p <- ggplot(clusters.raw[["4"]], aes(x=nms_d2, y=nms_d3, color=factor(cluster))) +
  geom_point(position=position_jitter(width=0.5, height=0.5))
print(p)

# Cluster 2 only
# temp2 because I don't know how many times I'll need to make these things
temp2 <- clusters.raw[["4"]][clusters.raw[["4"]]$cluster == 2, ]
p <- ggplot(temp2,
       aes(x=nms_d2, y=nms_d3)) +
  geom_point(position=position_jitter(width=0.5, height=0.5),
             colour = 'green3') +
  geom_smooth(method = 'lm')
print(p)

# 3d plot, just looking at nms28, for fun
with(temp2,
     scatter3d(nms_d2, nms28, nms_d3, surface=TRUE, grid=TRUE, residuals=FALSE)
)

# OTHER NMSONLY SUBSET PLOTS ====
# Visualizing individual domains with parallel coordinates and 3d plots,
# where dimensionality is low enough
# NOTE: This overwrites clus4.wide below, be careful
clus4.wide <- clusters.raw[["4"]][, c(NMS.30, 'cluster')]
clus4.wide.wocluster <- clus4.wide
clus4.wide.wocluster$cluster <- NULL

c2 <- clus4.wide[clus4.wide$cluster == 2, ]
c2.wocluster <- c2
c2.wocluster$cluster <- NULL

# 1 is black, 2 is red, 3 is green, 4 is blue
parcoord(clus4.wide.wocluster[, c('nms1', 'nms2')],
         col = as.numeric(clus4.wide$cluster))

# Cardiovascular
clus4.wide.cfactor <- clus4.wide
clus4.wide.cfactor$cluster <- as.factor(clus4.wide.cfactor$cluster)
ggplot(clus4.wide.cfactor[, c('cluster', 'nms1', 'nms2')],
       aes(x=nms1, y=nms2, color=cluster)) +
  geom_point(position=position_jitter(width=1, height=1), size=3)

# Hallucinations
with(clus4.wide.cfactor[, c('cluster', 'nms13', 'nms14', 'nms15')],
     scatter3d(nms13, nms14, nms15, surface=FALSE, ellipsoid = TRUE, grid = FALSE,
               groups = cluster)
)
dim(clus4.wide.cfactor[, c('cluster', 'nms13', 'nms14', 'nms15')])
ggplot(clus4.wide.cfactor[, c('cluster', 'nms1', 'nms2')],
       aes(x=nms1, y=nms2, color=cluster)) +
  geom_point(position=position_jitter(width=1, height=1), size=3)


# Looks like nms10 is the biggest difference - so let's construct a table
# of the differences we're looking for

# SLEEP/FATIGUE indivi results ====
sleep.fatigue <- clus4.wide[, c('nms3', 'nms4', 'nms5', 'nms6', 'cluster')]
disp.sf.results <- function() {
  cat('sleep/fatigue\n')
  cat('symp\tclus1\tclus2\tclus3\tclus4\t4 - 2\t2 - 3\n')
  for (symptom in c('nms3', 'nms4', 'nms5', 'nms6')) {
    cat(symptom)
    for (c in 1:4) {
      cat('\t',
          round(mean(sleep.fatigue[sleep.fatigue$cluster == c, symptom]), 3),
          sep='')
    }
    cat('\t',
      round(
        mean(sleep.fatigue[sleep.fatigue$cluster == 4, symptom]) -
        mean(sleep.fatigue[sleep.fatigue$cluster == 2, symptom]),
        3
      ), sep=''
    )
    cat('\t',
      round(
        mean(sleep.fatigue[sleep.fatigue$cluster == 2, symptom]) -
        mean(sleep.fatigue[sleep.fatigue$cluster == 3, symptom]),
        3
      ), sep=''
    )
    cat('\n')
  }
}
disp.sf.results()

# MOOD/COGNITION indivi results ====
mood.cognition <- clus4.wide[, c('nms7', 'nms8', 'nms9',
                                 'nms10', 'nms11', 'nms12', 'cluster')]
disp.mc.results <- function() {
  cat('mood/cognition\n')
  cat('symp\tclus1\tclus2\tclus3\tclus4\t4 - 2\t2 - 3\n')
  for (symptom in c('nms7', 'nms8', 'nms9', 'nms10', 'nms11', 'nms12')) {
    cat(symptom)
    for (c in 1:4) {
      cat('\t',
          round(mean(mood.cognition[mood.cognition$cluster == c, symptom]), 3),
          sep='')
    }
    # Display some other helpful statistics -
    # 2v4 - mean of clus4 minus mean of clus2. For especially pronounced
    # symptoms for nonmotor-dom, this will be low
    cat('\t',
      round(
        mean(mood.cognition[mood.cognition$cluster == 4, symptom]) -
        mean(mood.cognition[mood.cognition$cluster == 2, symptom]),
        3
      ), sep=''
    )
    # 2v3 - mean of clus2 minus mean of clus3. For pronounced nonmotor-dom,
    # this should be high. Chose clus3 over clus4 because clus3 (motor-dom) is
    # more progressed PD and thus has higher likelihood of having higher
    # nms
    cat('\t',
      round(
        mean(mood.cognition[mood.cognition$cluster == 2, symptom]) -
        mean(mood.cognition[mood.cognition$cluster == 3, symptom]),
        3
      ), sep=''
    )
    cat('\n')
  }
}
disp.mc.results()

# NMS30 CORRPLOTS ====
# All clusters
corrplot(cor(clus4.wide), method="ellipse", order="hclust")
# Cluster 2 only
c2.cor.wocluster <- clus4.wide[clus4.wide$cluster == 2,]
c2.cor.wocluster$cluster <- NULL
#  Ify you want nms.
# c2.cor.wocluster[, c('bradykin', 'rigidity', 'tremor', 'axial')] <- cor.nms.raw2
corrplot(cor(c2.cor.wocluster), method="ellipse", order="hclust", main="2")

# 1, 3, 4
c1.cor.wocluster <- clus4.wide[clus4.wide$cluster == 1,]
c1.cor.wocluster$cluster <- NULL
corrplot(cor(c1.cor.wocluster), method="ellipse", order="hclust", main="1")
c3.cor.wocluster <- clus4.wide[clus4.wide$cluster == 3,]
c3.cor.wocluster$cluster <- NULL
corrplot(cor(c3.cor.wocluster), method="ellipse", order="hclust", main="3")
c4.cor.wocluster <- clus4.wide[clus4.wide$cluster == 4,]
c4.cor.wocluster$cluster <- NULL
corrplot(cor(c4.cor.wocluster), method="ellipse", order="hclust", main="4")
# Regular correlation, again ====

# set up subsets, remove clusters
cor.nms.raw <- trainset.labeled
cor.nms.raw1 <- cor.nms.raw[trainset.labeled$cluster == 1, ]
cor.nms.raw1$cluster <- NULL
cor.nms.raw2 <- cor.nms.raw[trainset.labeled$cluster == 2, ]
cor.nms.raw2$cluster <- NULL
cor.nms.raw3 <- cor.nms.raw[trainset.labeled$cluster == 3, ]
cor.nms.raw3$cluster <- NULL
cor.nms.raw4 <- cor.nms.raw[trainset.labeled$cluster == 4, ]
cor.nms.raw4$cluster <- NULL
cor.nms.raw$cluster <- NULL

scatterplotMatrix(~nms_d1+nms_d2+nms_d3+nms_d4+nms_d5+nms_d6, cor.nms.raw2)
corrplot(cor(cor.nms.raw), method="ellipse", order="hclust", main="All")
corrplot(cor(cor.nms.raw1), method="ellipse", order="hclust", main="1")
corrplot(cor(cor.nms.raw2), method="ellipse", order="hclust", main="2")
corrplot(cor(cor.nms.raw3), method="ellipse", order="hclust", main="3")
corrplot(cor(cor.nms.raw4), method="ellipse", order="hclust", main="4")

corrplot(cor(cor.nms.raw2), method="ellipse", order="hclust", main="2")

# 2016 update - some better cluster validation ====
clres.stability <- clValid(obj = raw.filtered, nClust = 2:8, clMethods = "kmeans",
                 validation = "stability", maxitems=1000)
par(mfrow=c(2, 2), mar = c(4.5, 4.5, 4.5, 4.5), oma = c(0, 0, 0, 0))
plot(clres.stability, main = "", legend = FALSE, pch = 1)
if (SAVE.EXPLORE.PLOTS) {
  dev.copy(pdf, "../figures/stability-measures.pdf", width = 8, height = 6)
  dev.off()
}
optimalScores(clres.stability)
stab.measures <- as.data.frame(clres.stability@measures)
colnames(stab.measures) <- 2:8
stab.measures$measure <- rownames(stab.measures)
stab.measures.long <- melt(stab.measures, id.vars = c("measure"),
                           variable.name = "clusters",
                           value.name = c("value"))
ggplot(stab.measures.long) +
  geom_line(aes(x = clusters, y = value, group = measure, color = measure)) +
  theme_bw()
clres.internal <- clValid(obj = raw.filtered, nClust = 2:8, clMethods = "kmeans",
                 validation = "internal", maxitems=1000)
plot(clres.internal)
stab.measures
# Stability measures look good for 4, surprisingly
# http://bioinformatics.oxfordjournals.org/content/19/4/459.full.pdf

# GENDER STATS ====
# TODO: Calculate for all k
cluster.genders <- sapply(1:4, function(i) {
  clus4 <- clusters.raw.long[["4"]]
  mean(clus4[clus4$cluster == i & clus4$variable == "sex", ]$measurement)
})

# BIG VERSIONS OF CLUS4 ====
# ^^ 20160613 sidestep previous stuff. Now add bind to original data and call it clus4
clus4 <- cbind(raw.omitted, cluster = cl$cluster)
clus4.long <- gather(clus4, variable, measurement, age:scmmotcp)
# clus4.nmsonly <- clus4[, c(NMS.30, 'cluster')]
# Manually save these!
p <- ggplot(clus4, aes(x = factor(cluster), y = measurement, fill = factor(cluster))) +
  geom_boxplot() +
  guides(fill = FALSE) +
  facet_wrap( ~ variable, scales = "free")
print(p)

# Split by level of NMS impairment
# NOTE: 20151020: change upon changing cluster!!
low.nms <- clus4[clus4$cluster == 1 | clus4$cluster == 3, ]
high.nms <- clus4[clus4$cluster == 2 | clus4$cluster == 4, ]
p <- ggplot(low.nms, aes(x = factor(cluster), y = measurement, fill = factor(cluster))) +
  geom_boxplot() +
  guides(fill = FALSE) +
  facet_wrap( ~ variable, scales = "free")
print(p)
p <- ggplot(high.nms, aes(x = factor(cluster), y = measurement, fill = factor(cluster))) +
  geom_boxplot() +
  guides(fill = FALSE) +
  facet_wrap( ~ variable, scales = "free")
print(p)

# Significant differneces between population means ====
# Try axial (most demonstrative feature)
clus4.wide.st <- clusters.raw[["4"]]
# Why isn't this already a factor? Really confused
# Because it wasn't set in clusters.raw - if this is a bad thing lmk
clus4.wide.st$cluster <- as.factor(clus4.wide.st$cluster)
# Just NMS pls
clus4.wide.st <- clus4.wide.st[, ]
# Assuming 1st column is cluster (which it should be)
oneways <- lapply(colnames(clus4.wide.st[, -1]), function(col) {
  fm <- substitute(i ~ cluster, list(i = as.name(col)))
  oneway.test(fm, clus4.wide.st)
})
for (test in oneways) {
  if (test$p.value < 0.05) {
    cat('sig\n')
  } else {
    cat('INSIG:\n')
    cat(test$data.name, '\n')
  }
}

# Tukey's HSD test ====

tukeys <- lapply(colnames(clus4.wide.st[, -1]), function(col) {
  # Doesn't work the oneway way for some reason!
  fm <- as.formula(paste(col, '~ cluster'))
  TukeyHSD(aov(fm, clus4.wide.st))$cluster
})
names(tukeys) <- colnames(clus4.wide.st[, -1])

for (v in names(tukeys)) {
  test <- tukeys[[v]]
  # Check for nonsignificant, since there are more significant
  sigs <- test[test[, "p adj"] > 0.05, ]
  if (!identical(logical(0), as.logical(sigs))) {
    # Super hacky to figure out if null matrix without type error
    cat(v, ' insignificant differences', ':', '\n', sep='')
    if (class(sigs) == 'numeric') {  # If returned just a single vector, can't do anything
      # Print em all, don't know how to get around this
      print(test)
    }
    print(sigs)
  } else {
    cat(v, 'nothing', '\n')
  }
}

# Rank by most informative features ====
# Commented out because is throwing a java.lang weka error
# features.ranked <- information.gain(cluster ~ ., clus4.wide.st)
# Add to a separate column because rownames are annoying
# features.ranked$variable <- rownames(features.ranked)
# rownames(features.ranked) <- NULL
# info.gain.ranks <- features.ranked[with(features.ranked, order(-attr_importance)), c('variable', 'attr_importance')]
# print(info.gain.ranks)

# NOTES ====
# Note - with clustering, indeed the main feature (root of the tree) doesn't carry much
# information
# e.g. splits like
# 90/214/56/421
# 89/421/57/214
# 214/421/57/89

# but nms_d3 < 17.5 seems to be a main classifier here, eeven though it's not very
# practically informative
# also nms_d7

# WRITE TO ARFF FOR WEKA ====
clus4.wide <- clusters.raw[["4"]][, c(2:19, 1)]
# NOTE: Here is where the cluster/class renaming happens (.arff)
clus4.wide <- rename(clus4.wide, c("cluster" = "class"))
clus4.wide$class <- as.factor(clus4.wide$class)
head(clus4.wide)
# Rename to parkinsons because it makes more sense
parkinsons <- clus4.wide
if (assertthat::are_equal(ncol(clusters.raw[["4"]]), ncol(parkinsons)) && SAVE.ARFF) {
  foreign::write.arff(parkinsons, '../data/parkinsons-k4.arff')
}

# 1vA decision trees ====

# Check SAVE.OVA.DTREES constant up top!

# For each cluster
set.seed(0)
for (i in 1:4) {
  cat("Testing", i, "\n")
  trainset.ova <- trainset.labeled
  # Convert to as.numeric so I can put 0 in there
  trainset.ova$cluster <- as.numeric(trainset.ova$cluster)
  # added 20151020 - DON'T RENAME HERE (already renamed!!!)
  # rename.clusters(trainset.ova)  # Rename before doing weird things

  trainset.ova[trainset.ova$cluster != i, ]$cluster <- 0
  trainset.ova$cluster <- as.factor(trainset.ova$cluster)
  # Attach to unscaled,
  trainset.ova.unscaled <- cbind(raw.filtered.unscaled, cluster = trainset.ova$cluster)
  t.ova <- rpart(cluster ~ ., trainset.ova.unscaled)
  if (i == 2) {
    t.ova.2 <- t.ova
  }

  prp(t.ova, extra = 1, main = paste("Unpruned ", i, " vs all", sep=""),
      box.col = tail(my.palette, n = -1)[t.ova$frame$yval])
  # Prune it later
  cp.ova <- t.ova$cptable[which.min(t.ova$cptable[,"xerror"]),"CP"]
  printcp(t.ova)
  print(cp.ova)
  if (SAVE.OVA.DTREES) {
    dev.copy(pdf, paste('../figures/dtree-', i, 'va-unpruned.pdf', sep=''))
    dev.off()
  }
  pruned.t.ova <- prune(t.ova, cp = cp.ova)
  prp(pruned.t.ova, extra = 1, main = paste("Pruned ", i, " vs all", sep=""),
      box.col = tail(my.palette, n = -1)[pruned.t.ova$frame$yval])
  if (SAVE.OVA.DTREES) {
    dev.copy(pdf, paste('../figures/dtree-', i, 'va-pruned.pdf', sep=''))
    dev.off()
  }
}

# Visualization with parallel coordinates ====
library(MASS)  # Modern applied statistics with S

# If not done already in exploration, double check
c1 <- clus4.wide[clus4.wide$class == 1, ]
c2 <- clus4.wide[clus4.wide$class == 2, ]
c3 <- clus4.wide[clus4.wide$class == 3, ]
c4 <- clus4.wide[clus4.wide$class == 4, ]

# c1.woclass <- c1
# c1.woclass$class <- NULL
c2.woclass <- c2
c2.woclass$class <- NULL

# No, we care about c2 now

# Pointless, doesn't help.
# parcoord(c1.woclass, col = as.numeric(c1$class))

# 2 and 4 vs 0 decision trees ====
trainset.2va <- trainset.labeled
# Convert to as.numeric so I can put 0 in there
trainset.2va$cluster <- as.numeric(trainset.2va$cluster)
# added 20151020 - DON'T RENAME HERE (already renamed!!!)
# rename.clusters(trainset.2va)  # Rename before doing weird things

trainset.2va[trainset.2va$cluster != 2 & trainset.2va$cluster !=4 , ]$cluster <- 0
trainset.2va$cluster <- as.factor(trainset.2va$cluster)
# Attach to unscaled,
trainset.2va.unscaled <- cbind(raw.filtered.unscaled, cluster = trainset.2va$cluster)
t.2va <- rpart(cluster ~ ., trainset.2va.unscaled)

prp(t.2va, extra = 1, main = paste("Unpruned ", i, " vs all", sep=""),
    box.col = tail(my.palette, n = -1)[t.2va$frame$yval])
# Prune it later
cp.2va <- t.2va$cptable[which.min(t.2va$cptable[,"xerror"]),"CP"]
printcp(t.2va)
print(cp.2va)
if (SAVE.2VA.DTREES) {
  dev.copy(pdf, paste('../figures/dtree-2and4va-unpruned.pdf', sep=''))
  dev.off()
}
pruned.t.2va <- prune(t.2va, cp = cp.2va)
prp(pruned.t.2va, extra = 1, main = paste("Pruned 2 and 4 vs rest", sep=""),
    box.col = tail(my.palette, n = -1)[pruned.t.2va$frame$yval])
if (SAVE.2VA.DTREES) {
  dev.copy(pdf, paste('../figures/dtree-2and4va-pruned.pdf', sep=''))
  dev.off()
}


# 2 and 3 vs 0 dtrees ====
trainset.2va <- trainset.labeled
# Convert to as.numeric so I can put 0 in there
trainset.2va$cluster <- as.numeric(trainset.2va$cluster)
# added 20151020 - DON'T RENAME HERE (already renamed!!!)
# rename.clusters(trainset.2va)  # Rename before doing weird things

trainset.2va[trainset.2va$cluster != 2 & trainset.2va$cluster != 3 , ]$cluster <- 0
trainset.2va$cluster <- as.factor(trainset.2va$cluster)
# Attach to unscaled,
trainset.2va.unscaled <- cbind(raw.filtered.unscaled, cluster = trainset.2va$cluster)
t.2va <- rpart(cluster ~ ., trainset.2va.unscaled)

prp(t.2va, extra = 1, main = paste("Unpruned ", i, " vs all", sep=""),
    box.col = tail(my.palette, n = -1)[t.2va$frame$yval])
# Prune it later
cp.2va <- t.2va$cptable[which.min(t.2va$cptable[,"xerror"]),"CP"]
printcp(t.2va)
print(cp.2va)
if (SAVE.2VA.DTREES) {
  dev.copy(pdf, paste('../figures/dtree-2and3va-unpruned.pdf', sep=''))
  dev.off()
}
pruned.t.2va <- prune(t.2va, cp = cp.2va)
prp(pruned.t.2va, extra = 1, main = paste("Pruned 2 and 3 vs rest", sep=""),
    box.col = tail(my.palette, n = -1)[pruned.t.2va$frame$yval])
if (SAVE.2VA.DTREES) {
  dev.copy(pdf, paste('../figures/dtree-2and3va-pruned.pdf', sep=''))
  dev.off()
}
