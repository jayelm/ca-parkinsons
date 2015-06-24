# LOAD LIBRARIES ====
library(cluster)
library(rpart)
library(rpart.plot)
library(plyr)
library(fpc)
library(NbClust)
# library for TODO: BIC
library(mclust)
library(ggplot2)

# CONSTANTS ====
# Remember - PDFs can vary even if the (seeded) clusters don't
# So this should be false unless I've changed something about the
# kmeans analysis
SAVE.DTREES = FALSE
# TODO: Make k = 2, 3, 4 modifiable via constant

# LOAD DATA ====
source('./preprocessing.R')

# VISUALIZE WSS ERROR TO FIND OPTIMAL K ====
# NOTE: Doesn't work well, there isn't any elbow

wss <- (nrow(raw.filtered)-1)*sum(apply(raw.filtered, 2, var))
for (i in 2:15) {
  wss[i] <- sum(kmeans(raw.filtered, i)$withinss)
}

plot(1:15, wss, type="b",
     xlab="Number of Clusters", ylab="Within groups sum of squares")

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
# # 3 clusters
# plot(d_clust)

# NBCLUST ESTIMATION FOR OPTIMAL K (30 metrics) ====
# This is taking a long time
# nb <- NbClust(raw.filtered, distance = "euclidean",
#               min.nc=2, max.nc=15, method = "kmeans",
#               index = "alllong", alphaBeale = 0.1)
# hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,])))

# PAM ESTIMATION FOR OPTIMAL K ====
# Estimates 2
pam_sils <- c()
for (i in 1:15) {
  pam_sils <- c(pam_sils, pam(raw.filtered, k=i)$silinfo$avg.width)
}
plot(x=1:14, y=pam_sils, xlab="Clusters", ylab="Average Silhouette Width", type="b")

# GAP STATISTIC ESTIMATION ====
gaps <- clusGap(raw.filtered, kmeans, 14, B = 100)
plot(x=1:14, y=gaps$Tab[, "gap"], xlab="Clusters", ylab="Gap Statitsic", type="b")

# INITIAL KMEANS CLUSTERING ====
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

# KMEANS DTREE OBTAINING PLOTS ====
# FIXME: This can be vectorized with lapply
trees <- vector(mode = "list", length = 3)
names(trees) <- c("clusters2", "clusters3", "clusters4")
for (i in 2:4) {
  istr <- paste("clusters", i, sep="")
  # Remember - PDFs can vary even if the (seeded) clusters don't
  trees[[istr]] <-  kmeans_dtree(raw.filtered, i, save = SAVE.DTREES)
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
labeled <- vector(mode = "list", length = 3)
labeled <- lapply(1:3, function(i) cbind(raw.omitted, trees[[paste("clusters", i+1, sep="")]]$clustering$cluster))
lnames <- c("clusters2", "clusters3", "clusters4")
names(labeled) <- lnames
# I use 145 indexing because I'm not sure if tail() copies
for (name in lnames) {
  names(labeled[[name]])[145] <- "cluster"
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

# SUMMARY OF CLUSTER STATISTICS ====
# TODO: Find similarity among the variables "For the interpreter"
# Rounding constant
PRECISION = 3

interpreted <- c("age", "sex", "pdonset", "durat_pd", "cisitot")
for (i in 2:4) {
  name <- paste("clusters", i, sep="")
  current <- clusters[[name]]
  # Begin collecting clusters of csv
  clusters.csv <- c()
  for (c in 1:i) {
    cname <- paste(c)
    # WE ONLY CARE ABOUT THE INTERPRETED
    current.filtered <- current[[cname]][, interpreted]
    # Print out stuff first
    # NOTE: this format differs from what is written to csv (csv has cluster column)
    cat("k = ", i, ", ", "CLUSTER ", cname, '\n', sep="")
    filtered.description <- round(describe(current.filtered)[, c('mean', 'sd', 'min', 'max')],
                                  PRECISION)
    # Select a couple of statistics
    write.csv(filtered.description)
    # Attach varaible rownames, number of cluster (c)
    filtered.description <- cbind(variable=rownames(filtered.description), filtered.description)
    filtered.description <- cbind(cluster=rep(c, length(interpreted)), filtered.description)
    # Then add to csv list
    clusters.csv <- rbind(clusters.csv, filtered.description)
  }
  # Write cluster summaries as csv to data folder
  # NOTE: This (should) stays the same because of the seeding, unlike the PDFs,
  # so we don't need a boolean flag here
  write.csv(clusters.csv, file = paste("../data/kmeans-summaries-", i, ".csv", sep=""),
            row.names = FALSE)
}

# PLOT CLUSTER RESULTS ====
# TODO: I may just need to plot everything, not just variables for the interpreter
# Then ideally (but not necessary) I'd want to find some way to differentiate the variables
# "For the interpreter" and the clustered variables and serialize that into the csv
summaries <- lapply(2:4, function(i) {
  raw <- read.csv(paste("../data/kmeans-summaries-", i, ".csv", sep=""),
                  stringsAsFactors = TRUE)
  # Turn numeric clusters into factors
  raw$cluster = factor(raw$cluster)
  raw
})
# FIXME: By using standard numeric list names, I can't use summaries$2 or summaries[[2]],
# only summaries$"2" or summaries[["2"]]
# seems confusing and not standardized with the rest of the script
names(summaries) <- 2:4 # NOTE: These numbers become factors (strings)!
for (k in c("2", "3", "4")) {
  summary <- summaries[[k]]
  p <- ggplot(summary, aes(x = cluster, y = mean, fill=cluster)) +
    guides(fill = FALSE) +
    facet_wrap( ~ variable, scales = "free") +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), position = position_dodge(), width = 0.2)
  print(p)
  ggsave(paste("../figures/kmeans-summaries-", k, ".pdf", sep=""))
}

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