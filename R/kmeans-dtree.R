# ==== LOAD LIBRARIES ====
library(cluster)
library(rpart)
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

# ==== NBCLUST ESTIMATION FOR OPTIMAL K (30 metrics) ====
# This is taking a long time
# nb <- NbClust(raw.filtered, distance = "euclidean",
#               min.nc=2, max.nc=15, method = "kmeans",
#               index = "alllong", alphaBeale = 0.1)
# hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,])))

# ==== PAM ESTIMATION FOR OPTIMAL K ====
# Estimates 2
pamk.best <- pamk(raw.filtered)
cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
plot(pam(raw.filtered, pamk.best$nc), main="Optimal clustering")

# ==== INITIAL KMEANS CLUSTERING ====
splits = splitdf(raw.filtered)
str(splits)
# NOTE: We're not splitting yet! (Don't know what the learning task is!)
# so override this and use everything!
splits$trainset = raw.filtered

# FIXME: How are we identifying k = 4? (Arbitrary right now)
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
