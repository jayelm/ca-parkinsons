# ==== LOAD LIBRARIES ====
library(cluster)
library(rpart)
library(plyr)

# ==== UTILITY FUNCTIONS ====
splitdf <- function(dataframe, seed=NULL, trainfrac=0.7) {
  if (trainfrac<=0 | trainfrac>=1) stop("Training fraction must be between 0 and 1, not inclusive")
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, trunc(length(index)/(1/trainfrac)))
  trainset <- dataframe[trainindex, ]
  testset <- dataframe[-trainindex, ]
  list(trainset=trainset,testset=testset)
}

# ==== IMPORT DATA, PREPROCESSING ====
raw <- read.csv("../data/DATABASE_NMS Burden levels_15-4-2012.csv")
# Study is irrelevant
raw$study <- NULL
raw.omitted <- na.omit(raw)

# raw with only select nms and motor symptoms
raw.filtered <- raw.omitted[, c("nms_d1",
                                "nms_d2",
                                "nms_d3",
                                "nms_d4",
                                "nms_d5",
                                "nms_d6",
                                "nms_d7",
                                "nms_d8",
                                "nms_d9",
                                "tremor",
                                "bradykin",
                                "rigidity",
                                "axial",
                                "pigd")]

# raw.filtered <- raw.omitted[, c("tremor",
#                                 "bradykin",
#                                 "rigidity",
#                                 "axial",
#                                 "pigd")]

# ==== PLOTS ====

# cl_kmeans <- kmeans(raw.omitted, 4, nstart = 25)

wss <- (nrow(raw.filtered)-1)*sum(apply(raw.filtered, 2, var))
for (i in 2:15) {
  wss[i] <- sum(kmeans(raw.filtered, i)$withinss)
}

# http://cran.r-project.org/web/packages/NbClust/NbClust.pdf
# plot doesn't look good, no elbow
plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")


# ==== SANDBOX ====
splits = splitdf(raw.filtered)
str(splits)
# Interestingly enough this makes it worse. Then again I'm not sure if the splits
# were recalculated sooo...
# splits$trainset = raw.filtered

cl <- kmeans(splits$trainset, 4, nstart = 25)
# PCA -> 2d visuailze cl
# pc <- princomp(cl)
trainset.labeled <- cbind(splits$trainset, cl$cluster)
trainset.labeled <- rename(trainset.labeled, c("cl$cluster"="cluster"))
# Convert numerical cluster to factor
trainset.labeled$cluster <- as.factor(trainset.labeled$cluster)
t <- rpart(cluster ~ ., trainset.labeled)

# Print error statistics
printcp(t)
# Find resubstitution rate
pred.t <- table(predict(t, type="class"), trainset.labeled$cluster)
resub <- 1-sum(diag(pred.t))/sum(pred.t)
print(resub)
#
plot(t, uniform=TRUE, main="Decision Tree")
text(t, use.n=TRUE, all=TRUE, cex=.8)

# Prune the tree
# first find minimum cp
cp <- t$cptable[which.min(t$cptable[,"xerror"]),"CP"]
print(cp)
pruned.t <- prune(t, cp = cp)
# printcp(pruned.t)
plot(pruned.t, uniform=TRUE, main="Decision Tree (Pruned)")
text(pruned.t, use.n=TRUE, all=TRUE, cex=.8, pos=1)
# stopping point - pruned tree

# Note - with clustering, indeed the main feature (root of the tree) doesn't carry much
# information
# e.g. splits like
# 90/214/56/421
# 89/421/57/214
# 214/421/57/89
# but nms_d3 < 17.5 seems to be a main classifier here, eeven though it's not very
# practically informative
# also nms_d7
