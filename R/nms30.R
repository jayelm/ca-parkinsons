# Performing column clustering on only nms30. like nms-part2.R but graduated
# Depends on everything.wide data frame used in longitudinal.R
# NMS.30, NMS.30.NAMES from preprocessing
library(plyr)
library(FactoMineR)
library(rpart)
library(rpart.plot)
library(nFactors)
library(gplots) # For heatmaps
# Better dendrograms
library(dendextend)
library(dendsort)
library(dynamicTreeCut)
# Trying to get some color here
library(colorspace)
library(sparcl)
library(cluster)
library(fpc)
library(mclust)
library(RColorBrewer)
library(fmsb)
library(proxy)
library(reshape2)
library(ggplot2)
library(tidyr)
library(DMwR)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
# library(apcluster)

SAVE.NMS30.PLOTS <- TRUE

# Preprocess ====
NMS.30.MAP <- setNames(NMS.30.NAMES, NMS.30)

nms30 <- everything.wide[, NMS.30]
nms30 <- rename(nms30, NMS.30.MAP)
# Standardize
nms30.s <- as.data.frame(scale(nms30))

# With motor
nms30m <- nms30
nms30m[, MOTOR.SYMPTOMS] <- raw.omitted[, MOTOR.SYMPTOMS]

nms30m.scale <- scale(nms30m)
nms30m.s <- as.data.frame(nms30m.scale)

# Transposes (right now, don't need)
# nms30.t <- as.data.frame(t(nms30))
# nms30.st <- as.data.frame(t(nms30.s))

# Dimensionality on this is ridiculous. Not sure it makes sense.

# Factor analysis ====
# FactoMineR PCA
pca <- PCA(nms30.s)
# Obviously depression symptoms are on their own
# Plot variance explained
plot(pca$eig$`cumulative percentage of variance`,
     main="Cumulative variance explained",
     ylab="% variance",
     xlab="Component"
     )

# Look at PCA contributions
# HELP: Not sure how to interpret PCA
# heatmap.2(pca$var$contrib, Colv=F, col=cm.colors, dendrogram='row', trace='none')
# heatmap.2(pca$var$cor, Colv=F, col=cm.colors, dendrogram='row', trace='none')
corrplot(t(pca$var$cor))
if (SAVE.NMS30.PLOTS) {
  dev.copy(pdf, '../figures/nms30only-eigscorr.pdf', width=8, height=5)
  dev.off()
}
# Determine number of eigenvalues
dev.off()
ev <- eigen(cor(nms30.s)) # get eigenvalues
ap <- parallel(subject=nrow(nms30.s),var=ncol(nms30.s),
               rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
# Plot a couple of eigenvalue metrics. Notice how most metrics
# say 6.
plotnScree(nS)
if (SAVE.NMS30.PLOTS) {
  dev.copy(pdf, '../figures/nms30only-eigs.pdf', width=7.5, height=5)
  dev.off()
}

# Find optimal # clusters ====
wss <- (nrow(nms30.st)-1)*sum(apply(nms30.st, 2, var))
for (i in 2:15) {
  wss[i] <- sum(kmeans(nms30.st, i, nstart = 25)$withinss)
}

# Large heatmap, clustering for rows and cols ====
hm <- heatmap.2(as.matrix(nms30.s),
                hclustfun = function(x) hclust(x, method='average'),
                trace='none')
colD <- as.ggdend(hm$colDendrogram)
rowD <- as.ggdend(hm$rowDendrogram)

# motor
# NOTE: apcluster interferes with this!!
hm.m <- heatmap.2(as.matrix(nms30m.s[, 1:35]),
                  hclustfun = function(x) hclust(x, method='average'),
                  trace='none')
colD.m <- as.ggdend(hm.m$colDendrogram)

# Average is different
asdf <- hclust(dist(t(nms30m.s)), method='average')
plot(asdf)

# Plot column dendrogram, visualize dynamic clusters with cutreeDynamic
ctd <- cutree(hm$colDendrogram, 5)
par(mfrow = c(1, 1))
ColorDendrogram(as.hclust(hm$colDendrogram), y = ctd, labels = names(ctd),
                main = "NMS Hierarchical Clustering",
                xlab = "Symptom",
                sub = "",
                branchlength = 4)
if (SAVE.NMS30.PLOTS) {
  dev.copy(pdf, '../figures/nms30only-colhc.pdf', width=10, height=6)
  dev.off()
}

# Motor
ctd.m <- cutree(hm.m$colDendrogram, 4)
ctk.m <- hm.m$colDendrogram %>% set("branches_k_color", k = 4)
ColorDendrogram(as.hclust(hm.m$colDendrogram), y = ctd.m, labels = names(ctd.m),
                main = "Nonmotor + Motor Hierarchical Clustering",
                xlab = "Symptom",
                sub = "",
                branchlength = 5)
if (SAVE.NMS30.PLOTS) {
  dev.copy(pdf, '../figures/nms30m-colhc.pdf', width=10, height=6)
  dev.off()
}

# A slightly better plot

# Individuals: Determining k ====
set.seed(0)

# wss - 2
wss <- (nrow(nms30.s)-1)*sum(apply(nms30.s,2,var))
set.seed(0)
for (i in 2:15) wss[i] <- sum(kmeans(nms30.s, centers=i, nstart=25)$withinss)
plot(1:15, wss, type='b', xlab='Number of clusters', ylab='wss')

# motor - 2
wss.m <- (nrow(nms30m.s)-1)*sum(apply(nms30m.s,2,var))
set.seed(0)
for (i in 2:15) wss.m[i] <- sum(kmeans(nms30m.s, centers=i, nstart=25)$withinss)
plot(1:15, wss.m, type='b', xlab='Number of clusters', ylab='wss')

# pamk - 2
pamk.best <- pamk(nms30m.s)
cat("number of clusters estimated by optimum average silhoeutte width:", pamk.best$nc, "\n")
plot(pam(nms30.s, pamk.best$nc))

# motor - 2
pamk.best.m <- pamk(nms30m.s)
cat("number of clusters estimated by optimum average silhoeutte width:", pamk.best.m$nc, "\n")
plot(pam(nms30.s, pamk.best.m$nc))

# mclust - probably will take a while (it does)
MCLUST_MAX = 10
d_clust <- Mclust(as.matrix(nms30.s), G=1:MCLUST_MAX)
m.best <- dim(d_clust$z)[2]
cat("model based optimal number of clusters:", m.best,"\n")
plot(d_clust)

# motor
MCLUST_MAX = 15
d_clust.m <- Mclust(as.matrix(nms30m.s), G=1:MCLUST_MAX)
# OR, the VEV model ONLY
d_clust.vev <- Mclust(as.matrix(nms30m.s), modelNames = c('VEV'), G=1:MCLUST_MAX)
# How about d_clust.vev at 4? Still high BIC but maybe more info
d_clust.vev <- Mclust(as.matrix(nms30m.s), modelNames = c('VEV'), G=4)
m.best.m <- dim(d_clust.m$z)[2]
cat("model based optimal number of clusters:", m.best.m, "\n")
plot(d_clust.m)
# notice that VEI - diagonal, varying volume, equal shape is the best, peaks at
# 11 clusters,
# but VEV - ellipsoidal, equal shape hits a decent max at 3 clusters.
# STRATEGY - plot parts of intresting 11 clusters in VEI case.
# VEV - display k = 3 case
# NOT INTERESTING!!! SO BORING....
# GETS MAXIMIZED WELL ENOUGH AT 3 OR 4


# library(cluster)
# useless
set.seed(21)
# set.seed(15)
cg <- clusGap(nms30.s, kmeans, 14, B = 500)
plot.new()
print(cg, method = "Tibs2001SEmax")
par(mar = c(4.3, 4.7, 0.5, 0.5), ps = 18)
plot(cg, xlab=expression(k), ylab=expression(Gap(k)), xaxt = 'n')
# Embellishments
bestK <- 6 # Figure this out with method Tibs2001SEmax
gapk <- as.numeric(cg$Tab[bestK, "gap"])
gapk1.se <- as.numeric(cg$Tab[bestK + 1, "gap"] - cg$Tab[bestK + 1, "SE.sim"])
segments(x0 = bestK, y0 = 0, x1 = bestK, y1 = gapk, lty = 2)
segments(x0 = 0, y0 = gapk, x1 = bestK, y1 = gapk, lty = 2)
segments(x0 = bestK + 1, y0 = 0, x1 = bestK + 1, y1 = gapk1.se, lty = 2)
segments(x0 = 0, y0 = gapk1.se, x1 = bestK + 1, y1 = gapk1.se, lty = 2)
axis(1, at = c(1:14))

dev.copy(pdf, "../figures/gap-statistic-6.pdf", width = 8, height = 5)
dev.off()

print(cg, method = "Tibs2001SEmax", SE.factor = 1)

# Testing
# set.seed(21)
set.seed(15)
testing <- kmeans(nms30.s, 7, nstart = 1000)
heatmap.2(t(testing$centers), dendrogram ='col', Rowv = FALSE, Colv = TRUE,
          col = rev(colorRampPalette(brewer.pal(11, "RdBu"))(1000)),
          trace = 'none')
mclust::adjustedRandIndex(testing$cluster, cl$cluster)
table(cl$cluster[testing$cluster == 1])
table(cl$cluster[testing$cluster == 6])
table(cl$cluster[testing$cluster == 5])
table(cl$cluster[testing$cluster == 3])
table(cl$cluster[testing$cluster == 2])
table(cl$cluster[testing$cluster == 4])
table(cl$cluster[testing$cluster == 2])

# Why are there 2s in this section:
weirdos <- nms30.s[names(cl$cluster[testing$cluster == 6] == 2), ]
heatmap.2(t(as.matrix(weirdos)), trace = 'none',
          col = rev(colorRampPalette(brewer.pal(11, "RdBu"))(1000)),
          dendrogram = 'col', Rowv = FALSE, Colv = TRUE
          )

# pam_sils <- c()
# for (i in 1:15) {
#   pam_sils <- c(pam_sils, pam(nms30.s, k=i)$silinfo$avg.width)
# }
# plot(x=1:14, y=pam_sils, xlab="Clusters", ylab="Average Silhouette Width", type="b")
# if (SAVE.NMS30.PLOTS) {
#   dev.copy(pdf, "../figures/nms30-asw.pdf")
#   dev.off()
# }

# cluster metrics are totally, absolutely, meaningless
clres.stability <- clValid(obj = nms30.s, nClust = 2:8, clMethods = "kmeans",
                 validation = "stability", maxitems=1000)
stab.measures <- as.data.frame(clres.stability@measures)
colnames(stab.measures) <- 2:8
stab.measures <- melt(t(stab.measures))
stab.measures <- rename(stab.measures, c('Var1' = 'n', 'Var2' = 'measure'))
ggplot(stab.measures, aes(x=n, y=value, color=measure, shape=measure)) +
  geom_point() +
  geom_line() +
  theme_bw()

clres.internal <- clValid(obj = nms30.s, nClust = 2:8, clMethods = "kmeans",
                 validation = "internal", maxitems=1000)
internal.measures <- as.data.frame(clres.internal@measures)
colnames(internal.measures) <- 2:8
internal.measures <- melt(t(internal.measures))
internal.measures <- rename(internal.measures, c('Var1' = 'n', 'Var2' = 'measure'))
ggplot(internal.measures, aes(x=n, y=value, color=measure, shape=measure)) +
  geom_point() +
  geom_line() +
  theme_bw()

# KMEANS ====

nms30.s <- nms30m.s # UNCHECK IF YOU WANT MOTOR
set.seed(0)
K = 11
cl <- kmeans(nms30.s, centers = K, nstart = 25)
# Reattach
nms30.k4 <- nms30
# nms30.k4$cluster <- cl$cluster
# MESS AROUND WITH THIS. BUT BEFORE YOU DO SO, make sure to set K appropriately
nms30.k4$cluster <- d_clust.m$classification
# Counts
table(nms30.k4$cluster)
nms30.k4[, ALL.BUT.NMS] <- everything.wide[, ALL.BUT.NMS]

# Wide -> Long and reorder factors
nms30.k4.long <- gather(nms30.k4, variable, measurement,
                        c(`d1-1-lightheaded`:`d9-30-sweating`, age:axial))

cisitot <- nms30.k4.long[nms30.k4.long$variable == 'cisitot', ]
cisitot.means <- sapply(factor(1:K), function(i) {
  mean(cisitot[cisitot$cluster == i, ]$measurement)
})
nms30.k4.long$cluster <- factor(nms30.k4.long$cluster,
                                levels = order(cisitot.means))

p <- ggplot(nms30.k4.long, aes(x = factor(cluster), y = measurement, fill = factor(cluster))) +
  geom_boxplot() +
  guides(fill = FALSE) +
  facet_wrap( ~ variable, scales = "free")
print(p)

# Get back unscaled centers ====
# Expression heatmap
# Order should be 3 (mild), 1 (avg? wut), 2 (depression), 4 (severe)

cl.centers <- cl$centers[match(c(3, 2, 1, 4), rownames(cl$centers)), ]
rownames(cl.centers) <- 1:4

# 9 domains, and motor
gch <- gg_color_hue(10)
heatmap.2(t(cl.centers), Rowv = FALSE, Colv = FALSE, dendrogram = 'none', trace = 'none',
          col = colorRampPalette(c('green', 'black', 'red'))(n = 1000),
          RowSideColors = c(rep(gch[1], 2), rep(gch[2], 4), rep(gch[3], 6), rep(gch[4], 3), rep(gch[5], 3),
                            rep(gch[6], 3), rep(gch[7], 3), rep(gch[8], 2), rep(gch[9], 4), rep(gch[10], 4)),
          xlab = 'Symptom', ylab = 'Subtype', key.xlab = 'z-score',
          cexCol = 2, cexRow = 1.5, srtCol = 0,
          margins = c(8, 24),
          density.info = 'none'
          )
if (SAVE.NMS30.PLOTS) {
  dev.copy(pdf, '../figures/nms30m-k4-heatmap.pdf', width=10, height=14)
  dev.off()
}

# Write to table
rounded.centers <- round(t(unscale(cl$centers, nms30m.scaleinfo)), 2)
write.table(x = rounded.centers,
            file = '../data/nms30m-k4-centers.csv',
            # Settings make setting up LaTeX table easier
            sep = ' & ', eol=' \\\\\n', quote = FALSE)

# If using mclust, then calculate means of each symptom
nms30.k4.s <- as.data.frame(scale(nms30.k4))
nms30.k4.s$cluster <- nms30.k4$cluster
mclust.means <- as.data.frame(t(sapply(1:11, function(i) {
  apply(nms30.k4.s[nms30.k4.s$cluster == i, ], 2, mean)
})))
# Clean up
mclust.means$cluster <- NULL
mclust.means$sex <- NULL
# Order by cisitot (ish) since it's too much work to do otherwise
mclust.means.counts <- as.numeric(table(d_clust.m$classification))[order(mclust.means$cisitot)]
mclust.means <- mclust.means[match(order(mclust.means$cisitot), rownames(mclust.means)), ]
# Scale
rownames(mclust.means) <- 1:11
mclust.means[mclust.means > 2] <- 2
mclust.means[mclust.means < -2] <- -2
# mclust.means.s <- scale(mclust.means)

heatmap.2(t(mclust.means), Rowv = FALSE, Colv = FALSE, dendrogram = 'none', trace = 'none',
          # breaks = c(seq(-2, -0.6, length = 10), seq(-0.6, 0.6, length = 10), seq(0.6, 2, length = 10)),
          col = colorRampPalette(c('green', 'black', 'red'))(n = 1000), symkey = FALSE,
#           RowSideColors = c(rep(gch[1], 2), rep(gch[2], 4), rep(gch[3], 6), rep(gch[4], 3), rep(gch[5], 3),
#                             rep(gch[6], 3), rep(gch[7], 3), rep(gch[8], 2), rep(gch[9], 4), rep(gch[10], 4)),
          xlab = 'Symptom', ylab = 'Subtype', key.xlab = 'z-score',
          cexCol = 1.4, cexRow = 1.4, srtCol = 0,
          margins = c(8, 24),
          density.info = 'none'
          )

# Decision trees assuming cl cluster
nms30.dtree <- nms30.k4
nms30.dtree$cluster <- factor(cl$cluster)
# Wait, rename these??
# Remember it's 3, 2, 1, 4
# Nums are 507 97 249 49
nms30.dtree$cluster <- revalue(nms30.dtree$cluster, c('3'='1', '2'='2', '1'='3', '4'='4'))
nms30.dtree$cluster <- factor(nms30.dtree$cluster, levels = c('1', '2', '3', '4'))
tree.nms30 <- rpart(cluster ~ ., nms30.dtree)
cp <- tree.nms30$cptable[which.min(tree.nms30$cptable[,"xerror"]),"CP"]
pruned.tree.nms30 <- prune(tree.nms30, cp = cp)
prp(pruned.tree.nms30, extra = 1, varlen=0,
    box.col = gg_color_hue(4)[pruned.tree.nms30$frame$yval])
if (SAVE.NMS30.PLOTS) {
  dev.copy(pdf, '../figures/nms30m-dtree-all.pdf', width=12, height=10)
  dev.off()
}

# one vs all ====
gg_color_list <- function(i) {
  c(gg_color_hue(i)[i], '#FFFFFF')
}
ova.counts <- c(509, 97, 249, 49)
for (ONE in 1:4) {
  set.seed(0)
  nms30.ova <- nms30.dtree
  nms30.ova$newcluster <- NA
  nms30.ova[nms30.ova$cluster == ONE, ]$newcluster <- T
  nms30.ova[nms30.ova$cluster != ONE, ]$newcluster <- F
  this.count = length(nms30.ova[nms30.ova$newcluster == T, ]$newcluster)
  if (this.count != ova.counts[ONE]) stop("Lengths not equal")
  nms30.ova$cluster <- factor(nms30.ova$newcluster, levels=c(T, F))
  nms30.ova$newcluster <- NULL
  tree.nms30.ova <- rpart(cluster ~ ., nms30.ova)
  cp <- tree.nms30.ova$cptable[which.min(tree.nms30.ova$cptable[,"xerror"]),"CP"]
  if (ONE == 4 || ONE == 2) {
    pruned.tree.nms30.ova <- tree.nms30.ova
  } else {
    pruned.tree.nms30.ova <- prune(tree.nms30.ova, cp = cp)
  }
  prp(pruned.tree.nms30.ova, extra = 1, varlen=0,
      box.col = gg_color_list(ONE)[pruned.tree.nms30.ova$frame$yval],
      main = paste(ONE, ' vs all decision tree', sep = ''))
  if (SAVE.NMS30.PLOTS) {
    dev.copy(pdf, paste('../figures/nms30m-dtree-', ONE, 'va.pdf', sep=''), width=10, height=10)
    dev.off()
  }
}

# Mahalnobis distnace, hc, AP ====
nms30.s.dist <- dist(nms30.s, method = 'Mahalanobis')
nms30.s.nDM<- negDistMat(nms30.s, r=2)
hc <- hclust(nms30.s.dist)
plot(hc)
ctd <- cutreeDynamic(hc, distM = as.matrix(nms30.s.dist), method = 'hybrid')
ctd.k4 <- cutree(hc, 4)


# Use Mahalanobis distance to do affinity propagation
apc <- apcluster(as.matrix(nms30.s.dist), details=TRUE)
# Plot exemplars
nms30.apc <- nms30
nms30.apc$cluster <- NA
for (i in 1:length(apc@clusters)) {
  nms30.apc[apc@clusters[[i]], ]$cluster <- i
}

nms30.apc[, ALL.BUT.NMS] <- everything.wide[, ALL.BUT.NMS]

# Wide -> Long and reorder factors
nms30.apc.long <- gather(nms30.apc, variable, measurement,
                        c(`d1-1-lightheaded`:`d9-30-sweating`, age:axial))

cisitot <- nms30.apc.long[nms30.apc.long$variable == 'cisitot', ]
cisitot.means <- sapply(factor(1:K), function(i) {
  mean(cisitot[cisitot$cluster == i, ]$measurement)
})
nms30.apc.long$cluster <- factor(nms30.apc.long$cluster,
                                levels = order(cisitot.means))

p <- ggplot(nms30.apc.long, aes(x = factor(cluster), y = measurement, fill = factor(cluster))) +
  geom_boxplot() +
  guides(fill = FALSE) +
  facet_wrap( ~ variable, scales = "free")
print(p)

# Tukey's HSD
oneways <- lapply(colnames(nms30.k4[, -31]), function(col) {
  fm <- substitute(i ~ cluster, list(i = as.name(col)))
  oneway.test(fm, nms30.k4)
})
for (test in oneways) {
  if (test$p.value < 0.05) {
    cat('sig\n')
  } else {
    cat('INSIG:\n')
    cat(test$data.name, '\n')
  }
}

# Factor clusters

# ova.counts <- c(509, 97, 249, 49)
nms30.k4.hsd <- nms30.k4
nms30.k4.hsd$cluster <- factor(nms30.k4.hsd$cluster)
nms30.k4.hsd$cluster <- revalue(nms30.k4.hsd$cluster, c('3'='1', '2'='2', '1'='3', '4'='4'))
nms30.k4.hsd$cluster <- factor(nms30.k4.hsd$cluster, levels = c('1', '2', '3', '4'))
table(nms30.k4.hsd$cluster)
tukeys <- lapply(colnames(nms30.k4.hsd[, -31]), function(col) {
  # Doesn't work the oneway way for some reason!
  fm <- eval(substitute(i ~ cluster, list(i = as.name(col))))
  TukeyHSD(aov(fm, nms30.k4.hsd))
  # tryCatch(TukeyHSD(aov(fm, nms30.k4.hsd)), error = function(e) print(col))
})
names(tukeys) <- colnames(nms30.k4[, -31])

for (v in names(tukeys)) {
  test <- tukeys[[v]]$cluster
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
