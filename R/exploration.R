# NOTE: This is just bits and pieces of sandbox code.
# Not really meant to be sourced.

# Libraries ====
library(corrplot)
library(ggplot2)

SAVE.COR.PARS = FALSE

# Need clus4.wide from kmeans.dtree (the one with the class modified)

# Random exploration ====
c1 <- clus4.wide[clus4.wide$class == 1, ]
c2 <- clus4.wide[clus4.wide$class == 2, ]
c3 <- clus4.wide[clus4.wide$class == 3, ]
c4 <- clus4.wide[clus4.wide$class == 4, ]

# OR alternative cluster usage if the full kmeans-dtree script hasn't been run
if (FALSE) {
  c1 <- clus4.wide[clus4.wide$cluster == 1, ]
  c2 <- clus4.wide[clus4.wide$cluster == 2, ]
  c3 <- clus4.wide[clus4.wide$cluster == 3, ]
  c4 <- clus4.wide[clus4.wide$cluster == 4, ]
}

# Plot w/o cluster
par(mfrow=c(2, 2), oma=c(0, 0, 0, 0))
corrplot(cor(c1[, -19]), method = 'ellipse', title = 'Correlation 1',
         mar = c(1, 0, 1, 0), diag = F, type = 'lower')
corrplot(cor(c2[, -19]), method = 'ellipse', title = 'Correlation 2',
         mar = c(1, 0, 1, 0), diag = F, type = 'lower')
corrplot(cor(c3[, -19]), method = 'ellipse', title = 'Correlation 3',
         mar = c(1, 0, 1, 0), diag = F, type = 'lower')
corrplot(cor(c4[, -19]), method = 'ellipse', title = 'Correlation 4',
         mar = c(1, 0, 1, 0), diag = F, type = 'lower')
par(mfrow=c(1, 1))
if (SAVE.CORR.PLOT) {
  dev.copy(pdf, '../figures/corrplots.pdf')
  dev.off()
}

# Check out relationship between bradykinesia and cisitot -
# Correlation in c1 but not c2, c3, c4
# TODO
clus4.wide$class <- factor(clus4.wide$class, levels = c(1, 2, 3, 4))
p <- ggplot(clus4.wide, aes(x = bradykin, y = cisitot)) +
  facet_wrap( ~ class, scales = 'free') +
  geom_point() +
  stat_smooth(method = 'lm') +
  theme_bw()
print(p)

if (SAVE.COR.PARS) {
  ggsave("../figures/bradykin-v-cisitot.pdf")
}
# Same for rigidity - more severe PD has larger dynamics of
# brady/rigidity expression?
ggplot(clus4.wide, aes(x = rigidity, y = cisitot)) +
  facet_wrap( ~ class, scales = 'free') +
  geom_point() +
  stat_smooth(method = 'lm') +
  theme_bw()

if (SAVE.COR.PARS) {
  ggsave("../figures/rigidity-v-cisitot.pdf")
}

# Relationship between rigidity and nms_d6
p <- ggplot(clus4.wide, aes(x = bradykin, y = nms_d6)) +
  facet_wrap( ~ class, scales = 'free') +
  geom_point() +
  stat_smooth(method = 'lm') +
  theme_bw()
print(p)

if (SAVE.COR.PARS) {
  ggsave("../figures/bradykin-v-nms_d6.pdf")
}
# Same for rigidity - more severe PD has larger dynamics of
# brady/rigidity expression?
ggplot(clus4.wide, aes(x = rigidity, y = nms_d6)) +
  facet_wrap( ~ class, scales = 'free') +
  geom_point() +
  stat_smooth(method = 'lm') +
  theme_bw()

if (SAVE.COR.PARS) {
  ggsave("../figures/rigidity-v-nms_d6.pdf")
}

# NOTE: This hsan't been updated for new group names!
# Check chisq independence of sex and nms_d8 for c1 (severe group)
# Being male or female affects sexual function reports, in general
# (Assuming male here)

# CORRELATION TEST RESULTS ====
# (pearson moment)

print.cor <- function(ct, comp, n) {
  # ct$data.name preserves variable input
  cat(n, ": ", comp[[1]], " and ", comp[[2]], "\t", ct$conf.int, "\t",
      ct$p.value, "\n", sep="")
}

cor.comparisons <- list(
  c('bradykin', 'cisitot'),
  c('rigidity', 'cisitot'),
  c('bradykin', 'nms_d6'),
  c('rigidity', 'nms_d6')
)

test.cors <- function() {
  cat("Vars\tConfidence Interval\tPVal\n")
  for (comp in cor.comparisons) {
    cortest <- cor.test(c1[[comp[1]]], c1[[comp[2]]])
    print.cor(cortest, comp, "c1")
  }
  for (comp in cor.comparisons) {
    cortest <- cor.test(c2[[comp[1]]], c2[[comp[2]]])
    print.cor(cortest, comp, "c2")
  }
  for (comp in cor.comparisons) {
    cortest <- cor.test(c3[[comp[1]]], c3[[comp[2]]])
    print.cor(cortest, comp, "c3")
  }
  for (comp in cor.comparisons) {
    cortest <- cor.test(c4[[comp[1]]], c4[[comp[2]]])
    print.cor(cortest, comp, "c4")
  }
}

test.cors()

# Bayes nets ===
# Try learning bns within the clusters
c2.df <- as.data.frame(as.matrix(c2))
c4.df <- as.data.frame(as.matrix(c4))
c3.df <- as.data.frame(as.matrix(c3))
c1.df <- as.data.frame(as.matrix(c1))

# iamb w/o class
plot(iamb(c2.df[, -19]), main = '2')
plot(iamb(c4.df[, -19]), main = '4')
plot(iamb(c3.df[, -19]), main = '3')
plot(iamb(c1.df[, -19]), main = '1')

# hc w/o class
# Too sparse - need more data
plot(hc(c2.df[, -19]))
plot(hc(c4.df[, -19]))
plot(hc(c3.df[, -19]))
plot(hc(c1.df[, -19]))

# Markov blanket feature selection
parkinsons.chr <- parkinsons
parkinsons$class <- as.factor(as.character(parkinsons$class))
# Discretize first
pd.discrete <- discretize(parkinsons, method = "quantile")
iamb(parkinsons)
tan_cl(parkinsons, class = "class")

# More correlation (from asdm) ====
pd.cor <- cor(parkinsons[, 1:18])
corrplot(pd.cor, method = 'ellipse', type = 'lower', diag = F)
# Remember pigd is gone!
# pd.cor['pigd', 'axial']
# pd.cor['pigd', 'cisitot']
pd.cor['axial', 'cisitot']
pd.cor['pdonset', 'age']
pd.cor['bradykin', 'rigidity']

pd.cor['nms_d2', 'nms_d3']

pd <- parkinsons
pd[, 1:18] <- scale(pd[, 1:18])
pd.iamb <- iamb(pd)
plot(pd.iamb)
pd.tan_cl <- tan_cl(pd, class = "class")
pd.hc <- hc(pd)
pd.cpt <- cptable(pd)
as.grain(pd.cpt)
pd.whitelist <- data.frame(
  from = c('axial', 'tremor', 'bradykin', 'rigidity'),
  to = c('cisitot', 'cisitot', 'cisitot', 'cisitot')
)
pd.hc.wl <- hc(pd, whitelist = pd.whitelist)
pd.hc <- si.hiton.pc(pd)
plot(pd.hc)
pd.hc.mb <- mb(pd.hc, 'class')
pd.hc.mb
length(mb(pd.hc.mb, 'class'))

pd.woclass <- pd[, 1:18]
pd.woclass.hc <- hc(pd[, 1:18])
plot(pd.woclass.hc)
mb(pd.woclass.hc, 'cisitot')

# PCA ====
pc <- princomp(pd[, 1:18])
plot(pc)
