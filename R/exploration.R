# Libraries ====
library(corrplot)
library(ggplot2)

# Need clus4.wide from kmeans.dtree (the one with the class modified)

# Random exploration ====
c1 <- clus4.wide[clus4.wide$class == 1, ]
c2 <- clus4.wide[clus4.wide$class == 2, ]
c3 <- clus4.wide[clus4.wide$class == 3, ]
c4 <- clus4.wide[clus4.wide$class == 4, ]

# Plot w/o cluster
par(mfrow=c(2, 2))
corrplot(cor(c2[, -20]), method = 'pie', title = 'Correlation 2', diag = F, type = 'lower')
corrplot(cor(c4[, -20]), method = 'pie', title = 'Correlation 4', diag = F, type = 'lower')
corrplot(cor(c3[, -20]), method = 'pie', title = 'Correlation 3', diag = F, type = 'lower')
corrplot(cor(c1[, -20]), method = 'pie', title = 'Correlation 1', diag = F, type = 'lower')
par(mfrow=c(1, 1))

par(mfrow=c(2, 2))

# Check out relationship between bradykinesia and cisitot -
# Correlation in c1 but not c2, c3, c4
clus4.wide$class <- factor(clus4.wide$class, levels = c(2, 4, 3, 1))
ggplot(clus4.wide, aes(x = bradykin, y = cisitot)) +
  facet_wrap( ~ class, scales = 'free') +
  geom_point() +
  stat_smooth(method = 'lm') +
  theme_bw()

# Same for rigidity - more severe PD has larger dynamics of
# brady/rigidity expression?
ggplot(clus4.wide, aes(x = rigidity, y = cisitot)) +
  facet_wrap( ~ class, scales = 'free') +
  geom_point() +
  stat_smooth(method = 'lm') +
  theme_bw()

# Check chisq independence of sex and nms_d8 for c1 (severe group)
# Being male or female affects sexual function reports, in general
# (Assuming male here)
chisq.test(c2$rigidity, c2$cisitot, simulate.p.value = T)
chisq.test(c4$pdonset, c4$age, simulate.p.value = T)
cor.test(c4$bradykin, c4$durat_pd)


# Bayes nets ===
# Try learning bns within the clusters
c2.df <- as.data.frame(as.matrix(c2))
c4.df <- as.data.frame(as.matrix(c4))
c3.df <- as.data.frame(as.matrix(c3))
c1.df <- as.data.frame(as.matrix(c1))

# iamb w/o class
plot(iamb(c2.df[, -20]), main = '2')
plot(iamb(c4.df[, -20]), main = '4')
plot(iamb(c3.df[, -20]), main = '3')
plot(iamb(c1.df[, -20]), main = '1')

# hc w/o class
# Too sparse - need more data
plot(hc(c2.df[, -20]))
plot(hc(c4.df[, -20]))
plot(hc(c3.df[, -20]))
plot(hc(c1.df[, -20]))
