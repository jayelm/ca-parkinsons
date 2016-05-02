# Cleaning up and presenting data
# Because I need to learn Sweave at some poin
library(doBy)
library(xtable)
library(reshape)

nmsd.symptoms <- c(
  NMS.D.NAMES,
  "Tremor",
  "Bradykinesia",
  "Rigidity",
  "Axial",
  'Cluster'
)

nmsd.extra <- c(
  'Age',
  'Sex',
  'PD_Onset',
  'PD_Duration',
  'CISI_Total',
  'ldopa',
  'Surgery',
  'Cluster'
)

# Depends on valid cluster assignment, assert that the distribution is
# table(clus4.wide$cluster)
#   1   2   3  4
# 406 189 221 88

# ON NONMOTOR DOMAINS ====
# Somewhat oddly, the correct clustering vector lies in trees$clusters4$clustering$cluster
present <- rename(raw.omitted, PUB.MAP)
present.full <- rename(raw.omitted.full, PUB.MAP)
present$Cluster <- trees$clusters4$clustering$cluster
present.full$Cluster <- trees$clusters4$clustering$cluster

# Funcs for latex
mean.sd <- function(data, sig = 2) {
  paste(round(mean(data), sig), " (", round(sd(data), sig), ")", sep = "")
}
to.latex <- function(df, file = NULL) {
  summary.t <- t(summaryBy(. ~ Cluster, df, FUN = function(x) mean.sd(x, sig = 1), # Only 1 decimal
                           keep.names = TRUE))
  # Get rid of "cluster" row
  summary.t <- summary.t[-which(rownames(summary.t) == 'Cluster'), ]
  xt <- xtable(summary.t,
               sanitize.colnames.function = function(x) gsub(pattern = '\\_', replacement = '/', x = x))
  print(xt, type = "latex", file = file, booktabs = TRUE)
}

# Redo ANOVA + bonferroni correction.
# Apparently Tukey takes care of multiple comparisons but make sure that's not a
# setting or grouping you need to actually make happen in R.

to.latex(present[, nmsd.symptoms], "../writeup/manuscript/include/nmsd_summaries.tex")
to.latex(present.full[, nmsd.extra], "../writeup/manuscript/include/nmsd_extra.tex")
