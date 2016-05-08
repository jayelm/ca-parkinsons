# Cleaning up and presenting data
# Because I need to learn Sweave at some poin
library(doBy)
library(xtable)
library(reshape)
library(plyr)
library(infotheo)

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

# V-measure ====
# always useful...
v.measure <- function(a, b) {
  mi <- mutinformation(a, b)
  entropy.a <- entropy(a)
  entropy.b <- entropy(b)
  if (entropy.a == 0.0) {
    homogeneity <- 1.0
  } else {
    homogeneity <- mi / entropy.a
  }
  if (entropy.b == 0.0) {
    completeness <- 1.0
  } else {
    completeness <- mi / entropy.b
  }
  if (homogeneity + completeness == 0.0) {
    v.measure.score <- 0.0
  } else {
    v.measure.score <- (2.0 * homogeneity * completeness
                        / (homogeneity + completeness))
  }
  # Can also return homogeneity and completeness if wanted
  c(homogeneity, completeness, v.measure.score)
}

# ON NONMOTOR DOMAINS ====
# Somewhat oddly, the correct clustering vector lies in trees$clusters4$clustering$cluster
present <- reshape::rename(raw.omitted, PUB.MAP)
present.full <- reshape::rename(raw.omitted.full, PUB.MAP)
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

# Publication-ready dendrogram ====
remove.domain.n <- function(s) {
  splitted <- strsplit(s, '-')[[1]]
  if (length(splitted) == 1) {
    s
  } else {
    splitted[3]
  }
}
rid.of.middle <- function(s) {
  splitted <- strsplit(s, '-')[[1]]
  if (length(splitted) == 1) {
    s  # No -, leave it alone
  } else {
    # Get rid of that middle one
    # Could get rid of preceding d as well
    domain <- splitted[1]
    symp <- splitted[3]
    paste(domain, symp, sep = '-')
  }
}

labels.wo.d <- unname(sapply(labels(hm.m$colDendrogram), rid.of.middle))
# TODO: Capitalize map if necessary
par(mar=c(3, 4, 0.5, 0))
hm.m$colDendrogram %>%
  set("labels", labels.wo.d) %>%
  sort(type = "nodes") %>%
  set("branches_lwd", 2) %>%
  hang.dendrogram(hang_height = 3) %>%
  set("branches_k_color", k = 5) %>%
  # Here, "1" is motor, "2" is nonmotor (sorting by nodes is convenient here)
  color_labels(col = c(rep("blue", 4), rep("black", 30))) %>%
  plot(ylim = c(10, 45))#, xlab = "Symptom", ylab = "Height")

dev.copy(pdf, "../figures/nms30m-colhc-pub.pdf", width = 15, height = 7)
dev.off()

# Pub-ready boxplot graph ====
dev.off()
clus <- clusters.raw.long[["4"]]
# Get rid of extra
clus.pub <- clus
clus.pub <- clus.pub[clus.pub$variable != "sex", ]
clus.pub$variable <- sapply(clus.pub$variable, function(s) paste("\n", PUB.MAP.N[as.character(s)][[1]], "\n", sep = ""))
clus.pub$variable <- factor(clus.pub$variable)
clus.pub$variable <- factor(clus.pub$variable, levels(clus.pub$variable)[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 16, 17, 13, 10, 15, 14)])
# Add types. Need to do after factors have been reorganized for some reason
clus.pub$Type <- ""
clus.pub[clus.pub$variable %in%
           sapply(NMS.D, function(s) paste("\n", PUB.MAP.N[as.character(s)][[1]], "\n", sep = "")), ]$Type <- "Nonmotor (analyzed)"
clus.pub[clus.pub$variable %in% factor(c("\nAxial\n", "\nRigidity\n", "\nBradykinesia\n", "\nTremor\n")), ]$Type <- "Motor (analyzed)"
clus.pub[!(clus.pub$Type %in% c("Nonmotor (analyzed)", "Motor (analyzed)")), ]$Type <- "Other"
p <- ggplot(clus.pub, aes(x = factor(cluster), y = measurement, fill = factor(cluster))) +
  geom_boxplot() +
  guides(fill = FALSE) +
  facet_wrap( ~ variable, scales = "free") +
  xlab("") +
  ylab("") +
  theme_pub() +
  theme(strip.background = element_blank(), strip.text = element_text(lineheight = 0.4))
print(p)

dummy <- ggplot(clus.pub, aes(x = factor(cluster), y = measurement)) +
  facet_wrap( ~ variable, scales = "free") +
  geom_rect(aes(fill = Type), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_minimal() +
  theme(strip.text = element_text(lineheight = 0.4, size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  labs(fill = "Variable Type") +
  scale_fill_manual(values = c("#d7f6ff", "#d2c5ff", "#fed4f4")) +
  theme(legend.position = c(0.5, 0.1))
dummy

# Terribly complicated way to add colors
# http://stackoverflow.com/questions/19440069/ggplot2-facet-wrap-strip-color-based-on-variable-in-data-set
library(gtable)

g1 <- ggplotGrob(p)
g2 <- ggplotGrob(dummy)

gtable_select <- function (x, ...) 
{
  matches <- c(...)
  x$layout <- x$layout[matches, , drop = FALSE]
  x$grobs <- x$grobs[matches]
  x
}

panels <- grepl(pattern="panel", g2$layout$name)
strips <- grepl(pattern="strip_t", g2$layout$name)
legends <- grepl(pattern="guide-box", g2$layout$name)
g2$layout$t[panels] <- g2$layout$t[panels] - 1
g2$layout$b[panels] <- g2$layout$b[panels] - 1

new_strips <- gtable_select(g2, panels | strips | legends)
grid.newpage()
grid.draw(new_strips)

gtable_stack <- function(g1, g2){
  g1$grobs <- c(g1$grobs, g2$grobs)
  g1$layout <- transform(g1$layout, z= z-max(z), name="g2")
  g1$layout <- rbind(g1$layout, g2$layout)
  g1
}
## ideally you'd remove the old strips, for now they're just covered
new_plot <- gtable_stack(g1, new_strips)
grid.newpage()
grid.draw(new_plot)

dev.copy(pdf, "../figures/kmeans-summaries-4-pub.pdf", width = 14, height = 10)
dev.off()

# Anova with bonferroni correction on c1====
clus4.wide.st <- clusters.raw[["4"]]
# Why isn't this already a factor? Really confused
# Because it wasn't set in clusters.raw - if this is a bad thing lmk
clus4.wide.st$cluster <- as.factor(clus4.wide.st$cluster)
# Just NMS
clus4.wide.st <- clus4.wide.st[, c(NMS.D, "axial", "rigidity", "bradykin", "tremor", "cluster")]
# Assuming 1st column is cluster (which it should be)
oneways <- lapply(colnames(clus4.wide.st[, -which(colnames(clus4.wide.st) %in% c("cluster"))]), function(col) {
  fm <- substitute(i ~ cluster, list(i = as.name(col)))
  oneway.test(fm, clus4.wide.st)
})
for (test in oneways) {
  if (test$p.value < (0.05 / length(oneways))) { # BONFERRONI CORRECTION!
    cat('sig\n')
  } else {
    cat('INSIG:\n')
    cat(test$data.name, '\n')
  }
}

# Redo tukey's for sanity
tukeys <- lapply(colnames(clus4.wide.st[, -which(colnames(clus4.wide.st) %in% c("cluster"))]), function(col) {
  # Doesn't work the oneway way for some reason!
  fm <- as.formula(paste(col, '~ cluster'))
  TukeyHSD(aov(fm, clus4.wide.st))$cluster
})
names(tukeys) <- colnames(clus4.wide.st[, -which(colnames(clus4.wide.st) %in% c("cluster"))])

# Now output the insignificant diferences only
. <- sapply(names(tukeys), function(name) {
  tkdf <- as.data.frame(tukeys[name][[1]])
  sigs <- rownames(tkdf[tkdf[["p adj"]] >= 0.05, ])
  if (length(sigs) > 0) {
    cat(name, ": ", sep = "")
    cat(sigs, "\n")
  }
})

# Redo for age/sex/pdonset/duratpd/cisitot on c1 ====
clus4.wide.st <- clusters.raw[["4"]]
# Why isn't this already a factor? Really confused
# Because it wasn't set in clusters.raw - if this is a bad thing lmk
clus4.wide.st$cluster <- as.factor(clus4.wide.st$cluster)
clus4.wide.st <- clus4.wide.st[, c("age", "sex", "pdonset", "durat_pd", "cisitot", "cluster")]
# Assuming 1st column is cluster (which it should be)
oneways <- lapply(colnames(clus4.wide.st[, -which(colnames(clus4.wide.st) %in% c("cluster"))]), function(col) {
  fm <- substitute(i ~ cluster, list(i = as.name(col)))
  oneway.test(fm, clus4.wide.st)
})
for (test in oneways) {
  if (test$p.value < (0.05 / length(oneways))) { # BONFERRONI CORRECTION!
    cat('sig\n')
  } else {
    cat('INSIG:\n')
    cat(test$data.name, '\n')
  }
}
tukeys <- lapply(colnames(clus4.wide.st[, -which(colnames(clus4.wide.st) %in% c("cluster"))]), function(col) {
  # Doesn't work the oneway way for some reason!
  fm <- as.formula(paste(col, '~ cluster'))
  TukeyHSD(aov(fm, clus4.wide.st))$cluster
})
names(tukeys) <- colnames(clus4.wide.st[, -which(colnames(clus4.wide.st) %in% c("cluster"))])
. <- sapply(names(tukeys), function(name) {
  tkdf <- as.data.frame(tukeys[name][[1]])
  sigs <- rownames(tkdf[tkdf[["p adj"]] >= 0.05, ])
  if (length(sigs) > 0) {
    cat(name, ": ", sep = "")
    cat(sigs, "\n")
  }
})

# From nms30 ====
# Assert that we have
# c(509, 97, 249, 49)
nms30.present <- nms30.k4.hsd
# nms30.k4.hsd$cluster <- factor(nms30.k4.hsd$cluster)
# nms30.k4.hsd$cluster <- revalue(nms30.k4.hsd$cluster, c('3'='1', '2'='2', '1'='3', '4'='4'))
# nms30.k4.hsd$cluster <- factor(nms30.k4.hsd$cluster, levels = c('1', '2', '3', '4'))
nms30.present <- reshape::rename(nms30.present, PUB.MAP)
nms30.present <- reshape::rename(nms30.present, NMS.30.LONG.SHORT.MAP)

nms30.extra.cols <- c("Age", "Sex", "PD_Onset", "PD_Duration", "CISI_Total", "Cluster")
to.latex(nms30.present[, c(NMS.30.NAMES.PUB, MOTOR.PUB, "Cluster")],
         "../writeup/manuscript/include/nms30_summaries.tex")
to.latex(nms30.present[, nms30.extra.cols],
         "../writeup/manuscript/include/nms30_extra.tex")

# nms30 same drill, anova + tukey ====
# NOTE: cluster is captalized here since I'm using the PUB df
oneways <- lapply(colnames(nms30.present[, -which(colnames(nms30.present) %in% c("Cluster"))]), function(col) {
  fm <- substitute(i ~ Cluster, list(i = as.name(col)))
  oneway.test(fm, nms30.present)
})
for (test in oneways) {
  if (test$p.value < (0.05 / length(oneways))) { # BONFERRONI CORRECTION!
    cat('sig\n')
  } else {
    cat('INSIG:\n')
    cat(test$data.name, '\n')
  }
}

# Redo tukey's for sanity
tukeys <- lapply(colnames(nms30.present[, -which(colnames(nms30.present) %in% c("Cluster"))]), function(col) {
  # Doesn't work the oneway way for some reason!
  fm <- as.formula(paste(col, '~ Cluster'))
  TukeyHSD(aov(fm, nms30.present))$Cluster
})
names(tukeys) <- colnames(nms30.present[, -which(colnames(nms30.present) %in% c("Cluster"))])

# Now output the insignificant diferences only
. <- sapply(names(tukeys), function(name) {
  tkdf <- as.data.frame(tukeys[name][[1]])
  sigs <- rownames(tkdf[tkdf[["p adj"]] >= 0.05, ])
  if (length(sigs) > 0) {
    cat(name, ": ", sep = "")
    cat(sigs, "\n")
  }
})

# Gender binomial tests ====
print.proportions <- function(mat) {
  cat("Sex (\\% Male) ")
  sapply(1:dim(mat)[1], function(i) {
    v <- mat[i, ]
    cat("& ", round(v[1] / (v[1] + v[2]), 2) * 100, " ", sep = "")
  })
  cat("\\\\\n")
}
combs.1to4 <- combn(1:4, 2)

# For nmsd
present.sex <- table(present[c("Cluster", "Sex")])
print.proportions(present.sex)
chisq.test(present.sex)
# Welp, pairwise prop test is a much easier way to do this
# pairwise.prop.test(present.sex, p.adjust = "bonferroni")
. <- apply(combs.1to4, MARGIN = 2, FUN = function(comb) {
  pt <- prop.test(nms30.sex[comb, ])
  if (pt$p.value < (0.05 / dim(combs.1to4)[2])) { # Bonferroni correction
    cat("SIG:\n")
    cat("Cluster ", comb[1], " and ", comb[2], "\n", sep = "")
    print(pt)
  }
})

# For nms30
nms30.sex <- table(nms30.present[c("Cluster", "Sex")])
print.proportions(nms30.sex)
chisq.test(nms30.sex)
# pairwise.prop.test(nms30.sex, p.adjust = "bonferroni")
. <- apply(combs.1to4, MARGIN = 2, FUN = function(comb) {
  pt <- prop.test(nms30.sex[comb, ])
  if (pt$p.value < (0.05 / dim(combs.1to4)[2])) { # Bonferroni correction
    cat("SIG:\n")
    cat("Cluster ", comb[1], " and ", comb[2], "\n", sep = "")
    print(pt)
  }
})

# Correct nms30 heatmap ====
hm.nms30.raw.scaled <- nms30.present
# Gender no need
hm.nms30.raw.scaled$Sex <- NULL
# Nullify cluster then reattach once you've scaled
hm.nms30.raw.scaled$Cluster <- NULL
hm.nms30.raw.scaled <- as.data.frame(scale(hm.nms30.raw.scaled))
hm.nms30.raw.scaled$Cluster <- nms30.present$Cluster
hm.nms30.data <- summaryBy(. ~ Cluster, hm.nms30.raw.scaled, keep.names = T)
hm.nms30.data$Cluster <- NULL
# Re-add the domain number to the first 30
names(hm.nms30.data)[1:30] <- sapply(NMS.30.NAMES, rid.of.middle)
hm.nms30.data.t <- as.data.frame(t(hm.nms30.data))
# Reorder
hm.nms30.data.t <- hm.nms30.data.t[rownames(hm.nms30.data.t)[c(1:30, 35:38, 31:34)], ]
plot.new()
heatmap.2(as.matrix(hm.nms30.data.t), Rowv = FALSE, Colv = FALSE, dendrogram = 'none', trace = 'none',
          # cellnote = as.matrix(hm.nms30.data.t),
          col = colorRampPalette(c('green', 'black', 'red'))(n = 1000),
          # RowSideColors = c(rep(gch[1], 2), rep(gch[2], 4), rep(gch[3], 6), rep(gch[4], 3), rep(gch[5], 3),
          #                   rep(gch[6], 3), rep(gch[7], 3), rep(gch[8], 2), rep(gch[9], 4), rep(gch[10], 4)),
          xlab = 'Cluster', key.xlab = 'z-score',
          cexCol = 1.5, cexRow = 1.2, srtCol = 0,
          margins = c(5, 18),
          # Draw lines to separate categories
          rowsep = c(2, 6, 12, 15, 18, 21, 24, 26, 30, 34),
          sepcolor = "#cccccc",
          sepwidth = c(0.1, 0.1),
          lmat = rbind(c(0,3),c(2,1),c(0,4)),
          lwid = c(0.3,2),
          lhei = c(0.1,4,1),
          keysize = 0.5,
          key.par = list(mar = c(7, 8, 3, 12)),
          density.info = 'none'
          )
if (TRUE) {
  # Always write, for now
  # NOTE: I crop this in preview afterwards because it still has some
  # dead space
  dev.copy(pdf, "../figures/nms30-hm-pub.pdf", width = 7, height = 10)
  dev.off()
}

# Compute homogeneity/completeness/v-measure ====
v.measure(nms30.present$Cluster, present$Cluster)
# It's actually quite poor!
c(0.3988046,0.4661785,0.4298677)