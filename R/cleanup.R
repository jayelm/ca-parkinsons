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

# Publication-ready dendrogram ====
labels.wo.d <- unname(sapply(labels(hm.m$colDendrogram), function(s) {
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
}))
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

# Pub-ready boxplot graph
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
