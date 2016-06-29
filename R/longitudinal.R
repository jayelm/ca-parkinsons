# Longitudinal analysis among clusters. Assumes we have the data
# structures ran in kmeans-dtree.R

library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)
library(RColorBrewer)
library(reshape2)
library(plyr)
library(corrplot)

library(ggthemes)
source('./pubtheme.R')

SAVE.LONG.PLOTS <- TRUE

# durat_pd distribution for all and subclus ====
summary(raw.omitted$durat_pd)
table(raw.omitted$durat_pd)

qplot(raw.omitted$durat_pd)
qplot(raw.omitted[trainset.labeled$cluster == 1, ]$durat_pd)
qplot(raw.omitted[trainset.labeled$cluster == 2, ]$durat_pd)
qplot(raw.omitted[trainset.labeled$cluster == 3, ]$durat_pd)
qplot(raw.omitted[trainset.labeled$cluster == 4, ]$durat_pd)

# Bind durat_pd to 30 nms, CLUSTER PARAM===
everything.wide <- raw.omitted
everything.wide$durat_pd <- raw.omitted$durat_pd
everything.wide$cluster <- trainset.labeled$cluster
nms.d <- c("nms_d1", "nms_d2", "nms_d3", "nms_d4", "nms_d5", "nms_d6",
           "nms_d7", "nms_d8", "nms_d9")

# Change this if you want to examine a different cluster!
# use NULL if everything, and comment out subsetting
CURR.CLUSTER <- NULL
# everything.wide <- everything.wide[everything.wide$cluster == CURR.CLUSTER, ]
summary(everything.wide)


# Either: Divide data into bins. ====

se <- function(x) sqrt(var(x)/length(x))
# Optimal - for N, seems like best bin size is N
# Code for ALL
BIN_SIZE = 20
binned <- data.frame(matrix(ncol = 0, nrow = BIN_SIZE))
breaks <- hist(everything.wide$durat_pd, BIN_SIZE, plot=F)$breaks

# Code for K - specify the breaks more specifically
# breaks <- seq(from = 0, to = max(everything.wide$durat_pd), by = 2)
# hist(everything.wide$durat_pd, breaks = breaks, plot=T)
# binned <- data.frame(matrix(ncol = 0, nrow = length(breaks) - 1))

# COMMON CODE
all.breaks <- cut(everything.wide$durat_pd, breaks)
rownames(binned) <- levels(all.breaks)

binned$counts <- as.numeric(table(all.breaks))

for (col in colnames(everything.wide)) {
  if (col %in% c("breaks", "durat_pd", "cluster")) {
    next
  }
  binned[, col] <- tapply(everything.wide[, col], all.breaks, mean)
  binned[, paste(col, "_sd", sep="")] <- tapply(everything.wide[, col], all.breaks, se)
}

# Remove rows with too few counts. We'll assume <5
binned <- binned[!(binned$counts < 5), ]

# Or: Divide data into bins, by cluster ====
# TODO: make this binned_multi?
breaks.m <- seq(from = 0, to = max(everything.wide$durat_pd), by = 2)
binned.m <- data.frame(matrix(ncol = 0, nrow = length(breaks.m) - 1))
for (k in 1:4) {
  all.breaks.k.m <- cut(everything.wide[everything.wide$cluster == k, ]$durat_pd, breaks.m)
  rownames(binned.m) <- levels(all.breaks.k.m)  # Redundant
  binned.m[[paste("counts_k", k, sep="")]] <- as.numeric(table(all.breaks.k.m))
}
binned.m$counts_total <- binned.m$counts_k1 + binned.m$counts_k2 + binned.m$counts_k3 + binned.m$counts_k4

# Keep only counts where all of the subclusters have >= 5
for (col in colnames(everything.wide)) {
  if (col %in% c("breaks", "durat_pd", "cluster")) {
    next
  }
  for (k in 1:4) {
    c_lvector <- everything.wide$cluster == k
    binned.m[, paste(col, "_k", k, sep="")] <- tapply(everything.wide[c_lvector, col], all.breaks[c_lvector], mean);
    binned.m[, paste(col, "_sd", "_k", k, sep="")] <- tapply(everything.wide[c_lvector, col], all.breaks[c_lvector], se)
  }
}

binned.m <- binned.m[!(binned.m$counts_k1 < 4 | binned.m$counts_k2 < 4 | binned.m$counts_k3 < 4 | binned.m$counts_k4 < 4), ]

# Plot MULTI obs: setup ====

# Emulate ggplot palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

# Melt to long format for multi.nms stacked bar chart ====
binned.copy.m <- binned.m
binned.copy.m$names <- rownames(binned.copy.m)
binned.wide.pbar.m <- melt(binned.copy.m[, c("counts_k1", "counts_k2", "counts_k3", "counts_k4", "names")], id.var="names")
# Rename factor, subtype
levels(binned.wide.pbar.m$variable) <- c('1', '2', '3', '4')
binned.wide.pbar.m <- rename(binned.wide.pbar.m, c('variable'='Subtype', 'value'='Counts'))
rm(binned.copy.m)
binned.wide.pbar.m$names <- factor(binned.wide.pbar.m$names, levels=rownames(binned.m))
p.bar.m <- ggplot(binned.wide.pbar.m, aes(x = factor(names), y = Counts, fill = Subtype)) + 
  geom_bar(stat = "identity") +
  labs(x = "PD Duration") +
  # guides(fill = FALSE) +  # Guide already present in plot.nms
  theme_pub() +
  theme(plot.margin = unit(c(0.5,1,1,1), "cm"))

# Melt to long format for multi.nms ====
binned.copy.m <- binned.m
binned.copy.m$names <- rownames(binned.copy.m)
# Get variables we want
measures <- names(binned.copy.m)[grep('.*sd.*|(names)', names(binned.copy.m), invert=T)]
binned.wide.multi <- melt(
  binned.copy.m,
  id.var='names',
  variable.name='Subtype',
  value.name='Score',
  # Get rid of sds!!
  measure.vars=measures
)
# Get ymin, ymax for all of these variables
binned.wide.multi$ymin <- rep(NA, nrow(binned.wide.multi))
binned.wide.multi$ymax <- rep(NA, nrow(binned.wide.multi))
for (i in grep('(nms.*sd.*)|((tremor|cisitot).*sd.*)', names(binned.copy.m))) {
  # Get rid of _sd_ to get the real name
  name <- names(binned.copy.m)[[i]]
  vals <- binned.copy.m[[name]]
  orig <- gsub('_sd', '', name)
  ymins <- binned.wide.multi[binned.wide.multi$Subtype == orig, ]$Score - vals
  ymins[ymins < 0] <- 0
  ymaxs <- binned.wide.multi[binned.wide.multi$Subtype == orig, ]$Score + vals
  binned.wide.multi[binned.wide.multi$Subtype == orig, ]$ymin <- ymins
  binned.wide.multi[binned.wide.multi$Subtype == orig, ]$ymax <- ymaxs
}
# Assume 9 categories here.
binned.wide.multi$names <- rep(1:length(unique(binned.wide.multi$names)), nrow(binned.wide.multi) / 9)
rm(binned.copy.m)

# Create multi.nms, multi.nms.with.counts funcs ====
multi.nms <- function(nms_str, save=FALSE, means=FALSE) {
  colvec <- as.character(lapply(1:4, function(n) {paste(nms_str, '_k', n, sep='')}))
  bwm.subset <- binned.wide.multi[binned.wide.multi$Subtype %in% colvec, ]
  # Rename to subtypes - assumes 4 groups
  bwm.subset$Subtype <- rep(
    as.character(lapply(1:4, as.character)),
  each=9)
  pd <- position_dodge(width=0.2)
  p <- ggplot(bwm.subset,
    aes(x = names, y=Score, group=Subtype, color=Subtype, fill=Subtype)
  ) +
    geom_point(position=pd) +
    geom_smooth(se=FALSE) +
    geom_errorbar(aes(ymin=ymin, ymax=ymax), position=pd) +
    scale_x_continuous(breaks=1:9, labels=rownames(binned.m))
  # TODO: Melt to different format like above for p.bar; then position_dodge
  # http://stats.stackexchange.com/questions/103393/how-can-i-dodge-the-position-of-geom-point-in-ggplot2
  
  # Offset. Might need to tweak it
  mean.offset <- .1
#   if (means) {
#     abline.df <- data.frame(intercept=subclus.means, slope=rep(0, 4), colour=as.factor(1:4), text.x=rep(nrow(binned) - 0.5, 4))
#     p <- p +
#       geom_abline(aes(intercept = intercept, slope=slope, colour=colour), data=abline.df, linetype='dashed') +
#       annotate("text", x=abline.df$text.x, y=abline.df$intercept + mean.offset,
#                label=paste("µ = ", round(abline.df$intercept, 2), sep=""), size=6, colour=gg_color_hue(4))
#   }
#   
  # Global mean
  mean.nms <- mean(raw.omitted[[nms_str]])
  p <- p +
    geom_abline(aes(intercept=mean.nms, slope=0), linetype='dashed') +
    annotate("text", x=nrow(binned.m) - 0.5, y=mean.nms + mean.offset, label=paste("µ = ", round(mean.nms, 2), sep=""), size=6)
  
  # Formatting
  p <- p +
    theme_pub() +
    blank.theme +
    ggtitle(nms_str) +
    ylab('Score') +
    theme(plot.margin = unit(c(1,1,0,1), "cm")) +
    labs(colour="Subtype")
  
  if (save) {
    dev.copy(pdf, paste('../figures/longitudinal/', nms_str, '-multi-durat', '.pdf', sep=''),
             width=14, height=10)
    dev.off()
  }
  p
}

multi.nms.with.counts <- function(nms_str, save=FALSE, means=FALSE) {
  nms.plot <- multi.nms(nms_str, means=means)
  gA <- ggplot_gtable(ggplot_build(nms.plot))
  gB <- ggplot_gtable(ggplot_build(p.bar.m))
  maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
  gA$widths[2:5] <- as.list(maxWidth)
  gB$widths[2:5] <- as.list(maxWidth)
  grid.newpage()
  grid.arrange(
    arrangeGrob(gA, gB, ncol=1, heights=c(.8,.3))
  )
  if (save) {
    dev.copy(pdf, paste('../figures/longitudinal/', nms_str, '-multi-durat-counts', '.pdf', sep=''),
             width=14, height=10)
    dev.off()
  }
}


# Plot ALL obs: setup ====
blank.theme <- theme(axis.text.x = element_blank(),
                     axis.ticks = element_blank(),
                     axis.title.x = element_blank())

p.bar.solo <- ggplot(binned, aes(x=factor(1:nrow(binned)), y=counts)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=counts), position=position_dodge(width=0.9), vjust=-0.25) +
  labs(x="PD Duration") +
  theme_pub() +
  theme(plot.margin = unit(c(0.5,1,1,1), "cm")) +
  scale_x_discrete(breaks=1:nrow(binned), labels = rownames(binned))

plot.nms <- function(nms_str, save=FALSE) {
  aes_nms <- aes_string(x = "factor(1:nrow(binned))", y=nms_str)
  aes_sd <- aes_string(
    ymin=paste(nms_str, "-", nms_str, "_sd", sep=""),
    ymax=paste(nms_str, "+", nms_str, "_sd", sep="")
  )
  aes_smooth <- aes_string(x = "as.numeric(factor(1:nrow(binned)))", y=nms_str)
  mean.nms <- mean(raw.omitted[[nms_str]])
  # raw.omitted essentially the same as everything.wide
  mean.subclus <- ifelse(!is.null(CURR.CLUSTER), mean(everything.wide[everything.wide$cluster == CURR.CLUSTER, nms_str]), Inf)
  print(mean.subclus)
  p <- ggplot(binned, aes_nms) +
    geom_errorbar(aes_sd, width=0.25) +
    geom_point() +
    geom_smooth(aes_smooth, method="loess", se=F) +
    # Mean line
    geom_abline(aes(intercept=mean.nms, slope=0, colour='mean'), linetype='dashed') +
    annotate("text", x=nrow(binned) - 0.5, y=mean.nms + 0.5, label=paste("µ = ", round(mean.nms, 2), sep=""), size=6, colour='red') +
    # Subclus mean line
    # IF curr.cluster null, mean.subclus Inf
    geom_abline(aes(intercept=mean.subclus, slope=0, colour='subclusmean'), linetype='dashed') +
    # This one needs to be turned off, though
    # FIXME: Not working! Not sure why 
    # TODO: One graph, stacked bar chart, multiple horizontal lines that show timeline of disease
    # progression for the various subclusters. Would be super interesting
#     ifelse(!is.null(mean.subclus),
#            annotate("text", x=nrow(binned) - 0.5, y=mean.subclus + 0.5, label=paste("µ = ", round(mean.subclus, 2), sep=""), size=6, colour='cyan'),
#            annotate("text", x=0, y=0, label="")) +
    theme_pub() +
    blank.theme +
    ggtitle(nms_str) +
    theme(plot.title=element_text(vjust=9)) +
    theme(plot.margin = unit(c(1,1,0,1), "cm")) +
    guides(colour=FALSE) +
    scale_x_discrete(breaks=1:nrow(binned), labels = rownames(binned))
  if (save) {
    dev.copy(pdf, paste('../figures/longitudinal/', nms_str, '-durat', '.pdf', sep=''),
             width=14, height=10)
    dev.off()
  }
  p
}

plot.nms.with.counts <- function(nms_str, save=FALSE) {
  nms.plot <- plot.nms(nms_str)
  gA=ggplot_gtable(ggplot_build(nms.plot))
  gB=ggplot_gtable(ggplot_build(p.bar.solo))
  maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
  gA$widths[2:5] <- as.list(maxWidth)
  gB$widths[2:5] <- as.list(maxWidth)
  grid.newpage()
  grid.arrange(
    arrangeGrob(gA, gB, ncol=1, heights=c(.8,.3))
  )
  if (save) {
    dev.copy(pdf, paste('../figures/longitudinal/', nms_str, '-durat-counts', '.pdf', sep=''),
             width=14, height=10)
    dev.off()
  }
}

plot.nms.resid <- function(nms_str, save=FALSE) {
  par(mfrow=c(1, 1))
  par(mar=rep(3, 4))
  fm <- substitute(i ~ n, list(i = as.name(nms_str)))
  nms.linmod <- lm(fm, binned)
  plot(nms.linmod$residuals, main=paste(nms_str, "linear residuals", sep=" "), xaxt='n')
  grid()
  axis(1, at=1:13, labels=rownames(binned), cex.axis=0.8)
  abline(a=0, b=0)
  if (save) {
    dev.copy(pdf, paste('../figures/longitudinal/', nms_str, '-resids', '.pdf', sep=''),
             width=14, height=10)
    dev.off()
  }
  nms.linmod
}

plot.nms.seg.resid <- function(nms_str, save=FALSE) {
  nms.linmod <- plot.nms.resid(nms_str)
  nms.segmod <- segmented(nms.linmod, seg.Z = ~ n)
  plot(nms.segmod$residuals, main=paste(nms_str, "segmented residuals", sep=" "), xaxt='n')
  grid()
  axis(1, at=1:13, labels=rownames(binned), cex.axis=0.8)
  abline(a=0, b=0)
  if (save) {
    dev.copy(pdf, paste('../figures/longitudinal/', nms_str, '-resids-seg', '.pdf', sep=''),
             width=14, height=10)
    dev.off()
  }
  nms.segmod
}

plot.nms.seg <- function(nms_str, save=FALSE) {
  nms.segmod <- plot.nms.seg.resid(nms_str)
  par(mar=rep(6, 4))
  binned$n <- 1:nrow(binned)
  plot(binned$n, binned[[nms_str]], main=nms_str,
       pch=19,
       xlab="PD Duration (years)", xaxt='n',
       ylab=nms_str)
  axis(1, at=1:13, labels=rownames(binned), cex.axis=0.8)
  sdminus.capped <- binned[[nms_str]] - binned[[paste(nms_str, "_sd", sep="")]]
  sdplus.capped <- binned[[nms_str]] + binned[[paste(nms_str, "_sd", sep="")]]
  arrows(binned$n, sdminus.capped, binned$n, sdplus.capped,
         length=0.05, angle=90, code=3)
  plot(nms.segmod, add=T)
  abline(a=mean(raw.omitted[[nms_str]]), 0, col="red", lty=2)
  abline(a=mean(raw.omitted[[nms_str]]), 0, col="red", lty=2)
  cat("Mean for ", nms_str, ":", mean(raw.omitted[[nms_str]]), "\n", sep="")
  cat("Estimated breakpoints: \n")
  print(nms.segmod$psi)
  grid()
  if (save) {
    dev.copy(pdf, paste('../figures/longitudinal/', nms_str, '-seg-all.pdf', sep=''), width=14, height=10)
    dev.off()
  }
}

# First look at with durat_pd to find interesting things ====
par(mfrow=c(1, 1))
binned.no.sd <- binned[, grep("_sd$", colnames(binned), invert=TRUE)]
binned.no.sd <- binned.no.sd[, -which(names(binned.no.sd) %in% c('counts', 'pdonset', 'age', 'sex'))]
if (!is.null(binned$n)) binned$n <- NULL
corrplot(cor(binned.no.sd), order="hclust")
if (SAVE.LONG.PLOTS) {
  dev.copy(pdf, '../figures/longitudinal/corr-binned.pdf', width=10, height=10)
  dev.off()
}

# Correlation with PD Duration BINNING ====
durat.cor <- function(arr) cor(1:nrow(binned), arr)
correlations <- sapply(binned.no.sd, durat.cor)
correlations <- correlations[which(names(correlations) != "counts")]  # Lose counts
correlations <- sort(correlations)  # Sort ascending
# Get a vector of nms_d{1-9}
# and nms{1-30} according to what they are
is_d <- grepl("d", names(correlations))
is_d[which(is_d == TRUE)] <- "nms_d{1-9}"
is_d[which(is_d == FALSE)] <- "nms{1-30}"
is_d[grep("nms.*", names(correlations), invert=TRUE)] <- "motor"
correlations.df <- data.frame(
  names=names(correlations),
  r=correlations,
  variable=is_d
)

correlations.df$names <- factor(correlations.df$names,
                                levels=names(sort(correlations)))
ggplot(correlations.df, aes(x=names, y=r, fill=variable)) +
  geom_bar(stat="identity", position="identity") +
  geom_text(aes(label=round(r, 2)), position=position_dodge(width=0.9), vjust=2 * (correlations.df$r < 0) - .5) +
  scale_y_continuous(limits = c(0, 1)) +
  ylab("r") +
  xlab("")+ 
  guides(guides(fill=guide_legend(title="Variable Type"))) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Correlation with PD duration\n")
# Only save with different names!
if (SAVE.LONG.PLOTS) {
  ggsave('../figures/longitudinal/pd-durat-cor.pdf', width=14, height=7)
}

# Durat cor correlation NO BINNING ====
everything.wide <- cbind(raw.omitted, cluster = cl$cluster)
durat.cor.everything <- function(symptom) cor(everything.wide$durat_pd, as.numeric(everything.wide[[symptom]]))
durat.cor.test <- function(symptom) cor.test(everything.wide$durat_pd, everything.wide[[symptom]])
correlations.everything <- sapply(names(everything.wide),
                                  durat.cor.everything)
# Get rid of cluster, sex, durat_pd
# If include PD_onset, then need to change up scale
to.remove <- c("cluster", "sex", "durat_pd")
correlations.everything <- correlations.everything[!names(correlations.everything) %in% to.remove]
names(correlations.everything) <- sapply(names(correlations.everything), function(v) c(NMS.D.MAP.PUB.N, NMS.NUM.TO.PUB, MISC.MAP)[[v]])
correlations.everything <- sort(correlations.everything)  # Sort ascending
# Pretty meaningless - no negative correlations!!
is_d.e <- grepl("d", names(correlations.everything))
is_d.e[which(is_d.e == TRUE)] <- "Nonmotor symptom"
is_d.e[which(is_d.e == FALSE)] <- "Nonmotor domain"
is_d.e[which(names(correlations.everything) %in% MOTOR.PUB)] <- "Motor symptom"
is_d.e[which(names(correlations.everything) %in% c("CISI_PD_total", "Age", "PD_onset"))] <- "Other"
is_d.e <- factor(is_d.e, levels = c("Nonmotor domain", "Nonmotor symptom", "Motor symptom", "Other"))
is_d.cols <- c(brewer.pal(8, "Set2")[5], brewer.pal(8, "Pastel2")[5], brewer.pal(8, "Set2")[c(6, 7)])

correlations.df.e <- data.frame(
  names=names(correlations.everything),
  r=correlations.everything,
  variable=is_d.e
)

correlations.df.e$names <- factor(correlations.df.e$names,
                                  levels=names(sort(correlations.everything)))
# correlations.test.e <- lapply(c(ALL.SYMPTOMS, NMS.30), durat.cor.test)
# names(correlations.test.e) <- c(ALL.SYMPTOMS, NMS.30)
ggplot(correlations.df.e, aes(x=names, y=r, fill=variable)) +
  geom_bar(stat="identity", position="identity") +
  # geom_text(aes(label=round(r, 2)), position=position_dodge(width=0.9), vjust=2 * (correlations.df.e$r < 0) - .5) +
  scale_y_continuous(limits = c(-1, 1)) +
  ylab("r\n") +
  xlab("Variable") +
  guides(guides(fill=guide_legend(title="Variable Type"))) +
  scale_fill_manual(values = is_d.cols) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Correlation with PD duration\n")
# Only save with different names!
if (TRUE) {
  ggsave('../figures/cor-unbinned.pdf', width=16, height=8)
}


# SEGMENTED regression ====
library(segmented)

# nms10 depression ====
# Add n
binned$n <- 1:nrow(binned)
linear.mod <- lm(nms10 ~ n, binned)
segmented.mod <- segmented(linear.mod, seg.Z = ~ n)

# Plot it (more specific for depression) ====
par(mar=rep(6, 4))
binned$n <- 1:nrow(binned)
plot(binned$n, binned$nms10, main="Segmented regression: nms10 (depression)",
     pch=19,
     xlab="PD Duration (years)", xaxt='n',
     ylab="NMS10 Severity x Frequency (Max = 12)", ylim=c(0, 5))
axis(1, at=1:13, labels=rownames(binned), cex.axis=0.8)
sdminus.capped <- binned$nms10 - binned$nms10_sd
sdplus.capped <- binned$nms10 + binned$nms10_sd
arrows(binned$n, sdminus.capped, binned$n, sdplus.capped,
       length=0.05, angle=90, code=3)
plot(segmented.mod, add=T)
abline(a=mean(raw.omitted$nms10), 0, col="red", lty=2)
grid()
text(12.3, 2.38, labels="µ = 2.17", col="red")
if (SAVE.LONG.PLOTS) {
  dev.copy(pdf, '../figures/longitudinal/nms10-seg-all-fancy.pdf', width=14, height=10)
  dev.off()
}
# Display residuals ====
par(mfrow=c(2, 1), mar=rep(3, 4))
plot(linear.mod$residuals, ylim=c(-.75, .75), main="Linear residuals", xaxt='n')
grid()
axis(1, at=1:13, labels=rownames(binned), cex.axis=0.8)
abline(a=0, b=0)
plot(segmented.mod$residuals, ylim=c(-.75, .75), main="Segmented residuals", xaxt='n')
abline(a=0, b=0)
grid()
axis(1, at=1:13, labels=rownames(binned), cex.axis=0.8)
if (SAVE.LONG.PLOTS) {
  dev.copy(pdf, '../figures/longitudinal/nms10-resids.pdf', width=14, height=10)
  dev.off()
}

# Multiplot func ====
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Stop binning. Just plot normal stuff ====
# ggplot(everything.wide, aes(x=durat_pd))
facts.of.int <- c("\nAnxiety\n" = "nms9", "\nDepression\n" = "nms10",
                  "\nCISI Total\n" = "cisitot", "\nTremor\n" = "tremor")

el <- melt(everything.wide, id.var = c("cluster", "durat_pd"))
# Filter only less than 30s, for lack of
el <- el[el$durat_pd < 30, ]
el.sub <- el[el$variable %in% facts.of.int, ]
# No better way to switch the names?
el.sub$variable <- revalue(factor(el.sub$variable), setNames(names(facts.of.int), facts.of.int))
  
mean_vals <- data.frame(
  variable = names(facts.of.int),
  value = sapply(names(facts.of.int), function(f) {
    mean(el.sub[el.sub$variable == f, ]$value)
  })
)
ggplot(el.sub, aes(x=durat_pd, y=value, color = cluster)) +
  geom_point() +
  geom_hline(aes(yintercept = value), mean_vals, linetype='dashed') +
  geom_smooth(aes(color = "Overall"), se = FALSE, color = "black") +
  geom_smooth(se = FALSE) +
  geom_jitter(width = 0.7, height = 0.7) +
  facet_wrap(~ variable, nrow = 2, ncol = 2, scales = "free_y") +
  theme_bw() +
  ylab("Symptom Score\n") +
  xlab("\nPD Duration") +
  labs(color = "Cluster") +
  theme_pub() +
  theme(strip.text = element_text(lineheight = 0.5))